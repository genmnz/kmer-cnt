// kc-c7: a faster k-mer counter derived from kc-c4.
//
// Changes relative to kc-c4 (see README for measured speedups):
//   1. A dedicated producer thread decompresses the input with igzip (ISA-L,
//      SIMD-accelerated) *and* parses reads, concurrently with counting.  In
//      kc-c4 decompression (zlib) and parsing sit in a pipeline stage that
//      competes with counting; here they overlap it almost entirely.
//   2. K-mer *extraction* is parallelised.  In kc-c4 extraction runs in a
//      single pipeline thread while only insertion is multi-threaded; here both
//      extraction and insertion run across all worker threads.  Each worker
//      extracts into its own per-prefix buffers, then the per-prefix hash
//      tables are filled in parallel (one prefix per worker, so no locking).
//   3. The hash tables are pre-sized from the (known) uncompressed length so
//      insertion never pays for incremental rehash-resizes, and insertion
//      software-prefetches the hash bucket for a look-ahead k-mer to hide the
//      DRAM latency of the random-access probes.
//
// Input is read whole, then decompressed/parsed on the producer thread into
// zero-copy per-read spans (single-line FASTA/FASTQ, matching kseq for
// sequencing reads) while the main thread dispatches blocks of reads to the
// worker pool for counting.
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include "ketopt.h"
#include "kthread.h"
#include "khashl.h"
#include <isa-l.h>

#define KC_BITS 10
#define KC_MAX ((1<<KC_BITS) - 1)
#define kc_c4_eq(a, b) ((a)>>KC_BITS == (b)>>KC_BITS)
#define kc_c4_hash(a) ((a)>>KC_BITS)
KHASHL_SET_INIT(, kc_c4_t, kc_c4, uint64_t, kc_c4_hash, kc_c4_eq)

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123, everything else to 4
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// bijective 64-bit mix (MurmurHash3 finalizer): distinct k-mers stay distinct,
// so counts are unaffected, but bucket/table distribution is well balanced.
static inline uint64_t mix64(uint64_t x)
{
	x ^= x >> 33;
	x *= 0xff51afd7ed558ccdULL;
	x ^= x >> 33;
	x *= 0xc4ceb9fe1a85ec53ULL;
	x ^= x >> 33;
	return x;
}

typedef struct {
	int p;        // prefix length: 1<<p hash tables
	kc_c4_t **h;  // 1<<p hash tables
} kc_c4x_t;

static kc_c4x_t *c4x_init(int p)
{
	int i;
	kc_c4x_t *h;
	CALLOC(h, 1);
	MALLOC(h->h, 1<<p);
	h->p = p;
	for (i = 0; i < 1<<p; ++i)
		h->h[i] = kc_c4_init();
	return h;
}

typedef struct { // a growable buffer of hashed k-mers
	int64_t n, m;
	uint64_t *a;
} buf_t;

// ---- parsed input ----
typedef struct {
	const char *s; // pointer into the decompressed buffer
	int32_t l;     // length of the read
} read_t;

// ---- per-block parallel context ----
typedef struct {
	int k, p, n_thread;
	read_t *reads;      // all reads
	int64_t r0, r1;     // [r0, r1) reads in the current block
	int n_chunk;        // number of read chunks in the block
	buf_t *buf;         // n_thread * (1<<p) buffers, row-major [tid*P + prefix]
	kc_c4x_t *h;
} blk_t;

static inline void buf_push(buf_t *b, uint64_t y)
{
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
	}
	b->a[b->n++] = y;
}

static void extract_worker(void *data, long chunk, int tid) // parallel k-mer extraction
{
	blk_t *bk = (blk_t*)data;
	int k = bk->k, p = bk->p;
	uint64_t mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2, pmask = (1ULL<<p) - 1;
	buf_t *tbuf = &bk->buf[(int64_t)tid << p]; // this worker's per-prefix buffers
	int64_t nr = bk->r1 - bk->r0;
	int64_t c0 = bk->r0 + nr * chunk / bk->n_chunk;
	int64_t c1 = bk->r0 + nr * (chunk + 1) / bk->n_chunk;
	int64_t r;
	for (r = c0; r < c1; ++r) {
		const char *seq = bk->reads[r].s;
		int len = bk->reads[r].l, i, l;
		uint64_t x0 = 0, x1 = 0;
		for (i = l = 0; i < len; ++i) {
			int c = seq_nt4_table[(uint8_t)seq[i]];
			if (c < 4) {
				x0 = (x0 << 2 | c) & mask;
				x1 = x1 >> 2 | (uint64_t)(3 - c) << shift;
				if (++l >= k) {
					uint64_t y = x0 < x1? x0 : x1;
					uint64_t h = mix64(y);
					buf_push(&tbuf[h & pmask], h);
				}
			} else l = 0, x0 = x1 = 0;
		}
	}
}

#define KC_PREFETCH 12 // how far ahead to prefetch hash buckets

static void insert_worker(void *data, long prefix, int tid) // parallel insertion, one prefix per call
{
	blk_t *bk = (blk_t*)data;
	int p = bk->p, t, nt = bk->n_thread;
	kc_c4_t *g = bk->h->h[prefix];
	for (t = 0; t < nt; ++t) {
		buf_t *b = &bk->buf[((int64_t)t << p) + prefix];
		int64_t j, n = b->n;
		const uint64_t *a = b->a;
		for (j = 0; j < n; ++j) {
			khint_t itr;
			int absent;
			// prefetch the bucket for a future k-mer to hide DRAM latency;
			// bits is stable because tables are pre-sized (no resize here).
			int64_t jp = j + KC_PREFETCH;
			if (jp < n) {
				khint_t hb = (khint_t)(a[jp]>>p) * 2654435769U >> (32 - g->bits);
				__builtin_prefetch(&g->keys[hb], 1, 1);
			}
			itr = kc_c4_put(g, a[j]>>p<<KC_BITS, &absent);
			if ((kh_key(g, itr)&KC_MAX) < KC_MAX) ++kh_key(g, itr);
		}
	}
}

// ---- decompression ----
static uint8_t *slurp(const char *fn, size_t *len)
{
	FILE *fp = strcmp(fn, "-")? fopen(fn, "rb") : stdin;
	uint8_t *buf = 0;
	size_t cap = 0, n = 0;
	if (fp == 0) return 0;
	cap = 1<<24; buf = (uint8_t*)malloc(cap);
	for (;;) {
		size_t r;
		if (n == cap) { cap <<= 1; buf = (uint8_t*)realloc(buf, cap); }
		r = fread(buf + n, 1, cap - n, fp);
		n += r;
		if (r == 0) break;
	}
	if (fp != stdin) fclose(fp);
	*len = n;
	return buf;
}

// ---- streaming decompression (igzip): decompress on a producer thread while
// the main thread parses and counts already-available bytes.
// Resumable parser: parse whole records from out[pos, limit) into the
// pre-allocated $reads array (capacity $rcap, never exceeded because every
// stored read consumes >= k+2 bytes).  A record is only committed when fully
// present; if $final is 0 an incomplete trailing record is left for the next
// call.  Single-line FASTA/FASTQ (matches kseq for sequencing reads).
static inline size_t find_nl(const char *o, size_t from, size_t lim)
{
	const char *q = (const char*)memchr(o + from, '\n', lim - from);
	return q? (size_t)(q - o) : lim;
}

typedef struct {
	uint8_t *comp; size_t clen; // input (gzip or plain)
	int is_gzip, k;
	char *out; size_t out_cap;   // output buffer (pre-allocated)
	read_t *reads; int64_t rcap; // pre-allocated read array
	size_t avail;                // bytes decompressed & safe to read (guarded)
	int64_t n_reads;             // reads parsed so far (guarded)
	int done;                    // decompression + parsing finished (guarded)
	int err;
	pthread_mutex_t mtx;
	pthread_cond_t cv;
} decomp_t;

// parse as many complete records as possible from [pos,limit); returns new pos.
static size_t parse_some(decomp_t *d, size_t pos, size_t limit, int final, int64_t *pn)
{
	const char *out = d->out;
	read_t *reads = d->reads;
	int k = d->k;
	int64_t n = *pn;
	while (pos < limit) {
		char c = out[pos];
		if (c == '>' || c == '@') {
			size_t he = find_nl(out, pos, limit);
			if (he == limit) break;                 // header not terminated yet
			size_t ss = he + 1, se, next;
			if (ss > limit) break;
			se = find_nl(out, ss, limit);
			int seq_term = (se < limit);
			if (!seq_term && !final) break;          // wait for the sequence line
			if (c == '@') {                          // FASTQ: also need + and qual lines
				if (!seq_term) break;
				size_t pe = find_nl(out, se + 1, limit);
				if (pe == limit) break;
				size_t qe = find_nl(out, pe + 1, limit);
				if (qe == limit && !final) break;
				next = qe < limit? qe + 1 : qe;
			} else {
				next = seq_term? se + 1 : se;
			}
			int l = (int)(se - ss);
			if (l >= k && n < d->rcap) {
				reads[n].s = out + ss; reads[n].l = l; ++n;
			}
			pos = next;
		} else { // stray/continuation line -> skip
			size_t e = find_nl(out, pos, limit);
			if (e == limit && !final) break;
			pos = e < limit? e + 1 : e;
		}
	}
	*pn = n;
	return pos;
}

// producer: decompress (igzip) and parse concurrently with the consumer's counting
static void *decomp_thread(void *arg)
{
	decomp_t *d = (decomp_t*)arg;
	size_t produced = 0, ppos = 0; // decompressed bytes, parse position
	int64_t n = 0;
	if (!d->is_gzip) { // plain: whole buffer already available
		produced = d->clen;
		ppos = parse_some(d, 0, produced, 1, &n);
		pthread_mutex_lock(&d->mtx);
		d->avail = produced; d->n_reads = n; d->done = 1;
		pthread_cond_broadcast(&d->cv); pthread_mutex_unlock(&d->mtx);
		return 0;
	}
	struct inflate_state st;
	const size_t CHUNK = 4u<<20;
	isal_inflate_init(&st);
	st.crc_flag = ISAL_GZIP;
	st.next_in = d->comp; st.avail_in = d->clen;
	for (;;) {
		size_t room = d->out_cap - produced, want = room < CHUNK? room : CHUNK;
		uint32_t before;
		int r, fin;
		if (want == 0) break;
		st.next_out = (uint8_t*)d->out + produced;
		st.avail_out = (uint32_t)want;
		before = st.avail_out;
		r = isal_inflate(&st);
		produced += before - st.avail_out;
		fin = (r != ISAL_DECOMP_OK) || (st.block_state == ISAL_BLOCK_FINISH) ||
		      (st.avail_in == 0 && before - st.avail_out == 0);
		ppos = parse_some(d, ppos, produced, fin, &n); // parse newly available bytes
		pthread_mutex_lock(&d->mtx);
		d->avail = produced; d->n_reads = n;
		if (r != ISAL_DECOMP_OK) d->err = 1;
		if (fin) d->done = 1;
		pthread_cond_broadcast(&d->cv);
		pthread_mutex_unlock(&d->mtx);
		if (fin) break;
	}
	return 0;
}

typedef struct {
	uint64_t c[256];
} buf_cnt_t;

typedef struct {
	const kc_c4x_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_hist(void *data, long i, int tid)
{
	hist_aux_t *a = (hist_aux_t*)data;
	uint64_t *cnt = a->cnt[tid].c;
	kc_c4_t *g = a->h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k)) {
			int c = kh_key(g, k) & KC_MAX;
			++cnt[c < 255? c : 255];
		}
}

static void print_hist(const kc_c4x_t *h, int n_thread)
{
	hist_aux_t a;
	uint64_t cnt[256];
	int i, j;
	a.h = h;
	CALLOC(a.cnt, n_thread);
	kt_for(n_thread, worker_hist, &a, 1<<h->p);
	for (i = 0; i < 256; ++i) cnt[i] = 0;
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i < 256; ++i)
			cnt[i] += a.cnt[j].c[i];
	free(a.cnt);
	for (i = 1; i < 256; ++i)
		printf("%d\t%ld\n", i, (long)cnt[i]);
}

// dispatch one block of reads [r0,r1) for parallel extraction + insertion
static void run_block(blk_t *bk, int64_t r0, int64_t r1, int n_thread, int p, int P)
{
	int tid, pre;
	bk->r0 = r0; bk->r1 = r1;
	bk->n_chunk = n_thread * 8;
	for (tid = 0; tid < n_thread; ++tid)
		for (pre = 0; pre < P; ++pre)
			bk->buf[((int64_t)tid<<p) + pre].n = 0;
	kt_for(n_thread, extract_worker, bk, bk->n_chunk);
	kt_for(n_thread, insert_worker, bk, P);
}

static kc_c4x_t *count_file(const char *fn, int k, int p, int64_t block_len, int n_thread)
{
	size_t clen, isize;
	uint8_t *comp = slurp(fn, &clen);
	kc_c4x_t *h;
	blk_t bk;
	int P = 1<<p, tid, pre, i;
	decomp_t d;
	pthread_t pth;
	int64_t r_lo = 0, rmin;

	if (comp == 0) return 0;
	memset(&d, 0, sizeof(d));
	d.comp = comp; d.clen = clen; d.k = k;
	d.is_gzip = (clen >= 2 && comp[0] == 0x1f && comp[1] == 0x8b);
	if (d.is_gzip) {
		// ISIZE (uncompressed length mod 2^32) from the gzip footer
		isize = (size_t)comp[clen-4] | (size_t)comp[clen-3]<<8 | (size_t)comp[clen-2]<<16 | (size_t)comp[clen-1]<<24;
		if (isize < clen) isize = clen * 4 + 64;
		d.out = (char*)malloc(isize);
		d.out_cap = isize;
	} else {
		d.out = (char*)comp; d.out_cap = clen; isize = clen;
	}
	// every stored read consumes at least (k+2) bytes, so this bounds the count
	rmin = k + 2 < 4? 4 : k + 2;
	d.rcap = (int64_t)(isize / rmin) + 16;
	MALLOC(d.reads, d.rcap);
	pthread_mutex_init(&d.mtx, 0);
	pthread_cond_init(&d.cv, 0);

	// pre-size tables from the (known) uncompressed size estimate
	h = c4x_init(p);
	{
		int64_t per = (int64_t)isize >> (p + 2);
		khint_t cap = 256;
		while (cap < per) cap <<= 1;
		for (i = 0; i < P; ++i) kc_c4_resize(h->h[i], cap);
	}

	bk.k = k; bk.p = p; bk.n_thread = n_thread; bk.h = h; bk.reads = d.reads;
	CALLOC(bk.buf, (int64_t)n_thread * P);

	pthread_create(&pth, 0, decomp_thread, &d);

	for (;;) {
		int64_t nready; int dn;
		pthread_mutex_lock(&d.mtx);
		// wake when a full block of reads is ready, or when parsing is done
		while (!d.done && d.n_reads - r_lo < block_len / (k + 8) + 1)
			pthread_cond_wait(&d.cv, &d.mtx);
		nready = d.n_reads; dn = d.done;
		pthread_mutex_unlock(&d.mtx);

		// dispatch as many full blocks of reads as we can
		for (;;) {
			int64_t base = 0, r_hi = r_lo;
			while (r_hi < nready && base < block_len) base += d.reads[r_hi++].l;
			if (r_hi == r_lo) break;
			if (base < block_len && !dn) break; // wait for more reads
			run_block(&bk, r_lo, r_hi, n_thread, p, P);
			r_lo = r_hi;
		}
		if (dn && r_lo >= nready) break;
	}

	pthread_join(pth, 0);
	if (d.err) fprintf(stderr, "WARNING: decompression reported an error\n");

	for (tid = 0; tid < n_thread; ++tid)
		for (pre = 0; pre < P; ++pre)
			free(bk.buf[((int64_t)tid<<p) + pre].a);
	free(bk.buf);
	free(d.reads);
	if (d.is_gzip) free(d.out);
	free(comp);
	pthread_mutex_destroy(&d.mtx);
	pthread_cond_destroy(&d.cv);
	return h;
}

int main(int argc, char *argv[])
{
	kc_c4x_t *h;
	int i, c, k = 31, p = KC_BITS, n_thread = 0;
	int64_t block_size = 32000000;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:p:b:t:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'p') p = atoi(o.arg);
		else if (c == 'b') block_size = atoll(o.arg);
		else if (c == 't') n_thread = atoi(o.arg);
	}
	if (n_thread <= 0) { // default: one worker per core, plus one to overlap the decompressor
		long nc = sysconf(_SC_NPROCESSORS_ONLN);
		n_thread = nc > 0? (int)nc + 1 : 4;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c7 [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", p);
		fprintf(stderr, "  -b INT     block size [%ld]\n", (long)block_size);
		fprintf(stderr, "  -t INT     number of worker threads [cores+1]\n");
		return 1;
	}
	if (p < KC_BITS) {
		fprintf(stderr, "ERROR: -p should be at least %d\n", KC_BITS);
		return 1;
	}
	h = count_file(argv[o.ind], k, p, block_size, n_thread);
	if (h == 0) { fprintf(stderr, "ERROR: failed to open/parse '%s'\n", argv[o.ind]); return 1; }
	print_hist(h, n_thread);
	for (i = 0; i < 1<<p; ++i) kc_c4_destroy(h->h[i]);
	free(h->h); free(h);
	return 0;
}
