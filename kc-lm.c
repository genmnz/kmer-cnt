// kc-lm: a low-memory k-mer counter for constrained / low-end devices.
//
// Same result as kc-c4 / kc-c7 (identical count histogram of canonical <=32-mers)
// but designed to run in a small, *bounded* RAM footprint no matter how large the
// input is, by offloading intermediate data to storage. It trades speed for a
// tiny, predictable memory profile — it will not swap-thrash or OOM on a machine
// with little RAM.
//
// Strategy: external (disk-partitioned) counting.
//   Pass 1  — stream the input (zlib, no whole-file buffering), extract canonical
//             k-mers, and append each (as a 64-bit hashed value) to one of P
//             temporary partition files chosen by the low bits of the hash. Only
//             small per-partition write buffers live in RAM.
//   Pass 2  — process one partition file at a time: read it back, count distinct
//             values in a hash table that holds just that partition (≈ D/P
//             entries), fold into the histogram, free it, delete the file.
//
// Peak RAM ≈ max(P * write_buffer, one partition's hash table) — a few tens of MB
// by default, independent of input size. Temp disk ≈ 8 bytes * (number of k-mers).
// Single-threaded on purpose: predictable, no per-thread memory multiplier.
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "ketopt.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khashl.h"
// map: hashed-kmer -> occurrence count
KHASHL_MAP_INIT(, kclm_t, kclm, uint64_t, uint32_t, kh_hash_dummy, kh_eq_generic)

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
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

static inline uint64_t mix64(uint64_t x) // bijective; keeps distinct keys distinct
{
	x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
	x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
	x ^= x >> 33;
	return x;
}

typedef struct { // a buffered writer for one partition
	FILE *fp;
	uint64_t *buf;
	int n, m;
	char *path;
} part_t;

static void part_flush(part_t *w)
{
	if (w->n) { fwrite(w->buf, 8, w->n, w->fp); w->n = 0; }
}

static inline void part_push(part_t *w, uint64_t y)
{
	w->buf[w->n++] = y;
	if (w->n == w->m) part_flush(w);
}

// extract canonical k-mers from one sequence and route them to partitions
static void emit_seq(part_t *parts, int p_bits, int k, int len, const char *seq)
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2, pmask = (1ULL<<p_bits) - 1;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) {
			x[0] = (x[0] << 2 | c) & mask;
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			if (++l >= k) {
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				uint64_t h = mix64(y);
				part_push(&parts[h & pmask], h);
			}
		} else l = 0, x[0] = x[1] = 0;
	}
}

// count k-mers of one sequence that fall in residue class $r modulo $R (multi-pass mode)
static void count_seq_slice(kclm_t *h, uint64_t R1, uint64_t r, int k, int len, const char *seq)
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) {
			x[0] = (x[0] << 2 | c) & mask;
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			if (++l >= k) {
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				uint64_t hh = mix64(y);
				if ((hh & R1) == r) { // only this pass's slice
					int absent;
					khint_t it = kclm_put(h, hh, &absent);
					if (absent) kh_val(h, it) = 1; else ++kh_val(h, it);
				}
			}
		} else l = 0, x[0] = x[1] = 0;
	}
}

// Multi-pass strategy: no temp files at all; re-read the input 2^m_bits times,
// each pass counting only the k-mers whose hash lands in one residue class, so
// the in-RAM table only ever holds ~D/2^m_bits entries.
static int run_multipass(const char *fn, int k, int m_bits, uint64_t *hist)
{
	uint64_t R = 1ULL << m_bits, R1 = R - 1, r;
	int i;
	for (i = 0; i < 256; ++i) hist[i] = 0;
	for (r = 0; r < R; ++r) {
		gzFile fp = gzopen(fn, "r");
		kseq_t *ks;
		kclm_t *h;
		khint_t it;
		if (fp == 0) { fprintf(stderr, "ERROR: cannot open '%s'\n", fn); return 1; }
		ks = kseq_init(fp);
		h = kclm_init();
		while (kseq_read(ks) >= 0)
			if ((int)ks->seq.l >= k)
				count_seq_slice(h, R1, r, k, ks->seq.l, ks->seq.s);
		for (it = 0; it < kh_end(h); ++it)
			if (kh_exist(h, it)) { uint32_t v = kh_val(h, it); ++hist[v < 255? v : 255]; }
		kclm_destroy(h);
		kseq_destroy(ks);
		gzclose(fp);
	}
	return 0;
}

int main(int argc, char *argv[])
{
	int c, k = 31, p_bits = 8, buf_kb = 64, i, P, m_bits = -1;
	const char *tmpdir = 0;
	ketopt_t o = KETOPT_INIT;
	gzFile fp;
	kseq_t *ks;
	part_t *parts;
	uint64_t hist[256];
	long pid;

	while ((c = ketopt(&o, argc, argv, 1, "k:p:b:d:M:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'p') p_bits = atoi(o.arg);
		else if (c == 'b') buf_kb = atoi(o.arg);
		else if (c == 'd') tmpdir = o.arg;
		else if (c == 'M') m_bits = atoi(o.arg);
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-lm [options] <in.fa[.gz]>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT   k-mer size (<=32) [%d]\n", k);
		fprintf(stderr, "  -M INT   no-temp-disk mode: re-read input 2^INT times, RAM ~ D/2^INT\n");
		fprintf(stderr, "           (default is disk-partition mode, which reads input once)\n");
		fprintf(stderr, "  -p INT   disk mode: log2 of the number of on-disk partitions [%d]\n", p_bits);
		fprintf(stderr, "  -b INT   disk mode: per-partition write buffer in KB [%d]\n", buf_kb);
		fprintf(stderr, "  -d STR   disk mode: directory for temporary files [$TMPDIR or .]\n");
		fprintf(stderr, "\nBoth modes give identical output in bounded RAM. Disk mode is faster but\n");
		fprintf(stderr, "writes ~8*(#kmers) bytes of temp files; -M uses no temp disk but re-reads\n");
		fprintf(stderr, "the input 2^INT times. Lower -p/-b (disk) or raise -M (multi-pass) for less RAM.\n");
		return 1;
	}
	if (k < 1 || k > 32) { fprintf(stderr, "ERROR: -k must be in [1,32]\n"); return 1; }
	if (m_bits >= 0) { // no-temp-disk multi-pass mode
		if (m_bits > 30) { fprintf(stderr, "ERROR: -M too large\n"); return 1; }
		if (run_multipass(argv[o.ind], k, m_bits, hist)) return 1;
		for (i = 1; i < 256; ++i) printf("%d\t%ld\n", i, (long)hist[i]);
		return 0;
	}
	if (p_bits < 0 || p_bits > 20) { fprintf(stderr, "ERROR: -p must be in [0,20]\n"); return 1; }
	if (tmpdir == 0) { tmpdir = getenv("TMPDIR"); if (tmpdir == 0) tmpdir = "."; }
	P = 1 << p_bits;
	pid = (long)getpid();

	// ---- Pass 1: stream input, partition k-mers to temp files ----
	parts = (part_t*)calloc(P, sizeof(part_t));
	for (i = 0; i < P; ++i) {
		char nm[4096];
		snprintf(nm, sizeof(nm), "%s/kclm.%ld.%d.tmp", tmpdir, pid, i);
		parts[i].path = strdup(nm);
		parts[i].fp = fopen(nm, "wb");
		if (parts[i].fp == 0) { fprintf(stderr, "ERROR: cannot create temp file '%s'\n", nm); return 1; }
		parts[i].m = (buf_kb * 1024) / 8; if (parts[i].m < 1) parts[i].m = 1;
		parts[i].buf = (uint64_t*)malloc(parts[i].m * 8);
	}

	if ((fp = gzopen(argv[o.ind], "r")) == 0) { fprintf(stderr, "ERROR: cannot open '%s'\n", argv[o.ind]); return 1; }
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0)
		if ((int)ks->seq.l >= k)
			emit_seq(parts, p_bits, k, ks->seq.l, ks->seq.s);
	kseq_destroy(ks);
	gzclose(fp);
	for (i = 0; i < P; ++i) { part_flush(&parts[i]); fclose(parts[i].fp); free(parts[i].buf); parts[i].buf = 0; }

	// ---- Pass 2: count each partition, one at a time ----
	for (i = 0; i < 256; ++i) hist[i] = 0;
	{
		int rm = (buf_kb * 1024) / 8; if (rm < 4096) rm = 4096;
		uint64_t *rbuf = (uint64_t*)malloc(rm * 8);
		for (i = 0; i < P; ++i) {
			kclm_t *h = kclm_init();
			FILE *rf = fopen(parts[i].path, "rb");
			size_t got;
			khint_t it;
			if (rf) {
				while ((got = fread(rbuf, 8, rm, rf)) > 0) {
					size_t j;
					for (j = 0; j < got; ++j) {
						int absent;
						it = kclm_put(h, rbuf[j], &absent);
						if (absent) kh_val(h, it) = 1;
						else ++kh_val(h, it);
					}
				}
				fclose(rf);
			}
			for (it = 0; it < kh_end(h); ++it)
				if (kh_exist(h, it)) {
					uint32_t v = kh_val(h, it);
					++hist[v < 255? v : 255];
				}
			kclm_destroy(h);
			remove(parts[i].path);
			free(parts[i].path);
		}
		free(rbuf);
	}
	free(parts);

	for (i = 1; i < 256; ++i) printf("%d\t%ld\n", i, (long)hist[i]);
	return 0;
}
