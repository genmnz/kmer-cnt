# Research questions for pushing the counter toward "90 % faster"

These are the bottlenecks I hit after taking `kc-c7` to ~35 % faster than the
previous best (`kc-c4`). Everything below is stated as a pure
computer-science / systems problem — no domain knowledge required to answer.

## Measured context (single machine, 4 physical cores, no SMT)

Cache: L1d 32 KiB/core, L2 1 MiB/core, L3 33 MiB shared. RAM 16 GiB.

Workload per run:
- Input: ~192 MB compressed (single gzip stream) → ~554 MB decompressed text.
- From the text we derive **N ≈ 350 M 64-bit keys** (a stream, produced on the fly).
- Of those, **D ≈ 53 M are distinct**; we need the exact multiplicity of each
  (i.e. a histogram: for each count c, how many distinct keys occur c times).
- Multiplicity is heavily skewed: ~46 M keys occur exactly once (noise), the
  rest cluster around a mode of ~58.

Current wall-time breakdown (4 cores, our best build):
| stage | wall (overlapped) | CPU-seconds |
|---|---|---|
| decompress single gzip stream (igzip, 1 thread) | ~1.5 s, hidden behind counting | 1.5 |
| key extraction (parallel) | ~0.9 s | ~3.9 |
| **aggregation / counting (parallel)** | **~1.7 s** | **~6.6** |
| total | **~3.5 s** | ~12 |

Perfect-scaling floor = CPU/4 ≈ 3.0 s; we are at ~3.5 s. To get dramatically
lower we must **reduce total CPU work**, not just overlap better. The dominant
cost is the counting step, which is memory-latency bound.

---

## Q1 (highest value): cache-efficient parallel count-aggregation of 64-bit keys

**Abstract problem.** A single pass produces a stream of `N ≈ 3.5×10^8` 64-bit
integers (`D ≈ 5.3×10^7` distinct, each 62 significant bits, well-mixed by a
bijective hash so they are effectively uniform random). Output: the exact
occurrence count of every distinct value (we only ever need the *histogram of
counts*, so the distinct values themselves need not be reported). Hardware as
above; 4 threads.

**What we do now.** Partition the stream by the low `p=10` bits into `2^10`
open-addressing tables, each ≈ 1 MB, and probe with a look-ahead software
prefetch:

```c
// keys pre-partitioned into per-bucket buffers a[0..n); table `t` per bucket
for (j = 0; j < n; ++j) {
    if (j + 12 < n) __builtin_prefetch(&t->slot[ hash(a[j+12]) & mask ]);
    insert_or_increment(t, a[j]);        // random access into a ~1 MB table
}
```

At `p=10` each table (~1 MB) does not fit L2, so probes hit L3/DRAM;
aggregation ≈ 2.1 s of CPU-bound-but-latency-limited work. Measured trade-off
as we increase the number of partitions `2^p` (so each table shrinks):

| p (partitions) | scatter/extract | aggregation |
|---|---|---|
| 10 (1 K) | 0.71 s | 2.14 s |
| 12 (4 K) | 0.94 s | 1.69 s |
| 14 (16 K) | 1.25 s | 1.26 s |
| 16 (64 K) | 1.54 s | 1.26 s |

So **small tables (high p) make aggregation cache-resident and ~1.7× cheaper,
but the single-level scatter into many partitions gets proportionally more
expensive** (random writes into thousands of buffer tails → cache/TLB
thrashing). The two effects cancel and total is flat.

**Question.** What is the optimal way to keep *both* the partition-scatter and
the final aggregation cache-resident?

Specifically:
1. Is a **two-level (recursive) radix partition** the right answer here — first
   scatter into a small number `R1` of partitions that fit cache, then within
   each, scatter again into `R2` sub-partitions whose aggregation tables fit
   L1/L2? What `R1`, `R2`, and pass count minimize total traffic for
   `N=3.5×10^8`, `D=5.3×10^7`, the cache sizes above? A worked cost model would
   be ideal.
2. For the scatter itself, do **software write-combining buffers** (accumulate
   a cache-line of keys per partition in a small staging area, flush with
   streaming/non-temporal stores `movntdq`) meaningfully beat plain scattered
   stores at fan-out of hundreds–thousands? Any known crossover fan-out?
3. Is there a fundamentally better structure than open addressing for
   *count-only* aggregation of uniform-random keys at these sizes (e.g. sorting
   networks / radix-sort-then-run-length, given sorted counting is sequential
   and bandwidth-bound rather than latency-bound)? For `N=3.5×10^8` 8-byte keys,
   does a 2–3 pass LSD radix sort + run-length beat prefetched hashing on the
   hardware above?

## Q2: shrinking the per-key slot to double cache residency

Each hash slot is 8 bytes (62-bit key with the count packed into spare low
bits). D ≈ 5.3×10^7 slots at load 0.5 ⇒ ~1 GB of tables — far larger than the
33 MB L3, which is *the* reason aggregation is latency-bound.

**Question.** Can we safely use **4-byte slots** (halving the footprint, so more
fits in cache) without ever merging two distinct keys (which would corrupt the
counts)? The keys are 62-bit and effectively random; within a table of `2^b`
buckets, quotienting removes `b` bits, leaving `62-b` bits to store. For
`b≈17`, that is still 45 bits > 32.

Concretely: is there a compact/quotient-hashing or Cuckoo-filter-style scheme
that guarantees **zero false merges** for `5.3×10^7` uniform 62-bit keys in
32-bit slots — or a hybrid where the common case (small count, short key
remainder) lives in a 4-byte inline slot and only rare collisions spill to an
8-byte overflow area? What is the expected spill rate and net cache win?

## Q3: parallel decompression of a *single* DEFLATE/gzip stream

Right now decompression (SIMD igzip, single thread, ~1.5 s) is *hidden* behind
counting, so it is not our current wall. But it is a hard serial floor: on a
machine with more cores (where counting gets cheaper) or a larger input, this
1.5 s single-threaded stage becomes the limit, and it is the one piece that
cannot currently use more than one core.

**Question.** What is the state of the art for **decompressing one contiguous
DEFLATE stream with multiple threads**? We are aware of speculative approaches
(scan for candidate block boundaries, decompress segments in parallel, then
fix up the back-reference window by a second pass / re-sync). For a ~192 MB
stream on 4–16 cores:
1. What real speedup is achievable, and how much is lost to the resync pass?
2. Is there a way to make the boundary search cheap/robust for arbitrary
   DEFLATE, or does it only work well for certain encoders?
3. If we controlled the *compressor* (we sometimes do), is emitting a
   block-independent format (independent DEFLATE blocks / concatenated members,
   BGZF-style) with an index the pragmatic win, and what block size trades
   compression ratio vs. parallel-decompress granularity best?

## Q4: eliminating the extract↔aggregate barrier (scheduling)

We process the stream in blocks; within a block we run a parallel *extract*
phase, a barrier, then a parallel *aggregate* phase, then the next block. The
barrier plus block-serial structure leaves cores idle: measured wall 3.42 s
vs. the CPU/4 floor of ~2.2 s on the counting-only (pre-decompressed) path.

**Question.** For a two-phase-per-block pipeline where phase B (aggregate) of
block *i* has a true dependency on all of phase A (extract) of block *i*, but
is independent of block *i+1*, what scheduling minimizes idle time on a fixed
pool of T threads?
1. Is **double-buffering** (extract block *i+1* while aggregating block *i*,
   splitting the T threads between the two phases) better than running each
   phase across all T threads with a barrier — given phase A and phase B have
   different per-thread costs (roughly 0.9 s vs 1.7 s of work) and phase B is
   memory-bound while phase A is compute-bound (so co-scheduling them might use
   the core better)?
2. What split of T threads between a memory-bound and a compute-bound phase
   maximizes throughput when they run concurrently and contend for memory
   bandwidth?

---

### What a great answer looks like

For Q1/Q2 especially: a concrete recipe (partition fan-out, pass count, slot
layout) with an order-of-magnitude cost estimate on the cache hierarchy above,
and ideally a small self-contained C snippet for the hot loop. We can wire any
of these into the counter and measure.
