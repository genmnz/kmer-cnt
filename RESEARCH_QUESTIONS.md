# Open problems: a memory-bound high-cardinality frequency-counting kernel

A brief for a performance-engineering / systems / algorithms researcher. You are
encouraged to search the web for any technique, paper, or library; to propose
exact *or* approximate methods; to suggest SIMD/GPU/other-hardware approaches;
and to propose new experiments. "Drop in this fast library" is a fine answer.
So is "your model is wrong, here is a better one."

All numbers were measured on one fixed machine; they are best-of-N with runs
interleaved against the baseline to cancel shared-host noise.

---

## 1. The kernel

- **Input:** one contiguous compressed stream `C` (standard DEFLATE),
  `|C| ≈ 1.9×10^8` bytes, that inflates to an in-memory buffer `T`,
  `|T| ≈ 5.5×10^8` bytes.
- **Producer:** a cheap rolling map scans `T` and emits a stream of
  `N ≈ 3.5×10^8` **64-bit keys** (each has ≤ 62 significant bits). The map is a
  few ALU ops per key and is *not* the bottleneck. We are free to compose the
  keys with **any bijection** `f: u64→u64` at no correctness cost (a bijection
  never merges two distinct keys), so we have total freedom in how keys are
  hashed and laid out. Treat the keys as effectively uniform-random over 62 bits.
- **Output:** the exact **frequency histogram** — for each frequency `c`, how
  many distinct keys occur exactly `c` times. We do **not** need to emit the keys
  themselves, only this count-of-counts. Frequencies saturate at 1023.
- **Scale:** `N ≈ 3.5×10^8` keys; `D ≈ 5.3×10^7` distinct.
- **Distribution (important):** heavily skewed. ~`4.6×10^7` keys occur **exactly
  once**; the remaining mass concentrates near a mode of ~58. So a few ×10^6
  "hot" keys account for most of the `N` occurrences, while most *distinct* keys
  occur once.

## 2. Hardware (the target)

- 4 physical cores, **no SMT**. ~2.8 GHz.
- L1d 32 KiB/core, L2 1 MiB/core, **L3 33 MiB shared**, 16 GiB RAM.
- This 4-core / no-SMT / 33 MiB-L3 shape drives everything below; answers that
  assume many cores or huge caches should say so.

## 3. Current pipeline and measured breakdown

The best current build is ~35% faster than the prior baseline. Stages:

| stage | what it does | wall (overlapped) | CPU-s |
|---|---|---|---|
| inflate `C→T` | one thread, SIMD | ~1.5 s, hidden | 1.5 |
| scan `T` | find the spans to map over | ~0.3 s, hidden | 0.3 |
| **produce** | rolling map → key → scatter to partitions | ~0.85 s | ~3.9 |
| **aggregate** | count each distinct key | **~1.6 s** | **~6.6** |
| histogram + I/O | scan tables, print | ~0.2 s | 0.2 |
| **total** | | **~3.5 s** | ~12 |

Design that already works and is assumed as the baseline:
- A **producer thread** inflates + scans concurrently with counting, so inflate
  is almost fully hidden (uncompressed vs compressed input differ by only
  ~0.34 s of wall).
- Keys are **partitioned** by their low `p=10` bits into `2^10` thread-local
  buffers, so aggregation is lock-free (one partition per worker).
- Aggregation inserts each partition's keys into that partition's
  open-addressing table (frequency packed into the key's spare low bits). Tables
  are **pre-sized** from the known `|T|` (no rehash-resizes) and insertion
  **software-prefetches** the bucket of a look-ahead key.

Hot loops (representative):

```c
// produce: rolling map -> 62-bit key -> optional bijection -> partition
for (;;) {
    uint64_t key = roll(next);      // few ALU ops
    uint64_t h   = mix(key);        // any bijection, free to choose
    buf[h & (P-1)].push(h);         // scatter into 2^p partitions
}

// aggregate: drain one partition's buffer into its open-addressing table
for (j = 0; j < n; ++j) {
    if (j+PF < n) prefetch(&slot[ index(a[j+PF]) ]);  // hide DRAM latency
    insert_or_increment(table, a[j]);                 // random probe, ~1 MB table
}
```

## 4. The roofline — what "as fast as possible" means

We believe the kernel is **memory bound**, so the ceiling is a roofline:
`t_floor ≈ essential_bytes_moved / effective_bandwidth`. The lever is that both
terms depend on the algorithm. In particular a machine has two very different
"memory speeds":

- **sequential bandwidth** `B_seq` — streaming contiguous bytes; high (~25 GB/s
  on this box);
- **random-access throughput** `B_rand` — chasing scattered 64-byte lines;
  latency-limited, *much* lower (we measured ~30M lines/s/core ≈ 1.9 GB/s;
  ~4 cores ≈ 7–8 GB/s).

The current aggregation does ~`3.5×10^8` random 64-byte probes ≈ **22 GB of
random traffic** at ~7 GB/s ⇒ ~2–3 s — this is the wall. A design that keeps the
per-partition tables **cache-resident** turns those probes into cache hits, so
DRAM sees only sequential streams: load `T` once (~0.55 GB) + stream the keys
once (~2.8 GB) + write the `D` table entries once (~1 GB) ≈ 4.4 GB near `B_seq`.

Two questions we want judged:
- **(4a)** What is the right essential-traffic lower bound for this exact
  problem, and which design gets closest on §2's hardware?
- **(4b)** What effective bandwidth should each candidate actually achieve, given
  the residual random component? Even a perfectly cache-resident hash still does
  `N` discrete table ops (~a few ns each): `3.5×10^8` × ~2 ns / 4 cores ≈ 0.18 s.
  Is ~0.2 s the real aggregation floor, or is even that optimistic?

If a cache-resident design streams ~4 GB near `B_seq` (~0.18 s) and inflate is
parallelized (§Q3), the total floor is plausibly **~0.3–0.5 s** — i.e. ~10× the
current baseline is physically motivated, *if* random DRAM access is eliminated.
Please sanity-check or correct this model.

---

## Q1 (highest value): cache-resident aggregation without paying for the scatter

Single-level partitioning into `2^p` buckets has a hard measured trade-off:

| p (partitions) | scatter (produce) | aggregate |
|---|---|---|
| 10 (1 K, ~1 MB tables) | 0.71 s | 2.14 s |
| 12 (4 K) | 0.94 s | 1.69 s |
| 14 (16 K, ~64 KB tables) | 1.25 s | 1.26 s |
| 16 (64 K) | 1.54 s | 1.26 s |

Small tables (high `p`) make aggregation cache-resident and ~1.7× cheaper, but
the single scatter into thousands of buckets gets proportionally more expensive
(random writes into thousands of buffer tails whose working set exceeds cache).
The two cancel; total is flat, then worse. **We have tried the obvious escapes
and each lost — see the measured results below.** So the real question stands:

**Is there any way to make the aggregation cache-resident whose extra cost is
less than the ~0.8 s it saves?** Concretely:
1. An optimal **multi-pass radix** fan-out and pass count that we have not found,
   ideally where the second partition pass is **fused into the first** so each
   key moves once. (Naive two-level, and sorting, both lost — below.)
2. A cheaper **cache-resident hash** than open addressing for uniform keys — or a
   layout where the hot working set provably fits L2 at `p=10`.
3. Using the **skew**: most *distinct* keys occur once, most *occurrences*
   are a few ×10^6 hot keys. Is a hot/cold split worthwhile — a tiny always-
   cache-resident structure absorbing the hot keys' many increments, and a
   separate cheap path for the singleton-heavy tail?

## Q2: shrink the per-slot footprint to double cache residency

Slots are 8 bytes (62-bit key with the frequency in spare low bits).
`D≈5.3×10^7` slots at load 0.5 ⇒ ~1 GB of tables ≫ 33 MiB L3 — *the* reason
aggregation is latency-bound. Can we use **4-byte slots** (halving footprint,
so far more fits in cache) with **zero merges** of distinct keys (a merge
corrupts the histogram)? Within `2^b` buckets, quotienting removes `b` bits,
leaving `62−b`; for `b≈17` that is still 45 bits > 32. Is there a compact- or
cuckoo-style layout guaranteeing no merges for `5.3×10^7` uniform 62-bit keys in
32-bit slots — or a **hybrid** (common case: short remainder + small frequency
inline in 4 bytes; rare collisions spill to an 8-byte side table)? Expected spill
rate and net cache win?

## Q3: parallelize inflating a single DEFLATE stream

Inflate (SIMD, one thread, ~1.5 s) is currently *hidden* behind counting, but it
is a hard serial floor: with more cores, or once aggregation is cheaper, it
becomes the wall, and it is the one stage stuck on a single core. What is the
state of the art for inflating **one contiguous DEFLATE stream with N threads**?
For a ~190 MB stream on 4–16 cores: realistic speedup, robustness, and any
drop-in library? And if we *control the producer of `C`* (sometimes we do): is
emitting an independently-decodable / indexed layout the pragmatic win, and what
block size best trades ratio vs. parallel granularity?

## Q4: remove the idle time between the two phases

Per block we run a parallel **produce** (compute-bound), a barrier, then a
parallel **aggregate** (memory-latency-bound). On 5 workers:
`wait 0.24 + produce 0.82 + aggregate 1.59 = 2.66 s` — sequential. Aggregation
is latency-bound (it sped up markedly 4→5 oversubscribed threads: idle cores
waiting on memory), while produce is compute-bound. Running produce(i+1) *during*
aggregate(i) should let compute fill the memory-stall cycles. What schedule
actually achieves that on a **fixed 4-core, no-SMT pool**? A single work-stealing
pool over combined tasks gave **no gain** (all threads run one phase at a time;
aggregation appears to saturate the memory system so concurrent produce does not
proceed for free). A static core split is untried but halves each phase's
parallelism and both still share one memory controller. Is there a split (and
prefetch-distance / MLP tuning) that wins, or is the memory system simply
saturated during aggregation (answer: leave Q4, spend effort on Q1/Q2)?

---

## 5. Where we appear near-maxed (please confirm or break)

- **Inflate throughput per thread** — already SIMD (~1 GB/s of output). Only
  parallelism (Q3) can help, not a faster serial codec (two were tried).
- **Produce compute** — ~0.85 s wall / ~3.9 CPU-s for `3.5×10^8` keys; the
  rolling map + scatter is cheap and a cheaper bijection made no measurable
  difference. Is there a vectorized/branch-free formulation of "rolling map →
  62-bit key → partition" that is materially faster?
- **Overlap of inflate with counting** — already ~fully hidden.

## 6. Approaches already implemented and measured (so the search isn't repeated)

Worked (cumulative, in the current build, ~35% faster overall):
- SIMD inflate instead of the stock codec: `C→T` 4.6 s → ~1.1–1.5 s.
- Producer thread doing inflate + scan overlapped with counting.
- Parallelize the produce stage (was single-threaded).
- Pre-size tables from the known `|T|` (no resizes): aggregate 3.0 s → ~2.1 s.
- Software-prefetch the look-ahead bucket: aggregate → ~1.6 s.

Also **worked** (from a second review round):
- **Tune the load factor to ~0.4** (was ~0.2). Counterintuitively, in the *real
  parallel* run *smaller* tables (higher load, up to ~0.4–0.5) win, even though a
  single-thread microbench shows the opposite. The parallel cost is dominated by
  table allocation / first-touch, shared-L3 pressure across workers, and the
  final histogram scan — all ∝ total table bytes — not by probe-chain length.
  This alone took the build from ~35% → ~44% faster than baseline; still exact.

Did **not** help (measured):
- **Custom single-array open-addressing table** (frequency in key, no occupancy
  bitmap): *slower* — the library's compact ~16 MB occupancy bitmap fits L3 and
  answers most probes without touching the 1 GB key array; a single array hits
  DRAM every probe.
- **More partitions (larger `p`)**: net flat then worse (Q1 table).
- **Sort-then-run-length** (partition, then LSD-radix-sort each partition, then
  count runs): **~6× slower** aggregation than hashing (8.3 s vs 1.3 s
  single-thread). The multi-pass sort moves far more bytes than the random
  hashing that prefetch already hides. Sorting is *not* the answer here.
- **Write-combining scatter**: at high fan-out (`2^14`) it cuts the scatter ~35%
  but the total still loses to plain `2^10`; at the production `2^10` fan-out it
  is a **regression** (~3.8 s → ~4.6 s end-to-end) because the ≤2^10 buffer tails
  already fit L2, so the staging is pure overhead (and is itself a scattered
  write). WCB only pays when the tail set exceeds cache.
- **Compact 4-byte slots**: (a) provably impossible to keep *exact* for 5.3×10^7
  uniform 62-bit keys (birthday bound needs ≥61 preserved bits); and (b) even the
  *approximate* form does **not** cut the dominant cost — aggregation is bounded
  by ~3.5×10^8 *random 64-byte cache-line* touches, and a narrower slot does not
  reduce cache-lines-touched-per-random-probe. It only densifies the singleton
  first-writes (~46M) from ~368 MB to ~186 MB, negligible against ~22 GB. So Q2's
  premise (smaller slot → less traffic) does not hold for random-access probing.
- **Lower load factor to fit L2** (the Q1×Q2 "residency" thesis): refuted — see
  the load-factor result above; bigger tables are *slower* in parallel.
- **2 MB pages (`MADV_HUGEPAGE`) for the tables**: ~4% on a 1 GB random-update
  microbenchmark — prefetch already hides the page-walk.
- **Pipelining / overlapping produce with aggregate (Q4)**: no gain. Aggregation
  is *latency*-bound and scales with thread count (2→5 threads: 2.74→1.63 s), so
  it needs *all* cores; a static core split starves it. Confirmed dead on 4 cores.
- **Hot/cold heavy-hitter split (Q1.3)**: not implemented — the reuse-distance
  argument (every "hot" key recurs ~N/58 ≈ 6×10^6 keys apart, ≫ any cache) says a
  front cache turns over completely between visits. Accepted as a dead end.
- **Block size**: 32 MB is the measured optimum (smaller adds overhead, larger
  loses producer/consumer overlap).
- **Skipping teardown frees before exit**: no measurable gain.

**Net:** aggregation is bounded by ~3.5×10^8 *random cache-line* probes; on
§2's hardware every attempt to reach a genuinely cache-resident regime either
breaks exactness or costs more in scatter/allocation than it saves. The load
factor and prefetch depth are the knobs that actually move it. The open crux is
whether the *count* of random probes can be reduced for exact counting, or the
random component converted to sequential without a net-losing extra pass.

## 7. What we're looking for

Concrete, testable direction ranked by expected payoff, ideally with a small
self-contained C snippet for any proposed hot loop and an order-of-magnitude
estimate against §2's caches and §4's roofline. Naming a specific
library/algorithm to drop in is welcome — we will wire it in and benchmark.
Novel or unconventional ideas welcome; so is "the roofline in §4 is wrong, here
is why." The bar to beat is **~1.6 s of aggregation** (the dominant term) and
the **~3.5 s total**; the stretch target is the roofline floor of ~0.3–0.5 s.
