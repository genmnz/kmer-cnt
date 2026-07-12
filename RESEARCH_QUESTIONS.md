# Open problems in a memory-bound streaming count-aggregation pipeline

A brief for a performance-engineering / systems / algorithms researcher. The
problem below is stated in fully abstract terms on purpose — it is a pure
computer-science kernel. **You are encouraged to sweep the web for any
technique, paper, or library; to propose lossless *or* approximate methods; to
suggest SIMD/GPU/other-hardware approaches; and to propose novel experiments.**
"Install a ridiculously fast library" is a fine answer. So is "your framing is
wrong, here is a better one."

Everything here was measured on one fixed machine; numbers are best-of-N,
runs interleaved against the baseline to cancel shared-host noise.

---

## 1. The problem (abstract)

- **Input:** one contiguous compressed byte stream `C` (a single DEFLATE / gzip
  member), `|C| ≈ 1.9×10^8` bytes, that inflates to a text `T`,
  `|T| ≈ 5.5×10^8` bytes. `T` is line-oriented: it alternates short *label*
  lines with *payload* lines; only payload lines matter. Payloads are strings
  over a **4-symbol alphabet**.
- **Key derivation:** slide a window of fixed width `W = 31` over each payload.
  Each window is mapped by a **fixed, cheap, reversible transform** to a
  **62-bit key** (`W` symbols × 2 bits). "Reversible" matters: we may compose
  the key with *any* bijection `f: u64→u64` we like for free, because a
  bijection never merges two distinct keys, so the final answer is unchanged —
  we have total freedom in how keys are hashed/laid out.
- **Output:** the exact **histogram of multiplicities** — for each count `c`,
  how many distinct keys occur exactly `c` times. Note we do **not** need to
  emit the keys themselves, only the count-of-counts. Counts saturate at 1023
  (a key seen ≥1023 times just reports 1023).
- **Scale:** `N ≈ 3.5×10^8` keys produced; `D ≈ 5.3×10^7` distinct.
- **Distribution (important):** heavily skewed. ~`4.6×10^7` keys occur **exactly
  once** (noise-like), and the remaining mass clusters around a mode of ~58.
  So a few ×10^6 "hot" keys account for most of the `N` occurrences, while most
  *distinct* keys are singletons.

The keys can be treated as effectively uniform-random over 62 bits after the
bijection.

## 2. Hardware (the machine we optimize for)

- 4 physical cores, **no SMT**. ~2.8 GHz.
- L1d 32 KiB/core, L2 1 MiB/core, **L3 33 MiB shared**, 16 GiB RAM.
- This 4-core / no-SMT / 33 MiB-L3 shape drives everything below. Answers that
  assume many cores or huge caches should say so.

## 3. Current pipeline and measured breakdown

Best current build ("kc-c7") is ~35% faster than the prior best. Stages:

| stage | what it does | wall (overlapped) | CPU-s |
|---|---|---|---|
| decompress `C→T` | one thread, SIMD inflate | ~1.5 s, hidden | 1.5 |
| parse | find payload spans | ~0.3 s, hidden | 0.3 |
| **extract** | window→key→scatter to buckets | ~0.85 s | ~3.9 |
| **aggregate** | count each distinct key | **~1.6 s** | **~6.6** |
| histogram + I/O | scan tables, print | ~0.2 s | 0.2 |
| **total** | | **~3.5 s** | ~12 |

Structure that already works and should be assumed as the baseline design:
- A **producer thread** inflates + parses concurrently with counting, so
  decompression is almost fully hidden (plain vs. compressed input differ by
  only ~0.34 s of wall).
- **Extraction** partitions keys by their low `p=10` bits into `2^10`
  thread-local buffers (so insertion can be lock-free, one partition per
  worker).
- **Aggregation** inserts each partition's keys into that partition's
  open-addressing table (count packed into the key's spare low bits). Tables are
  **pre-sized** from the known `|T|` (no rehash-resizes) and insertion
  **software-prefetches** the bucket of a look-ahead key.

The hot loops (abstract):

```c
// extract: rolling window -> reversible 62-bit key -> bijection -> bucket
for (each symbol s in payload) {
    x = ((x << 2) | code[s]) & WMASK;          // forward rolling value; a second
                                               // rolling value `twin` tracks the mirror window
    key = min(x, twin);                        // canonicalize: take the smaller of the pair
    h   = bijective_mix(key);                  // free to choose
    buf[h & (P-1)].push(h);                    // scatter into 2^p partitions
}

// aggregate: drain one partition's buffer into its open-addressing table
for (j = 0; j < n; ++j) {
    if (j+PF < n) prefetch(&slot[ index(a[j+PF]) ]);   // hide DRAM latency
    insert_or_increment(table, a[j]);                  // random probe into a ~1 MB table
}
```

## 4. The theoretical floor (roofline) — what "as fast as possible" means

We believe the kernel is **memory bound**, so the ceiling is a roofline:
`t_floor ≈ essential_bytes_moved / effective_bandwidth`, and the whole game is
(a) minimizing essential bytes and (b) converting **random** access (which runs
at *latency-limited* effective bandwidth, far below peak) into **sequential /
cache-resident** access (which runs near *peak* bandwidth).

Essential DRAM traffic of a *cache-optimal* design (order of magnitude):
- read `T` once: ~0.55 GB;
- stream the `N` keys through buffers once (if materialized): ~2×`N`×8 ≈ 5.6 GB
  — a **fused** extract+aggregate that never materializes keys could remove most
  of this, at the cost of locality;
- write the `D` distinct table entries once: ~1 GB.

So essential traffic is roughly **1.5–7 GB** depending on how much is fused vs.
materialized. The current code, by contrast, does **~350M random 64-byte
probes ≈ 22 GB of latency-bound traffic**, which is why aggregation is the wall.

**Two key questions we want your judgment on:**
- (4a) What is the right *essential-traffic* lower bound for this exact problem,
  and which design (fused vs. partition-then-aggregate vs. sort-then-run-length)
  gets closest to it on the hardware in §2?
- (4b) How close to *peak sequential bandwidth* can each candidate design run,
  i.e. what effective bandwidth should we expect, given the random component?

If a cache-optimal design moves ~4 GB near peak bandwidth (~25 GB/s here), the
aggregation floor is ~0.15 s; with parallel decompression (§Q3) the *total*
floor is plausibly ~0.4–0.5 s — i.e. **~90% below the original baseline is
physically motivated, not fantasy** — but only if random DRAM access is
eliminated. Please sanity-check or correct this model.

---

## Q1 (highest value): make the aggregation cache-resident without paying for the scatter

Single-level partitioning into `2^p` buckets faces a hard trade-off we measured:

| p (partitions) | scatter/extract | aggregate |
|---|---|---|
| 10 (1 K, ~1 MB tables) | 0.71 s | 2.14 s |
| 12 (4 K) | 0.94 s | 1.69 s |
| 14 (16 K, ~64 KB tables) | 1.25 s | 1.26 s |
| 16 (64 K) | 1.54 s | 1.26 s |

Small tables (high `p`) make aggregation cache-resident and ~1.7× cheaper, but
the single scatter into thousands of buckets gets proportionally more expensive
(random writes into thousands of buffer tails → cache/TLB thrash). The two
cancel; total is flat, then worsens.

**Question:** what keeps *both* the scatter and the aggregation cache-resident?
1. Optimal **multi-pass radix** fan-out `R1, R2, …` and pass count for
   `N=3.5×10^8`, `D=5.3×10^7`, and §2's caches? A worked cost model + a hot-loop
   snippet would be ideal. *(Naive `1024→×16` two-level was tried and lost — the
   extra pass over all keys cost more than the cache residency it bought; the
   open sub-question is whether the second scatter can be **fused into the
   first** so each key is moved once, or hidden with write-combining below.)*
2. Do **software write-combining buffers** (stage a full cache line per
   partition, flush with non-temporal `movnt` stores) beat plain scattered
   stores at fan-out of hundreds–thousands, and at what crossover fan-out?
3. Is **sort-then-run-length** (LSD radix sort of `N` 8-byte keys, then count
   consecutive runs) — which is sequential/bandwidth-bound rather than
   latency-bound — actually faster here than prefetched hashing? For `N=3.5×10^8`
   on §2's hardware, how many passes, and does it beat ~1.6 s?
4. Given the extreme skew (most *distinct* keys are singletons, most
   *occurrences* are a few ×10^6 hot keys), is a **hot/cold split** worthwhile —
   a tiny always-cache-resident structure absorbing the hot keys' many
   increments, and a separate cheap path for the singleton-heavy tail?

## Q2: shrink the per-slot footprint to double cache residency

Slots are 8 bytes (62-bit key with the count in spare low bits). `D≈5.3×10^7`
slots at load 0.5 ⇒ ~1 GB of tables ≫ 33 MiB L3 — *the* reason aggregation is
latency-bound.

**Question:** can we safely use **4-byte slots** (halving footprint, so far more
fits in cache) with **zero false merges** of distinct keys (a merge corrupts the
counts)? Within a table of `2^b` buckets, quotienting removes `b` bits, leaving
`62−b`; for `b≈17` that is still 45 bits > 32. Is there a compact/quotient- or
cuckoo-filter-style layout that guarantees no merges for `5.3×10^7` uniform
62-bit keys in 32-bit slots — or a **hybrid** (common case: short remainder +
small count inline in 4 bytes; rare collisions spill to an 8-byte side table)?
Expected spill rate and net cache win?

## Q3: parallel decompression of a *single* DEFLATE stream

Decompression (SIMD inflate, one thread, ~1.5 s) is currently *hidden* behind
counting, but it is a hard serial floor: with more cores, or once aggregation is
cheaper, it becomes the wall, and it is the one stage stuck on a single core.

**Question:** state of the art for decompressing **one contiguous DEFLATE stream
with N threads**? We know of speculative approaches (scan for candidate block
boundaries, inflate segments in parallel, fix up the back-reference window in a
second pass). For a ~190 MB stream on 4–16 cores: realistic speedup, cost of the
resync pass, robustness for arbitrary encoders? And if we *control the
compressor* (sometimes we do): is emitting a block-independent / indexed format
(independent blocks or concatenated members) the pragmatic win, and what block
size best trades ratio vs. parallel granularity? Any drop-in library?

## Q4: eliminate idle time between the two phases (scheduling)

Per block we run a parallel **extract** (compute-bound), a barrier, then a
parallel **aggregate** (memory-latency-bound). Phase timings on 5 workers:
`wait 0.24 + extract 0.82 + aggregate 1.59 = 2.66 s`, i.e. the phases are
*sequential*. Aggregation is latency-bound (it sped up markedly going 4→5
oversubscribed threads — idle cores waiting on DRAM), while extraction is
compute-bound. In principle, running extract(i+1) *during* aggregate(i) should
let compute fill the memory-stall cycles.

**Question:** what schedule actually achieves that on a **fixed pool of 4 cores,
no SMT**? We tried one work-stealing pool over combined
(aggregate-i ∪ extract-i+1) tasks → **no gain** (the pool runs all threads on
one phase at a time; aggregation seems to saturate the memory system so
concurrent extraction doesn't proceed "for free"). Untried: a **static split**
of cores between a memory-bound and a compute-bound phase — but that halves each
phase's parallelism and both still share one memory controller. Is there a split
(and prefetch-distance / MLP tuning) that wins, or is this memory system simply
saturated during aggregation (in which case the answer is "don't bother — attack
Q1/Q2 instead")?

---

## 5. Where we think we're near-maxed (please confirm or break)

These are the parts we believe are close to their floor on §2's hardware; we
want either confirmation or a counter-idea:
- **Decompression throughput per thread** — already SIMD inflate (~1 GB/s). Only
  parallelism (Q3) can help, not a faster serial codec (we tried two).
- **Extraction compute** — ~0.85 s wall / ~3.9 CPU-s for 3.5×10^8 windows; the
  rolling-window + reversible-key + scatter is cheap; a cheaper bijection made
  no measurable difference. SIMD encoding of the 2-bit symbols is possible but
  the rolling dependency + per-key scatter looks inherently serial-ish. Is there
  a vectorized/branch-free formulation of "rolling window → min-of-pair 62-bit key
  → partition" that is materially faster?
- **Overlap of decompression with counting** — already ~fully hidden.

## 6. Approaches already tried and rejected (with outcomes)

So the search is not repeated:
- **Custom single-array open-addressing table** (count in key, no used-bitmap):
  *slower* than the library table, because the library's compact ~16 MB
  used-bitmap fits L3 and answers most probes without touching the 1 GB key
  array; a single array touches DRAM on every probe.
- **More partitions (larger `p`)**: net flat then worse (Q1 table).
- **Naive two-level radix** (`1024→×16`): extra pass costs more than it saves.
- **2 MB huge pages (`MADV_HUGEPAGE`) for the tables**: ~4% on a 1 GB random
  RMW microbenchmark — software prefetch already hides the page-walk, so TLB was
  not the bottleneck.
- **Software-pipelining extract(i+1) with aggregate(i)** on one work-stealing
  pool: no gain (Q4).
- **Skipping teardown `free()`s before exit**: no measurable gain here.

## 7. What we're looking for

Concrete, testable direction, ranked by expected payoff, ideally with a small
self-contained C snippet for any proposed hot loop and an order-of-magnitude
estimate against §2's cache hierarchy and §4's roofline. Naming a specific
library/algorithm to drop in (radix sort, compact-hash, parallel-inflate, SIMD
codec, …) is welcome — we will wire it in and benchmark. Novel or unconventional
ideas welcome; so is "the roofline model in §4 is wrong, here's why." The bar to
beat is **~1.6 s of aggregation** (the dominant term) and the **~3.5 s total**;
the stretch target is the roofline floor of ~0.4–0.5 s.
