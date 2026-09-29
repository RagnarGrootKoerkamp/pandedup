# Pandedup

A tool for quick-and-dirty k-mer spectrum preserving deduplication of (human) pangenomes.

Takes as input an `.agc` or `.tar.gz` file, a parameter `k`, and the minimizer-window size `w`,
and outputs a `deduped.fa.zst` containing an SPSS (spectrum preserving string
set, or k-mer spectrum) of the input.
Thus, each k-mer in the input occurs in some contig of the output, and the
output does not contain new k-mers. However, the output is _not_ minimal: some
k-mers will occur twice, and k-mers (nor contigs) are _not_ deduplicated against
their reverse-complement.

In particular, **pandedup is designed as a quick-and-dirty first pass**, to be
used before running a De Bruijn graph (unitig) construction tool like ggcat or
Cuttlefish and then a compression using greedy matchtigs.

### HPRCv2 spectrum

A `k=64` k-mer spectrum of HPRCv2
([this](https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/B4174A5F-F20E-4DCF-8470-F8A907B640BC--HPRCv2_0.6.1_pr_agc_submission/HPRC_r2_assemblies_0.6.1.agc)
AGC file) can be found here: https://zenodo.org/records/21724558 (3.0GB).
Note that this version has a small fraction of "fake" k-mers that were
introduced by wrongly projecting `N` characters to `ACTG`.

A `k=63` version  without ambiguous bases can be found at
https://ragnargrootkoerkamp.nl/upload/hprcv2-k63-greedytigs.fa.zst (2.4GB).
This one was minimized using greedy matchtigs.

### Usage

Typical usage example on a server machine (64 cores; using 50-100GB of memory),
taking 4 minutes:

``` sh
pandedup hprcv2.agc -k 64 -w 100 -o hprcv2.spss.k64.fa.zst --threads 64
```

Decreasing `w` gives a smaller output, but requires inversely more memory.
Reducing the number of threads helps to reduce overall memory usage, as
each thread has a 3GB human genome in memory.
The command also works on my 64GB-memory laptop when using 6 threads.

If you run the output through ggcat anyway, using `w` much smaller than `100`
probably won't help the overall time.

### Ambiguous bases
Input sequences are _split_ on ambiguous IUPAC bases such as `N`, `Y`, and `R`.
Thus, the output contains only `ACGT` bases.

### Canonical
Pandedup can be run both in forward and canonical mode.
In practice, HPRCv2 does not have much duplication between the forward and
reverse-complement strand, and running in forward-only mode is sufficient to
obtain a sufficiently small output that can be fed to `ggcat`.

### Limitations
Although very unlikely, in theory it is possible that the output is missing
some k-mers when two phrases have the same 128-bit hash. In practice, this
is exceedingly unlikely.
