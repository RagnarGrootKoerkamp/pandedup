# Pandedup

A tool for building a k-mer spectrum of pangenomes.

Takes as input an AGC file, a parameter `k`, and the minimizer-window size `w`,
and outputs a `deduped.fa.zst` containing an SPSS (spectrum preserving string
set, or k-mer spectrum) of the input.
Thus, each k-mer in the input occurs in some contig of the output, and the
output does not contain new k-mers. However, the output is _not_ minimal: some
k-mers will occur twice, and k-mers (nor contigs) are _not_ deduplicated against
their reverse-complement.

**Warning:** This tool is still under development and has not yet been tested
thoroughly. In particular, although I do not expect it, it might be possible
that due to 128-bit hash collisions, some k-mers in the input are not present in
the output.

### HPRCv2

A `k=64` k-mer spectrum of HPRCv2 ([this](https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/B4174A5F-F20E-4DCF-8470-F8A907B640BC--HPRCv2_0.6.1_pr_agc_submission/HPRC_r2_assemblies_0.6.1.agc) AGC file) can be found here: https://zenodo.org/records/21724558

### Usage

Example:

``` sh
pandedup hprcv2.agc -k 64 -w 100 -o hprcv2.spss.k64.fa.zst --threads 64
```

Decreasing `w` gives a smaller output, but requires inversely more memory.
Reducing the number of threads can also help to reduce overall memory usage.


