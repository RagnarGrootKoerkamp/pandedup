use std::ops::Range;

use gxhash::gxhash128;
use simd_minimizers::packed_seq::AsciiSeq;

fn main() {
    let path = "/home/philae/git/eth/data/hprcv2.agc";
    use ragc_core::{Decompressor, DecompressorConfig};

    // Open an archive
    let config = DecompressorConfig::default();
    let mut decompressor = Decompressor::open(path, config).unwrap();

    // List available samples
    let samples = decompressor.list_samples();
    // println!("Found {} samples", samples.len());
    // eprintln!("samples: {samples:?}");

    let mut seen = std::collections::HashSet::new();
    let k = 40; // overlap
    let w = 500; // at most w apart
    let l = k + w - 1; // syncmer length
    let syncmers = simd_minimizers::closed_syncmers(k, w);

    let mut num_syncmers = 0;
    let mut taken = 0;
    let mut skipped = 0;
    let mut short = 0;
    let mut short_bp = 0;
    let mut prefix_bp = 0;
    let mut suffix_bp = 0;
    let mut output_bp = 0;
    let mut skipped_bp = 0;
    let mut input_bp = 0;

    let mut cnt = 0;
    let mut syncmer_poss = vec![];
    for sample in &samples {
        eprintln!("sample {}", sample);
        let contigs = decompressor.list_contigs(sample).unwrap();
        for contig in &contigs {
            cnt += 1;
            let seq = decompressor.get_contig(&sample, &contig).unwrap();
            input_bp += seq.len();

            eprintln!("{cnt:>3} ({:>10}): contig {contig}", seq.len());

            // Minimizer-based PFP instead?
            syncmer_poss.clear();
            syncmers.run(AsciiSeq(&seq), &mut syncmer_poss);

            if syncmer_poss.is_empty() {
                short += 1;
                short_bp += seq.len();
                continue;
            }

            // sanity check
            for &[p, q] in syncmer_poss.array_windows() {
                assert!(p < q, "{p} < {q}");
                assert!(q <= p + w as u32, "{q} <= {p} + {w}");
                // assert!(p as usize + w <= seq.len(), "{p} + {w} <= {}", seq.len());
                assert!(p as usize + l <= seq.len(), "{p} + {l} <= {}", seq.len());
            }
            // assert!(*syncmer_poss.last().unwrap() as usize + l >= seq.len());

            // Emit the prefix.
            let mut active = 0..syncmer_poss[0] as usize + k - 1;
            let mut push = |range: Range<usize>| {
                assert!(range.end >= active.end);
                if range.start <= active.end {
                    active.end = range.end;
                } else {
                    output_bp += active.len();
                    skipped_bp += range.start - active.end;
                    active = range;
                }
            };
            let prefix = &seq[..syncmer_poss[0] as usize + k - 1];
            prefix_bp += prefix.len();
            // Emit the syncmers.
            for &p in syncmer_poss.iter() {
                let p = p as usize;
                num_syncmers += 1;
                let syncmer = &seq[p..p + l];
                // output_bp += syncmer.len();
                let hash = gxhash128(syncmer, 0);
                if seen.insert(hash) {
                    // eprintln!("take {p} hash={hash}");
                    taken += 1;
                    push(p..p + l);
                } else {
                    // eprintln!("skip {p} hash={hash}");
                    skipped += 1;
                }
            }
            // Emit the suffix.
            let suffix_range = syncmer_poss[syncmer_poss.len() - 1] as usize + w..seq.len();
            let suffix = &seq[suffix_range.clone()];
            suffix_bp += suffix.len();
            push(suffix_range);
            output_bp += active.len();
        }
        eprintln!("taken:   {taken}");
        eprintln!("skipped: {skipped}");
        eprintln!("short:   {short}");
        eprintln!("num_syncmers: {num_syncmers}");
        eprintln!("prefix_bp:  {prefix_bp}");
        eprintln!("suffix_bp:  {suffix_bp}");
        eprintln!("input_bp:   {input_bp}");
        eprintln!("output_bp:  {output_bp}");
        eprintln!("skipped_bp: {skipped_bp}");
        eprintln!("short_bp:   {short_bp}");
    }
}
