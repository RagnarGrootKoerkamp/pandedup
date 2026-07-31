use clap::Parser;
use fxhash::{FxHashMap, FxHashSet};
use ragc_core::{Decompressor, DecompressorConfig};
use std::{
    io::{BufWriter, Write},
    ops::Range,
    path::PathBuf,
    sync::{atomic::AtomicUsize, Mutex, RwLock},
    time::Duration,
};

fn hasher(seq: &[u8]) -> u128 {
    // const SEED: i64 = 1983274983247984327;
    // gxhash::gxhash128(seq, SEED)
    xxhash_rust::xxh3::xxh3_128(seq)
}

use seq_hash::AntiLexHasher;
use simd_minimizers::packed_seq::AsciiSeq;

/// Build a non-minimal SPSS (spectrum-preserving string set, or k-mer spectrum) from an .agc file.
#[derive(clap::Parser)]
struct Args {
    /// Input .agc file.
    input: PathBuf,
    /// Output path. Defaults to `input.dedup.fa.zst`.
    #[clap(short)]
    output: Option<PathBuf>,
    /// Build a k-mer spectrum.
    #[clap(short, default_value = "64")]
    k: usize,
    /// Window-size for minimizer-phrases. Small w shrink the output but need more memory.
    #[clap(short, default_value = "100")]
    w: usize,
    /// Number of threads to use. Defaults to the number of logical cores.
    #[clap(short = 'j', long)]
    threads: Option<usize>,

    /// Dedup across reverse-complements.
    #[clap(long)]
    canonical: bool,

    /// Minimizer length for phrases.
    #[clap(long, default_value = "8")]
    mini_k: usize,
}

#[derive(Default, Clone, Copy, derive_more::AddAssign)]
struct Stats {
    input_bp: usize,
    total_phrases: usize,
    filtered_phrases: usize,
    output_bp: usize,
    unique_phrases: usize,
    output_contigs: usize,
}

fn main() {
    let args = Args::parse();
    let Args {
        input,
        output,
        k,
        w,
        threads,
        ..
    } = Args::parse();

    eprintln!("k: {k}   w: {w}");
    // Open an archive
    let config = DecompressorConfig { verbosity: 0 };
    let decompressor = &mut Decompressor::open(&input.to_string_lossy(), config).unwrap();

    // List available samples
    let samples = decompressor.list_samples();
    let samples: &Vec<_> = &samples
        .into_iter()
        .map(|sample| {
            let contigs = decompressor.list_contigs(&sample).unwrap();
            (sample, contigs)
        })
        .collect();

    println!("Found {} samples", samples.len());

    // TODO: zstd output
    if let Some(output) = &output {
        assert!(
            output.extension().unwrap() == "zst",
            "Output file must have .zst extension"
        );
    }
    let output_path = output.unwrap_or_else(|| input.with_extension("dedup.fa"));
    let buf_writer = BufWriter::with_capacity(1 << 20, std::fs::File::create(output_path).unwrap());
    let writer = &Mutex::new(zstd::Encoder::new(buf_writer, 0).unwrap().auto_finish());
    let global_stats = &Mutex::new(Stats::default());

    let seen: &[_; 256] = &std::array::from_fn(|_i| RwLock::new(FxHashSet::default()));

    // Process the first/reference sample separately.
    let reference = RwLock::new((vec![], FxHashMap::default()));
    process_sample(
        &args,
        decompressor,
        samples,
        seen,
        global_stats,
        writer,
        0,
        &reference,
    );

    let next = &AtomicUsize::new(1);

    std::thread::scope(|scope| {
        let threads = threads.unwrap_or_else(|| num_cpus::get_physical());
        for _t in 0..threads {
            scope.spawn(|| loop {
                let idx = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if idx >= samples.len() {
                    break;
                }
                process_sample(
                    &args,
                    decompressor,
                    samples,
                    seen,
                    global_stats,
                    writer,
                    idx,
                    &reference,
                );
            });
        }
    });
}

fn process_sample(
    args: &Args,
    decompressor: &Decompressor,
    samples: &Vec<(String, Vec<String>)>,
    seen: &[RwLock<FxHashSet<u128>>; 256],
    global_stats: &Mutex<Stats>,
    writer: &Mutex<
        zstd::stream::AutoFinishEncoder<
            '_,
            BufWriter<std::fs::File>,
            Box<dyn FnMut(Result<BufWriter<std::fs::File>, std::io::Error>) + Send>,
        >,
    >,
    idx: usize,
    reference: &RwLock<(Vec<u8>, FxHashMap<u128, usize>)>,
) {
    let Args {
        k,
        mini_k,
        canonical,
        w,
        ..
    } = *args;

    let mut local_stats = Stats::default();
    let mut t_read = Duration::ZERO;
    let mut t_minis = Duration::ZERO;
    let mut t_phrases = Duration::ZERO;
    let mut t_sort = Duration::ZERO;
    let mut t_lookups = Duration::ZERO;
    let mut t_sort2 = Duration::ZERO;
    let mut t_lock = Duration::ZERO;
    let mut t_output = Duration::ZERO;

    let (sample, contigs) = &samples[idx];

    let reference_guard = (idx > 0).then(|| reference.read().unwrap());
    let reference_vec = reference_guard.as_ref().map(|x| &x.0);
    let reference_map = reference_guard.as_ref().map(|x| &x.1);

    let mut build_reference_vec = vec![];
    let mut build_reference_map: FxHashMap<u128, usize> = FxHashMap::default();

    for contig in contigs {
        let i_start = std::time::Instant::now();

        let mut seq = decompressor.get_contig(&sample, &contig).unwrap();
        seq.iter_mut().for_each(|b| *b = b"ACGT"[(*b as usize) % 4]);
        local_stats.input_bp += seq.len();

        let i_read = std::time::Instant::now();
        t_read += i_read - i_start;

        let mut positions = vec![];

        let nthasher = AntiLexHasher::<false>::new(mini_k);
        if canonical {
            simd_minimizers::canonical_minimizers(mini_k, w).run(AsciiSeq(&seq), &mut positions);
        } else {
            // let scheme = simd_minimizers::closed_syncmers(k, w);
            // let scheme = simd_minimizers::minimizers(k, w);
            // Either::Right(simd_minimizers::minimizers(mini_k, w))
            simd_minimizers::minimizers(mini_k, w)
                .hasher(&nthasher)
                .run(AsciiSeq(&seq), &mut positions);
        };

        let i_minis = std::time::Instant::now();
        t_minis += i_minis - i_read;

        let rc_seq: Vec<_> = if canonical {
            let mut rc_seq = vec![];
            rc_seq.extend(seq.iter().rev().map(|bp| 3 - bp));
            rc_seq
        } else {
            vec![]
        };

        let mut phrases = vec![];

        let with_hash = |p, q| {
            let phrase = &seq[p..q];
            let hash = hasher(phrase);
            if canonical {
                let rc_phrase = &rc_seq[seq.len() - q..seq.len() - p];
                let rc_hash = hasher(rc_phrase);
                (p, q, hash + rc_hash)
            } else {
                (p, q, hash)
            }
        };

        // Prefix phrase.
        if positions[0] > 0 {
            phrases.push(with_hash(0, (positions[0] as usize + k).min(seq.len())));
        }
        let mut end_of_seen = 0usize;
        for &[p, q] in positions.array_windows::<2>() {
            // Skip non-forward minimizer pairs.
            if q < p {
                continue;
            }
            let p = p as usize;
            let q = (q as usize + k).min(seq.len());

            local_stats.total_phrases += 1;

            // Add extra buffer so that coming up minimizers don't mess things up.
            if q + w < end_of_seen {
                local_stats.filtered_phrases += 1;
                continue;
            }

            let (p, q, hash) = with_hash(p, q);

            // if q + w < end_of_seen {
            //     eprintln!("p {p} q {q} end_of_seen {end_of_seen}");
            //     assert!(reference_map.unwrap().contains_key(&hash));
            //     continue;
            // }

            if let Some(&pos) = reference_map.and_then(|map| map.get(&hash)) {
                let ref_seq = reference_vec.as_ref().unwrap()[pos..].as_ref();
                let seq = &seq[p..];
                // hash was seen before at given `pos` in `reference_vec`.
                // Linear scan to find the equal range, and skip it.
                assert_eq!(
                    &seq[..q - p],
                    &ref_seq[..q - p],
                    "UNEQUAL RANGES with hashes {} and {} (baseline {hash})",
                    hasher(&seq[..q - p]),
                    hasher(&ref_seq[..q - p])
                );
                let mut i = q - p;
                while i < seq.len().min(ref_seq.len()) && seq[i] == ref_seq[i] {
                    i += 1;
                }
                end_of_seen = p + i;
                // eprintln!(
                //     "Saw hash of {p}..{q} before at {pos}..{}; extend to {end_of_seen}",
                //     pos + q - p
                // );
            }

            phrases.push((p, q, hash));
        }
        // Suffix phrase.
        if (*positions.last().unwrap() as usize) + k < seq.len() {
            phrases.push(with_hash(*positions.last().unwrap() as usize, seq.len()));
        }

        let i_phrases = std::time::Instant::now();
        t_phrases += i_phrases - i_minis;

        let mut order = (0..256).map(|_| vec![]).collect::<Vec<_>>();
        for (i, p) in phrases.iter().enumerate() {
            let part = p.2 as u8;
            order[part as usize].push(i as u32);
        }
        let i_sort = std::time::Instant::now();
        t_sort += i_sort - i_phrases;

        for part in 0..=255 {
            if order[part as usize].is_empty() {
                continue;
            }

            // Read-only filter phase
            {
                let seen = seen[part as usize].read().unwrap();
                for &idx in &order[part as usize] {
                    let (p, _q, hash) = &mut phrases[idx as usize];
                    assert!(*hash as u8 == part);
                    if seen.contains(hash) {
                        *p = usize::MAX;
                    }
                }
            }

            // Try to write missing
            {
                let mut seen = seen[part as usize].write().unwrap();
                for &idx in &order[part as usize] {
                    let (p, _q, hash) = &mut phrases[idx as usize];
                    if *p != usize::MAX {
                        if !seen.insert(*hash) {
                            *p = usize::MAX;
                        }
                    }
                }
            }
        }
        let i_lookups = std::time::Instant::now();
        t_lookups += i_lookups - i_sort;

        phrases.retain(|(p, _q, _hash)| *p != usize::MAX);
        local_stats.unique_phrases += phrases.len();

        let i_sort2 = std::time::Instant::now();
        t_sort2 += i_sort2 - i_lookups;

        // Write new contigs.
        let mut writer = writer.lock().unwrap();
        let i_lock = std::time::Instant::now();
        t_lock += i_lock - i_sort2;

        // Output contigs

        let mut active = 0..0;

        let mut push = |range: Range<usize>| -> usize {
            assert!(range.start >= active.start);
            if range.start <= active.end {
                // End can decrease for non-forward canonical minimizers.
                active.end = active.end.max(range.end);
            } else {
                writer.write_all(b">\n").unwrap();
                writer.write_all(&seq[active.clone()]).unwrap();
                writer.write_all(b"\n").unwrap();

                if idx == 0 {
                    build_reference_vec.extend_from_slice(&seq[active.clone()]);
                    build_reference_vec.push(b'\n');
                }

                local_stats.output_contigs += 1;
                local_stats.output_bp += active.len();

                active = range.clone();
            }
            let ref_pos = build_reference_vec.len() + range.start - active.start;
            ref_pos
        };

        // Emit the phrases.
        for (p, q, hash) in phrases {
            let pos = push(p..q);
            if idx == 0 {
                assert!(build_reference_map.insert(hash, pos).is_none());
            }
        }
        push(usize::MAX..usize::MAX);
        assert!(active.len() == 0);

        let i_output = std::time::Instant::now();
        t_output += i_output - i_lock;
    }

    if idx == 0 {
        eprintln!(
            "Reference vec: {:8.2} Mbp",
            build_reference_vec.len() as f32 / 1e6
        );
        eprintln!(
            "Reference map: {:8.2} M phrases",
            build_reference_map.len() as f32 / 1e6
        );
        *reference.write().unwrap() = (build_reference_vec, build_reference_map);
    }

    eprintln!(
        "push sample {idx:>3} ({:3.1} Gbp {:3} ctg): \
                     read: {:5.2?}s minis: {:5.2?}s phrases: \
                     {:5.2?}s sort: {:5.2?}s lookups: {:5.2?}s sort: {:5.2?}s lock: {:5.2?}s output: {:5.2?}s",
        (local_stats.input_bp as f32) / 1e9,
        contigs.len(),
        t_read.as_secs_f32(),
        t_minis.as_secs_f32(),
        t_phrases.as_secs_f32(),
        t_sort.as_secs_f32(),
        t_lookups.as_secs_f32(),
        t_sort2.as_secs_f32(),
        t_lock.as_secs_f32(),
        t_output.as_secs_f32()
    );

    let mut global_stats = global_stats.lock().unwrap();
    *global_stats += local_stats;

    eprintln!(
        "  new bp:            {:>8.3} Mbp ({:3.1} bp/contig)",
        local_stats.output_bp as f32 / 1e6,
        local_stats.output_bp as f32 / local_stats.output_contigs as f32
    );
    eprintln!(
        "  new phrases:       {:>8.3} M   ({:3.1} /contig; {:4.1}%; {:4.1}% filtered away)",
        local_stats.unique_phrases as f32 / 1e6,
        local_stats.unique_phrases as f32 / local_stats.output_contigs as f32,
        local_stats.unique_phrases as f32 / local_stats.total_phrases as f32 * 100.0,
        local_stats.filtered_phrases as f32 / local_stats.total_phrases as f32 * 100.0,
    );
    eprintln!(
        "  unique phrases:    {:>8.3} M   ({:3.1}%)",
        global_stats.unique_phrases as f32 / 1e6,
        100.0 * global_stats.unique_phrases as f32 / global_stats.total_phrases as f32
    );

    eprintln!(
        "  num_contigs:       {:>8.3} M",
        global_stats.output_contigs as f32 / 1e6
    );
    eprintln!(
        "  output_bp:         {:>8.3} Gbp ({:3.1}%)",
        global_stats.output_bp as f32 / 1e9,
        100.0 * global_stats.output_bp as f32 / global_stats.input_bp as f32
    );
    eprintln!();
}
