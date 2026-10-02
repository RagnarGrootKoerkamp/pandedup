//! Run GGCAT on a FASTA file, optionally compressed with zstd.
#[path = "../src/timing.rs"]
mod timing;

use clap::{Parser, ValueEnum};
use ggcat_api::{
    DnaSequence, DnaSequencesFileType, DynamicSequencesStream, ExtraElaboration, GGCATConfig,
    GGCATInstance, GeneralSequenceBlockData, MessageLevel, SequenceInfo,
};
use std::{
    io::{IsTerminal, Write},
    path::{Path, PathBuf},
    sync::Arc,
};
use tracing::info;

#[derive(Parser)]
struct Args {
    /// Input .fa or .fa.zst file.
    input: PathBuf,
    /// Output FASTA file. Defaults to INPUT-kK-MODE.fa.
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// K-mer size.
    #[arg(short, long, default_value_t = 64)]
    k: usize,
    /// Number of threads. Defaults to the number of physical cores.
    #[arg(short = 'j', long)]
    threads: Option<usize>,
    /// GGCAT extra elaboration mode.
    #[arg(long, value_enum, default_value_t = Mode::Simplitigs)]
    mode: Mode,
}

#[derive(Clone, Copy, ValueEnum)]
enum Mode {
    None,
    UnitigLinks,
    GreedyMatchtigs,
    Eulertigs,
    Pathtigs,
    Simplitigs,
    FastEulertigs,
}

impl Mode {
    fn elaboration(self) -> ExtraElaboration {
        match self {
            Self::None => ExtraElaboration::None,
            Self::UnitigLinks => ExtraElaboration::UnitigLinks,
            Self::GreedyMatchtigs => ExtraElaboration::GreedyMatchtigs,
            Self::Eulertigs => ExtraElaboration::Eulertigs,
            Self::Pathtigs => ExtraElaboration::Pathtigs,
            Self::Simplitigs => ExtraElaboration::FastSimplitigs,
            Self::FastEulertigs => ExtraElaboration::FastEulertigs,
        }
    }
}

fn main() {
    tracing_subscriber::fmt()
        .compact()
        .with_target(false)
        .with_writer(std::io::stderr)
        .with_ansi(std::io::stderr().is_terminal())
        .init();
    let args = Args::parse();
    let timing = timing::StageTiming::start();
    let name = args.input.to_string_lossy();
    assert!(
        name.ends_with(".fa") || name.ends_with(".fa.zst"),
        "input must be .fa or .fa.zst"
    );
    let mut reader = needletail::parse_fastx_file(&args.input).expect("failed to read input FASTA");
    let mut sequences = Vec::new();
    while let Some(record) = reader.next() {
        sequences.push(record.expect("invalid FASTA record").seq().into_owned());
    }
    let input_bases = sequences.iter().map(Vec::len).sum();
    pandedup::log_file_stats("Read", &args.input, Some((sequences.len(), input_bases))).unwrap();
    let output = run(&args, sequences);
    pandedup::log_file_stats("Wrote", &output, None).unwrap();
    info!(
        "ggcat: {}, input {} bytes, output {} bytes ({})",
        timing.finish(),
        std::fs::metadata(&args.input).unwrap().len(),
        std::fs::metadata(&output).unwrap().len(),
        output.display()
    );
}

fn default_output_path(input: &Path, k: usize, mode: &str) -> PathBuf {
    let name = input.file_name().unwrap().to_string_lossy();
    let stem = name
        .strip_suffix(".fa.zst")
        .or_else(|| name.strip_suffix(".fa"))
        .expect("input must be .fa or .fa.zst");
    input.with_file_name(format!("{stem}-k{k}-{mode}.fa"))
}

struct InMemorySequences {
    blocks: Vec<Vec<Vec<u8>>>,
    bases: Vec<u64>,
}

impl InMemorySequences {
    fn new(sequences: Vec<Vec<u8>>) -> Self {
        let mut blocks = vec![Vec::new()];
        let mut bases = vec![0];
        for sequence in sequences {
            if bases.last().copied().unwrap() >= 64 * 1024 * 1024 {
                blocks.push(Vec::new());
                bases.push(0);
            }
            *bases.last_mut().unwrap() += sequence.len() as u64;
            blocks.last_mut().unwrap().push(sequence);
        }
        Self { blocks, bases }
    }
}

impl DynamicSequencesStream for InMemorySequences {
    fn read_block(
        &self,
        block: usize,
        _copy_ident_data: bool,
        partial_read_copyback: Option<usize>,
        callback: &mut dyn FnMut(DnaSequence<'_, &[u8]>, SequenceInfo),
    ) {
        for sequence in &self.blocks[block] {
            let Some(copyback) = partial_read_copyback else {
                callback(
                    DnaSequence {
                        ident_data: b">pandedup",
                        seq: sequence,
                        format: DnaSequencesFileType::FASTA,
                    },
                    SequenceInfo { color: None },
                );
                continue;
            };
            // Match GGCAT's FASTA reader: split long sequences with an
            // overlapping suffix and mark intermediate pieces as FASTQ.
            let chunk_size = (1 << 20).max(copyback.saturating_mul(2));
            let mut start = 0;
            while start < sequence.len() {
                let end = start.saturating_add(chunk_size).min(sequence.len());
                callback(
                    DnaSequence {
                        ident_data: b">pandedup",
                        seq: &sequence[start..end],
                        format: if end == sequence.len() {
                            DnaSequencesFileType::FASTA
                        } else {
                            DnaSequencesFileType::FASTQ
                        },
                    },
                    SequenceInfo { color: None },
                );
                if end == sequence.len() {
                    break;
                }
                start = end - copyback;
            }
        }
    }

    fn estimated_base_count(&self, block: usize) -> u64 {
        self.bases[block]
    }
}

fn run(args: &Args, sequences: Vec<Vec<u8>>) -> PathBuf {
    let (short_sequences, sequences): (Vec<_>, Vec<_>) = sequences
        .into_iter()
        .partition(|sequence| sequence.len() < args.k);
    let input = Arc::new(InMemorySequences::new(sequences));
    let input_bases: u64 = input.bases.iter().sum();
    let input_contigs: usize = input.blocks.iter().map(Vec::len).sum();
    info!(
        "ggcat input: {input_contigs} in-memory contigs, {input_bases} bases (from {})",
        args.input.display()
    );

    let threads = args.threads.unwrap_or_else(num_cpus::get_physical);
    let temp_dir = args
        .input
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let ggcat = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(temp_dir.to_path_buf()),
        memory: 50.0,
        prefer_memory: true,
        total_threads_count: threads,
        intermediate_compression_level: None,
        stats_file: None,
        messages_callback: Some(|level, message| match level {
            MessageLevel::Info => {}
            MessageLevel::Warning | MessageLevel::Error => eprintln!("ggcat: {message}"),
            MessageLevel::UnrecoverableError => panic!("ggcat: {message}"),
        }),
    })
    .expect("failed to initialize GGCAT");

    let label = args.mode.to_possible_value().unwrap().get_name().to_owned();
    let elaboration = args.mode.elaboration();
    let output = args
        .output
        .clone()
        .unwrap_or_else(|| default_output_path(&args.input, args.k, &label));
    let blocks = (0..input.blocks.len())
        .map(|index| {
            GeneralSequenceBlockData::Dynamic((
                input.clone() as Arc<dyn DynamicSequencesStream>,
                index,
            ))
        })
        .collect();
    let output = if input_contigs == 0 {
        std::fs::File::create(&output).expect("failed to create GGCAT output");
        output
    } else {
        ggcat
            .build_graph(
                blocks,
                output,
                None,
                args.k,
                threads,
                /* forward_only */
                false,
                None,
                false,
                1,
                elaboration,
                None,
            )
            .unwrap_or_else(|error| panic!("GGCAT {label} failed: {error:#}"))
    };
    if !short_sequences.is_empty() {
        let mut file = std::fs::OpenOptions::new()
            .append(true)
            .open(&output)
            .expect("failed to open GGCAT output for short contigs");
        for (index, sequence) in short_sequences.iter().enumerate() {
            writeln!(file, ">pandedup-short-{index}").unwrap();
            file.write_all(sequence).unwrap();
            file.write_all(b"\n").unwrap();
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_output_uses_input_stem_k_and_mode() {
        for input in ["x.fa", "x.fa.zst", "x.part.fa.zst"] {
            let expected = if input.starts_with("x.part") {
                "x.part-k31-greedy-matchtigs.fa"
            } else {
                "x-k31-greedy-matchtigs.fa"
            };
            assert_eq!(
                default_output_path(Path::new(input), 31, "greedy-matchtigs"),
                PathBuf::from(expected)
            );
        }
    }

    #[test]
    fn in_memory_stream_preserves_overlap_between_chunks() {
        let original = vec![b'A'; (1 << 20) + 10];
        let stream = InMemorySequences::new(vec![original.clone()]);
        let mut chunks = Vec::new();
        stream.read_block(0, false, Some(7), &mut |sequence, _| {
            chunks.push((
                sequence.seq.to_vec(),
                matches!(sequence.format, DnaSequencesFileType::FASTQ),
            ));
        });
        assert_eq!(chunks.len(), 2);
        assert!(chunks[0].1);
        assert!(!chunks[1].1);
        assert_eq!(&chunks[0].0[chunks[0].0.len() - 7..], &chunks[1].0[..7]);
        let mut reconstructed = chunks.remove(0).0;
        reconstructed.extend_from_slice(&chunks[0].0[7..]);
        assert_eq!(reconstructed, original);
    }
}
