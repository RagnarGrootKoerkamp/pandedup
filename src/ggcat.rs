use crate::Args;
use ggcat_api::{
    DnaSequence, DnaSequencesFileType, DynamicSequencesStream, ExtraElaboration, GGCATConfig,
    GGCATInstance, GeneralSequenceBlockData, MessageLevel, SequenceInfo,
};
use std::{path::Path, sync::Arc, time::Instant};

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

pub fn run(args: &Args, dedup_path: &Path, sequences: Vec<Vec<u8>>) {
    let input = Arc::new(InMemorySequences::new(sequences));
    let input_bases: u64 = input.bases.iter().sum();
    let input_contigs: usize = input.blocks.iter().map(Vec::len).sum();
    println!(
        "ggcat input: {input_contigs} in-memory contigs, {input_bases} bases (from {})",
        dedup_path.display()
    );

    let threads = args.threads.unwrap_or_else(num_cpus::get_physical);
    let temp_dir = dedup_path
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

    for (label, elaboration) in [
        ("none", ExtraElaboration::None),
        ("unitig-links", ExtraElaboration::UnitigLinks),
        ("simplitigs", ExtraElaboration::FastSimplitigs),
        ("eulertigs", ExtraElaboration::FastEulertigs),
        ("greedy-matchtigs", ExtraElaboration::GreedyMatchtigs),
    ] {
        let output = dedup_path
            .with_extension("")
            .with_extension(format!("ggcat-{label}.fa"));
        let blocks = (0..input.blocks.len())
            .map(|index| {
                GeneralSequenceBlockData::Dynamic((
                    input.clone() as Arc<dyn DynamicSequencesStream>,
                    index,
                ))
            })
            .collect();
        let start = Instant::now();
        let output = ggcat
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
            .unwrap_or_else(|error| panic!("GGCAT {label} failed: {error:#}"));
        println!(
            "ggcat {label}: {:.2?}, {} bytes ({})",
            start.elapsed(),
            std::fs::metadata(&output).unwrap().len(),
            output.display()
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
