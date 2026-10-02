//! Compute a masked superstring using kmercamel's bigreedy algorithm.
//!
//! See: https://github.com/OndrejSladky/kmercamel, and the paper:
//!
//! Ondřej Sladký, Pavel Veselý, and Karel Břinda:
//! Masked superstrings as a unified framework for textual k-mer set representations.
//! bioRxiv 2023.02.01.526717, 2023. https://doi.org/10.1101/2023.02.01.526717

use seq_hash::packed_seq::{self, Seq};
use tracing::info;
use voracious_radix_sort::{RadixSort, Radixable};

#[derive(Clone, Copy)]
struct SortEntry {
    key: u128,
    index: usize,
    reverse: usize,
}

impl PartialEq for SortEntry {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}

impl PartialOrd for SortEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.key.cmp(&other.key))
    }
}

impl Radixable<u128> for SortEntry {
    type Key = u128;

    fn key(&self) -> Self::Key {
        self.key
    }
}

struct Contig {
    /// The head of the contig, in fwd and rc direction
    kmer_in: [u128; 2],
    /// The tail of the contig, in fwd and rc direction
    kmer_out: [u128; 2],
    /// Neighbors at the forward head (0) and forward tail (1).
    nbs: [Option<Link>; 2],
}

struct Link {
    index: usize,
    reverse: bool,
    overlap: u8,
}

/// Input: unitigs, simplitigs=pathtigs, eulertigs (or (greedy) matchtigs).
///
/// Output: a superstring that contains all the input contigs, with newly
/// introduced "stitching" kmers masked out by lowercasing them.
///
/// If the input contains duplicate kmers, in case of (greedy) matchtigs, those
/// will be preserved in the output.
pub fn masked_superstring(k: usize, contigs: &Vec<Vec<u8>>) -> Vec<u8> {
    info!("masked_superstring: k={}, contigs={}", k, contigs.len());
    assert!((1..=64).contains(&k));

    // Each contig is stored as a fwd and rc kmer for the head, and one for the tail.
    // [[fwd in, fwd out], [rc in, rc out]]
    let mut ends: Vec<Contig> = vec![];
    for c in contigs {
        // TODO handle short contigs
        if c.len() < k {
            continue;
        }
        let seq = packed_seq::AsciiSeq(c);
        let head_fw = seq.read_kmer_u128(k, 0);
        let head_rc = seq.read_revcomp_kmer_u128(k, 0);
        let tail_fw = seq.read_kmer_u128(k, c.len() - k);
        let tail_rc = seq.read_revcomp_kmer_u128(k, c.len() - k);
        ends.push(Contig {
            kmer_in: [head_fw, tail_rc],
            kmer_out: [tail_fw, head_rc],
            nbs: [None, None],
        });
    }

    let mut total_merged = 0;
    let mut total_len = contigs.iter().map(|c| c.len()).sum::<usize>();

    for overlap in (0..=k - 1).rev() {
        info!("overlap: {}", overlap);
        let head_mask = (1u128 << (2 * overlap)) - 1;
        let tail_shift = 2 * (k - overlap);
        let mut heads = Vec::with_capacity(2 * ends.len());
        let mut tails = Vec::with_capacity(2 * ends.len());
        for (index, contig) in ends.iter().enumerate() {
            for reverse in 0..2 {
                // An RC head uses the physical forward tail, and vice versa.
                if contig.nbs[reverse].is_none() {
                    heads.push(SortEntry {
                        key: contig.kmer_in[reverse] & head_mask,
                        index,
                        reverse,
                    });
                }
                if contig.nbs[1 - reverse].is_none() {
                    tails.push(SortEntry {
                        key: contig.kmer_out[reverse]
                            .checked_shr(tail_shift as u32)
                            .unwrap_or(0),
                        index,
                        reverse,
                    });
                }
            }
        }
        let threads = std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1);
        info!("sorting {} heads..", heads.len());
        heads.voracious_mt_sort(threads);
        info!("sorting {} tails..", tails.len());
        tails.voracious_mt_sort(threads);

        // mergesort the two lists; make a connection whenever possible.
        info!("merging..");
        let mut i = 0;
        let mut j = 0;
        let mut merged = 0;
        while i < heads.len() && j < tails.len() {
            let head = heads[i];
            let tail = tails[j];
            if head.key == tail.key {
                // Skip if one of the free ends was already filled in the rc direction.
                if ends[head.index].nbs[head.reverse].is_some() {
                    i += 1;
                    continue;
                }
                if ends[tail.index].nbs[1 - tail.reverse].is_some() {
                    j += 1;
                    continue;
                }
                // Skip self-loops from an end into itself.
                if head.index == tail.index && head.reverse == 1 - tail.reverse {
                    i += 1;
                    continue;
                }

                let head_link = Link {
                    index: head.index,
                    reverse: head.reverse != 0,
                    overlap: overlap as u8,
                };
                let tail_link = Link {
                    index: tail.index,
                    reverse: tail.reverse != 0,
                    overlap: overlap as u8,
                };
                ends[head.index].nbs[head.reverse] = Some(tail_link);
                ends[tail.index].nbs[1 - tail.reverse] = Some(head_link);
                i += 1;
                j += 1;
                merged += 1;
            } else if head.key < tail.key {
                i += 1;
            } else {
                j += 1;
            }
        }
        total_merged += merged;
        total_len -= merged * overlap;
        info!(
            "overlap {overlap} merged {:>9} total merged {total_merged:>9} remaining {:>9} total len {total_len:>11}",
            merged,
            ends.len() - total_merged
        );
    }
    assert_eq!(total_merged, ends.len());
    // TODO break cycles

    todo!()
}
