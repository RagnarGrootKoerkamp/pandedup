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
    /// The previous contig and overlap
    prev: Option<(usize, u8)>,
    /// The next contig and overlap
    next: Option<(usize, u8)>,
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
    assert!(k <= 64);

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
        let head_rc = seq.read_kmer_u128(k, 0);
        let tail_fw = seq.read_kmer_u128(k, c.len() - k);
        let tail_rc = seq.read_kmer_u128(k, c.len() - k);
        ends.push(Contig {
            kmer_in: [head_fw, tail_rc],
            kmer_out: [tail_fw, head_rc],
            prev: None,
            next: None,
        });
    }

    let mut merged = 0;

    for overlap in (0..=k - 1).rev() {
        info!("overlap: {}", overlap);
        let head_mask = (1u128 << (2 * overlap)) - 1;
        let tail_shift = 2 * (k - overlap);
        // TODO handle rc
        let mut heads = ends
            .iter()
            .enumerate()
            .filter(|(_, e)| e.prev.is_none())
            .map(|(i, e)| SortEntry {
                key: e.kmer_in[0] & head_mask,
                index: i,
            })
            .collect::<Vec<_>>();
        let mut tails = ends
            .iter()
            .enumerate()
            .filter(|(_, e)| e.next.is_none())
            .map(|(i, e)| SortEntry {
                key: e.kmer_out[0] >> tail_shift,
                index: i,
            })
            .collect::<Vec<_>>();
        let threads = std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1);
        info!("sorting head..");
        heads.voracious_mt_sort(threads);
        info!("sorting tail..");
        tails.voracious_mt_sort(threads);

        // mergesort the two lists; make a connection whenever possible.
        info!("merging..");
        let mut i = 0;
        let mut j = 0;
        while i < heads.len() && j < tails.len() {
            let head_kmer = heads[i].key;
            let head_idx = heads[i].index;
            let tail_kmer = tails[j].key;
            let tail_idx = tails[j].index;
            if head_kmer == tail_kmer {
                // connect the two contigs
                ends[head_idx].prev = Some((tail_idx, overlap as u8));
                ends[tail_idx].next = Some((head_idx, overlap as u8));
                i += 1;
                j += 1;
                merged += 1;
            } else if head_kmer < tail_kmer {
                i += 1;
            } else {
                j += 1;
            }
        }
        info!("merged {:>9} remaining {}", merged, ends.len() - merged);
    }
    // TODO break cycles

    todo!()
}
