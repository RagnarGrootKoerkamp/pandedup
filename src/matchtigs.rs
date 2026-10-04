//! Greedy matching of unitig ends through the compacted de Bruijn graph.

use crate::{log_file_stats, mss::MssKey, timing::StageTiming};
use packed_seq::{PackedSeqVec, SeqVec, complement_char};
use std::fmt::{self, Display, Formatter};
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use tracing::info;
use voracious_radix_sort::{RadixSort, Radixable};

const DEAD: u32 = u32::MAX;

/// Decimal SI units with three significant digits for counts of at least 1k.
struct Compact(u128);

fn compact(value: impl TryInto<u64>) -> Compact {
    Compact(value.try_into().ok().expect("count exceeds u64") as u128)
}

impl Display for Compact {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        const SUFFIXES: [&str; 7] = ["", "k", "M", "G", "T", "P", "E"];
        let value = self.0;
        let mut unit = 0;
        let mut divisor = 1u128;
        while unit + 1 < SUFFIXES.len() && value >= divisor * 1_000 {
            unit += 1;
            divisor *= 1_000;
        }
        if unit == 0 {
            return write!(f, "{value}");
        }
        loop {
            let whole = value / divisor;
            let mut decimals = if whole < 10 {
                2
            } else if whole < 100 {
                1
            } else {
                0
            };
            loop {
                let factor = 10u128.pow(decimals);
                let rounded = (value * factor + divisor / 2) / divisor;
                if rounded >= 1_000 {
                    if decimals > 0 {
                        decimals -= 1;
                        continue;
                    }
                    if unit + 1 < SUFFIXES.len() {
                        unit += 1;
                        divisor *= 1_000;
                        break;
                    }
                }
                return match decimals {
                    2 => write!(
                        f,
                        "{}.{:02}{}",
                        rounded / 100,
                        rounded % 100,
                        SUFFIXES[unit]
                    ),
                    1 => write!(f, "{}.{}{}", rounded / 10, rounded % 10, SUFFIXES[unit]),
                    _ => write!(f, "{}{}", rounded, SUFFIXES[unit]),
                };
            }
        }
    }
}

#[derive(Clone, Copy, Default)]
struct VisitStats {
    direct_overlap: u64,
    nearest_receiver: u64,
    connect: u64,
    remove_source: u64,
}

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
#[doc(hidden)]
pub struct End<K> {
    kmer: K,
    id: u32,
}

const _: () = {
    assert!(std::mem::size_of::<End<u64>>() == 12);
    assert!(std::mem::size_of::<End<u128>>() == 20);
};

impl Radixable<u64> for End<u64> {
    type Key = u64;
    fn key(&self) -> Self::Key {
        self.kmer
    }
}

impl Radixable<u128> for End<u128> {
    type Key = u128;
    fn key(&self) -> Self::Key {
        self.kmer
    }
}

impl<K: MssKey> PartialEq for End<K> {
    fn eq(&self, other: &Self) -> bool {
        let left = self.kmer;
        let right = other.kmer;
        left == right
    }
}

impl<K: MssKey> PartialOrd for End<K> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let left = self.kmer;
        let right = other.kmer;
        Some(left.cmp(&right))
    }
}

#[derive(Clone, Copy, PartialEq, PartialOrd)]
struct Edge {
    from: u32,
    to: u32,
}

impl Radixable<u64> for Edge {
    type Key = u64;
    fn key(&self) -> Self::Key {
        (u64::from(self.from) << 32) | u64::from(self.to)
    }
}

/// CSR adjacency indexed by an oriented outgoing tail. An edge's destination
/// is the oriented unitig entered at its head; its low bit is the orientation.
struct Graph {
    edges: Vec<u32>,
    offsets: Vec<u32>,
}

impl Graph {
    fn outgoing(&self, id: u32) -> &[u32] {
        &self.edges[self.offsets[id as usize] as usize..self.offsets[id as usize + 1] as usize]
    }
}

fn graph<K: MssKey>(k: usize, seq: &PackedSeqVec, ranges: &[Range<usize>]) -> Graph
where
    End<K>: Radixable<K, Key = K>,
{
    let overlap = k - 1;
    let mut heads = Vec::with_capacity(2 * ranges.len());
    let mut tails = Vec::with_capacity(2 * ranges.len());
    info!("Building graph of {} unitigs..", compact(ranges.len()));
    for (index, range) in ranges.iter().enumerate() {
        let id = (index as u32) * 2;
        let (head, tail, rc_head, rc_tail): (K, K, K, K) = (
            K::read_kmer(seq, overlap, range.start),
            K::read_kmer(seq, overlap, range.end - overlap),
            K::read_revcomp_kmer(seq, overlap, range.start),
            K::read_revcomp_kmer(seq, overlap, range.end - overlap),
        );
        heads.push(End { kmer: head, id });
        heads.push(End {
            kmer: rc_tail,
            id: id + 1,
        });
        tails.push(End { kmer: tail, id });
        tails.push(End {
            kmer: rc_head,
            id: id + 1,
        });
    }
    info!("Sorting heads");
    heads.voracious_mt_sort(rayon::current_num_threads());
    info!("Sorting tails");
    tails.voracious_mt_sort(rayon::current_num_threads());

    info!("Matching unitig ends");
    let mut pairs = Vec::new();
    let (mut h, mut t) = (0, 0);
    while h < heads.len() && t < tails.len() {
        // Copy keys out of packed entries before borrowing them for comparison.
        let head_key = heads[h].kmer;
        let tail_key = tails[t].kmer;
        match head_key.cmp(&tail_key) {
            std::cmp::Ordering::Less => h += 1,
            std::cmp::Ordering::Greater => t += 1,
            std::cmp::Ordering::Equal => {
                let key = heads[h].kmer;
                let mut he = h;
                let mut te = t;
                while he < heads.len() {
                    let next = heads[he].kmer;
                    if next != key {
                        break;
                    }
                    he += 1;
                }
                while te < tails.len() {
                    let next = tails[te].kmer;
                    if next != key {
                        break;
                    }
                    te += 1;
                }
                for tail in &tails[t..te] {
                    for head in &heads[h..he] {
                        // The reverse orientation of the same endpoint does not
                        // traverse a unitig or introduce a new graph edge.
                        if tail.id != (head.id ^ 1) {
                            pairs.push(Edge {
                                from: tail.id,
                                to: head.id,
                            });
                        }
                    }
                }
                h = he;
                t = te;
            }
        }
    }
    info!("Sorting {} edges", compact(pairs.len()));
    pairs.voracious_mt_sort(rayon::current_num_threads());
    // let old_len = pairs.len();
    // info!("Dedup edges");
    // pairs.dedup();
    // info!("Dedup edges => {} left", pairs.len());
    // assert_eq!(old_len, pairs.len(), "duplicate graph edges");
    assert!(pairs.len() < u32::MAX as usize, "too many graph edges");
    info!("Building CSR adjacency");
    let mut offsets = vec![0; 2 * ranges.len() + 1];
    info!(
        "Filling offsets: {}B",
        compact(std::mem::size_of_val(offsets.as_slice()))
    );
    for edge in &pairs {
        offsets[edge.from as usize + 1] += 1;
    }
    for i in 1..offsets.len() {
        offsets[i] += offsets[i - 1];
    }
    info!("Mapping pairs");
    Graph {
        edges: pairs.into_iter().map(|edge| edge.to).collect(),
        offsets,
    }
}

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
struct Link {
    target: u32,
    overlap: u8,
}

const _: () = assert!(std::mem::size_of::<Link>() == 5);

impl Link {
    const EMPTY: Self = Self {
        target: DEAD,
        overlap: 0,
    };

    fn is_empty(self) -> bool {
        self.target == DEAD
    }
}

const UNREACHABLE: u8 = u8::MAX;

fn traversal_weight(ranges: &[Range<usize>], node: u32, k: usize) -> usize {
    ranges[(node / 2) as usize].len() - (k - 1)
}

/// Distance to the nearest still-free outgoing end for each oriented tail.
struct DistanceField {
    distance: Vec<u8>,
    worklist: Vec<(u32, u8)>,
}

impl DistanceField {
    fn new(nodes: usize) -> Self {
        Self {
            distance: vec![0; nodes],
            worklist: Vec::new(),
        }
    }

    fn replacement_distance(
        &self,
        node: u32,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        visited: &mut u64,
    ) -> u8 {
        let weight = traversal_weight(ranges, node, k);
        let predecessors = graph.outgoing(node ^ 1);
        *visited += 1 + predecessors.len() as u64;
        let mut best = k;
        for &reverse_predecessor in predecessors {
            let predecessor = reverse_predecessor ^ 1;
            best = best.min(self.distance[predecessor as usize] as usize + weight);
        }
        if best == k { UNREACHABLE } else { best as u8 }
    }

    /// Propagate increases after deleting a source. A neighbor only needs
    /// repair if its current label was tight through the changed node.
    fn remove_source(
        &mut self,
        source: u32,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        visited: &mut u64,
    ) {
        assert_eq!(
            self.distance[source as usize], 0,
            "active source has nonzero distance"
        );
        self.worklist.clear();
        self.distance[source as usize] =
            self.replacement_distance(source, graph, ranges, k, visited);
        debug_assert!(self.distance[source as usize] > 0);
        self.worklist.push((source, 0));
        while let Some((node, old_distance)) = self.worklist.pop() {
            let next_nodes = graph.outgoing(node);
            *visited += 1 + next_nodes.len() as u64;
            for &next in next_nodes {
                let old_next = self.distance[next as usize];
                if old_next != UNREACHABLE
                    && old_distance as usize + traversal_weight(ranges, next, k)
                        == old_next as usize
                {
                    let new_next = self.replacement_distance(next, graph, ranges, k, visited);
                    debug_assert!(new_next >= old_next);
                    if new_next > old_next {
                        self.distance[next as usize] = new_next;
                        self.worklist.push((next, old_next));
                    }
                }
            }
        }
    }

    /// Follow tight incoming edges to recover any source other than forbidden.
    /// Distances strictly decrease, so recursion is bounded by k.
    fn source(
        &self,
        node: u32,
        forbidden: u32,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        visited: &mut u64,
    ) -> Option<u32> {
        *visited += 1;
        let distance = self.distance[node as usize];
        if distance == UNREACHABLE {
            return None;
        }
        if distance == 0 {
            return (node != forbidden).then_some(node);
        }
        let weight = traversal_weight(ranges, node, k);
        if weight > distance as usize {
            return None;
        }
        let predecessors = graph.outgoing(node ^ 1);
        *visited += predecessors.len() as u64;
        for &reverse_predecessor in predecessors {
            let predecessor = reverse_predecessor ^ 1;
            let prior = self.distance[predecessor as usize];
            if prior != UNREACHABLE && prior as usize + weight == distance as usize {
                if let Some(source) = self.source(predecessor, forbidden, graph, ranges, k, visited)
                {
                    return Some(source);
                }
            }
        }
        None
    }
}

/// Reusable scratch when all nearest tight paths reach the receiver's own
/// physical end, which cannot be linked to itself.
struct ReverseSearch {
    distance: Vec<u8>,
    touched: Vec<u32>,
}

impl ReverseSearch {
    fn new(nodes: usize) -> Self {
        Self {
            distance: vec![UNREACHABLE; nodes],
            touched: Vec::new(),
        }
    }

    fn push(&mut self, node: u32, distance: usize, k: usize) {
        if distance >= k || distance >= self.distance[node as usize] as usize {
            return;
        }
        if self.distance[node as usize] == UNREACHABLE {
            self.touched.push(node);
        }
        self.distance[node as usize] = distance as u8;
    }

    fn nearest(
        &mut self,
        target: u32,
        limit: usize,
        field: &DistanceField,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        visited: &mut u64,
    ) -> Option<(u32, usize)> {
        let predecessors = graph.outgoing(target ^ 1);
        *visited += predecessors.len() as u64;
        for &reverse_predecessor in predecessors {
            self.push(reverse_predecessor ^ 1, 0, k);
        }
        let mut result = None;
        let mut best_distance = limit;
        for distance in 0..limit.min(k) {
            if distance >= best_distance {
                break;
            }
            let mut cursor = 0;
            while cursor < self.touched.len() {
                let node = self.touched[cursor];
                cursor += 1;
                *visited += 1;
                if self.distance[node as usize] as usize != distance {
                    continue;
                }
                let combined = distance + field.distance[node as usize] as usize;
                if combined < best_distance {
                    if let Some(source) = field.source(node, target ^ 1, graph, ranges, k, visited)
                    {
                        best_distance = combined;
                        result = Some((source, combined));
                    }
                }
                let next_distance = distance.saturating_add(traversal_weight(ranges, node, k));
                if next_distance < best_distance && next_distance < k {
                    let predecessors = graph.outgoing(node ^ 1);
                    *visited += predecessors.len() as u64;
                    for &reverse_predecessor in predecessors {
                        self.push(reverse_predecessor ^ 1, next_distance, k);
                    }
                }
            }
        }
        *visited += self.touched.len() as u64;
        for &node in &self.touched {
            self.distance[node as usize] = UNREACHABLE;
        }
        self.touched.clear();
        result
    }
}

fn tail_slot(id: u32) -> usize {
    (id ^ 1) as usize
}
fn head_slot(id: u32) -> usize {
    id as usize
}

fn connect(
    links: &mut [Link],
    from: u32,
    to: u32,
    overlap: usize,
    total_bases: &mut usize,
    visited: &mut u64,
) {
    debug_assert!(overlap <= u8::MAX as usize);
    links[tail_slot(from)] = Link {
        target: to,
        overlap: overlap as u8,
    };
    links[head_slot(to)] = Link {
        target: from ^ 1,
        overlap: overlap as u8,
    };
    *total_bases -= overlap;
    *visited += 2;
}

/// Find the cheapest still-free source for one receiving head. Distances in the
/// field are to outgoing tails, so entering the receiver costs nothing here.
fn nearest_receiver(
    target: u32,
    field: &DistanceField,
    reverse_search: &mut ReverseSearch,
    graph: &Graph,
    ranges: &[Range<usize>],
    k: usize,
    visited: &mut u64,
) -> Option<(u32, usize)> {
    let forbidden = target ^ 1;
    let mut best = None;
    let mut forbidden_distance = usize::MAX;
    let predecessors = graph.outgoing(target ^ 1);
    *visited += predecessors.len() as u64;
    for &reverse_predecessor in predecessors {
        let predecessor = reverse_predecessor ^ 1;
        let index = predecessor as usize;
        let distance = field.distance[index] as usize;
        if distance >= k {
            continue;
        }
        if let Some(source) = field.source(predecessor, forbidden, graph, ranges, k, visited) {
            if best.is_none_or(|(_, old_distance)| distance < old_distance) {
                best = Some((source, distance));
            }
        } else {
            forbidden_distance = forbidden_distance.min(distance);
        }
    }
    let limit = best.map_or(k, |(_, distance)| distance);
    if forbidden_distance <= limit {
        reverse_search
            .nearest(target, limit, field, graph, ranges, k, visited)
            .or(best)
    } else {
        best
    }
}

/// Greedy distance-ordered matching with one shortest-path distance per graph
/// node. A receiver's scheduled distance is a lower bound: deleting sources
/// can only increase it. Repair the field immediately after every link, so
/// each receiver needs to be examined at most once at a given distance.
fn match_ends(
    k: usize,
    ranges: &[Range<usize>],
    graph: &Graph,
    initial_bases: usize,
) -> (Vec<Link>, usize) {
    let nodes = 2 * ranges.len();
    let mut links = vec![Link::EMPTY; nodes];
    let mut field = DistanceField::new(nodes);
    let mut reverse_search = ReverseSearch::new(nodes);
    let mut receiver_distance = vec![0u8; nodes];
    let mut linked = 0usize;
    let mut remaining_receivers = nodes;
    let mut visits = VisitStats::default();
    let mut estimated_bases = initial_bases;
    info!(
        "Matching {} unitigs, {} ends, k {}: two ends per unitig; overlap = k minus one minus distance; examined = receiver checks; links made = connections; receiving ends remain = unlinked ends scheduled below k; free ends = initial ends minus twice cumulative links",
        compact(ranges.len()),
        compact(nodes),
        compact(k),
    );
    info!(
        "Before distance {}: {} receiving ends examined, {} links made, {} receiving ends remain, {} free ends, {} estimated bases",
        compact(0u64),
        compact(0u64),
        compact(0u64),
        compact(remaining_receivers),
        compact(nodes),
        compact(initial_bases),
    );
    info!(
        "Node visits count repeated inspections: direct_overlap scans graph neighbors at the initial distance; nearest_receiver includes its reverse search; connect touches two link slots; remove_source scans changed nodes and edges used to recompute distances"
    );

    // An edge in the unitig graph is a direct (k-1)-character overlap. Match
    // those edges without querying shortest-path labels or a reverse search.
    let mut direct_examined = 0usize;
    let mut direct_made = 0usize;
    for target in 0..nodes as u32 {
        if !links[head_slot(target)].is_empty() {
            continue;
        }
        direct_examined += 1;
        for &reverse_predecessor in graph.outgoing(target ^ 1) {
            visits.direct_overlap += 1;
            let source = reverse_predecessor ^ 1;
            if source == target ^ 1 || !links[tail_slot(source)].is_empty() {
                continue;
            }
            connect(
                &mut links,
                source,
                target,
                k - 1,
                &mut estimated_bases,
                &mut visits.connect,
            );
            linked += 1;
            direct_made += 1;
            for slot in [head_slot(target), tail_slot(source)] {
                receiver_distance[slot] = UNREACHABLE;
                remaining_receivers -= 1;
            }
            field.remove_source(source, graph, ranges, k, &mut visits.remove_source);
            field.remove_source(target ^ 1, graph, ranges, k, &mut visits.remove_source);
            break;
        }
    }
    // With no free direct edge left, all remaining overlaps cost at least one
    // character. Schedule those receivers for the next distance.
    for target in 0..nodes {
        if receiver_distance[target] == UNREACHABLE {
            continue;
        }
        if k > 1 {
            receiver_distance[target] = 1;
        } else {
            receiver_distance[target] = UNREACHABLE;
            remaining_receivers -= 1;
        }
    }
    info!(
        "After distance {}: {} receiving ends examined, {} links made, {} receiving ends remain, {} free ends, {} estimated bases",
        compact(0u64),
        compact(direct_examined),
        compact(direct_made),
        compact(remaining_receivers),
        compact(nodes - 2 * linked),
        compact(estimated_bases),
    );
    info!(
        "After distance {} node visits: direct_overlap {} (total {}), nearest_receiver {} (total {}), connect {} (total {}), remove_source {} (total {})",
        compact(0u64),
        compact(visits.direct_overlap),
        compact(visits.direct_overlap),
        compact(0u64),
        compact(0u64),
        compact(visits.connect),
        compact(visits.connect),
        compact(visits.remove_source),
        compact(visits.remove_source),
    );

    for distance in 1..k {
        let before = visits;
        let mut examined = 0usize;
        let mut made = 0usize;
        for target in 0..nodes as u32 {
            let index = target as usize;
            if receiver_distance[index] as usize != distance || !links[head_slot(target)].is_empty()
            {
                continue;
            }
            examined += 1;
            let Some((source, actual_distance)) = nearest_receiver(
                target,
                &field,
                &mut reverse_search,
                graph,
                ranges,
                k,
                &mut visits.nearest_receiver,
            ) else {
                receiver_distance[index] = UNREACHABLE;
                remaining_receivers -= 1;
                continue;
            };
            assert!(actual_distance >= distance, "receiver distance decreased");
            if actual_distance > distance {
                receiver_distance[index] = actual_distance as u8;
                continue;
            }
            debug_assert!(links[tail_slot(source)].is_empty());
            let overlap = k - 1 - distance;
            connect(
                &mut links,
                source,
                target,
                overlap,
                &mut estimated_bases,
                &mut visits.connect,
            );
            linked += 1;
            made += 1;
            // A link occupies both physical ends. An end may already have
            // been removed from the receiver schedule as unreachable.
            for slot in [index, tail_slot(source)] {
                if receiver_distance[slot] != UNREACHABLE {
                    receiver_distance[slot] = UNREACHABLE;
                    remaining_receivers -= 1;
                }
            }
            // Both physical ends are independent sources. Repairing them now
            // exposes alternatives to receivers later in this same scan.
            field.remove_source(source, graph, ranges, k, &mut visits.remove_source);
            field.remove_source(target ^ 1, graph, ranges, k, &mut visits.remove_source);
        }
        info!(
            "After distance {}: {} receiving ends examined, {} links made, {} receiving ends remain, {} free ends, {} estimated bases",
            compact(distance),
            compact(examined),
            compact(made),
            compact(remaining_receivers),
            compact(nodes - 2 * linked),
            compact(estimated_bases),
        );
        info!(
            "After distance {} node visits: direct_overlap {} (total {}), nearest_receiver {} (total {}), connect {} (total {}), remove_source {} (total {})",
            compact(distance),
            compact(0u64),
            compact(visits.direct_overlap),
            compact(visits.nearest_receiver - before.nearest_receiver),
            compact(visits.nearest_receiver),
            compact(visits.connect - before.connect),
            compact(visits.connect),
            compact(visits.remove_source - before.remove_source),
            compact(visits.remove_source),
        );
    }

    let outgoing: Vec<_> = (0..nodes as u32)
        .filter(|&id| links[tail_slot(id)].is_empty())
        .collect();
    assert_eq!(outgoing.len() % 2, 0, "odd number of unmatched ends");
    for pair in outgoing.chunks_exact(2) {
        let target = pair[1] ^ 1;
        connect(
            &mut links,
            pair[0],
            target,
            0,
            &mut estimated_bases,
            &mut visits.connect,
        );
    }
    info!(
        "Final node visits (including closure links): direct_overlap {}, nearest_receiver {}, connect {}, remove_source {}",
        compact(visits.direct_overlap),
        compact(visits.nearest_receiver),
        compact(visits.connect),
        compact(visits.remove_source),
    );
    (links, estimated_bases)
}

fn append(
    output: &mut Vec<u8>,
    seq: &PackedSeqVec,
    range: &Range<usize>,
    reverse: bool,
    overlap: usize,
    k: usize,
) {
    let mut bases = seq.slice(range.clone()).unpack();
    if reverse {
        bases.reverse();
        for base in &mut bases {
            *base = complement_char(*base);
        }
    }
    if !output.is_empty() {
        let mask = (k - 1 - overlap).min(bases.len() - overlap);
        for base in &mut bases[overlap..overlap + mask] {
            *base = base.to_ascii_lowercase();
        }
    }
    output.extend_from_slice(&bases[overlap..]);
}

/// Build a masked superstring from unitigs. Newly created crossing k-mers are
/// lowercase, as in `mss::masked_superstring`.
pub fn masked_superstring<K: MssKey>(
    k: usize,
    seq: PackedSeqVec,
    ranges: Vec<Range<usize>>,
) -> Vec<u8>
where
    End<K>: Radixable<K, Key = K>,
{
    assert!(k > 0 && k <= K::BITS / 2);
    let initial_bases: usize = ranges.iter().map(Range::len).sum();
    let (ranges, short): (Vec<_>, Vec<_>) = ranges.into_iter().partition(|r| r.len() >= k);
    assert!(ranges.len() <= (u32::MAX as usize / 2), "too many unitigs");
    let graph = graph::<K>(k, &seq, &ranges);

    let (links, estimated_bases) = match_ends(k, &ranges, &graph, initial_bases);

    info!("Reconstruct output");
    let mut output = Vec::new();
    let mut done = vec![false; ranges.len()];
    for start in 0..ranges.len() {
        if done[start] {
            continue;
        }
        let mut id = (start as u32) * 2;
        let mut overlap = 0;
        loop {
            let index = (id / 2) as usize;
            assert!(!done[index], "link cycle revisits unitig before closing");
            done[index] = true;
            append(&mut output, &seq, &ranges[index], id & 1 != 0, overlap, k);
            let link = links[tail_slot(id)];
            assert!(!link.is_empty(), "unmatched unitig end");
            if (link.target / 2) as usize == start {
                break;
            }
            id = link.target;
            overlap = link.overlap as usize;
        }
    }
    for range in &short {
        if !range.is_empty() {
            append(&mut output, &seq, range, false, 0, k);
        }
    }
    info!(
        "Length summary: {} estimated bases after links, {} output bases",
        compact(estimated_bases),
        compact(output.len()),
    );
    output
}

/// Read unitigs and write one masked superstring in FASTA format.
pub fn run(input: &Path, output: Option<&Path>, k: usize) -> PathBuf {
    let timing = StageTiming::start();
    info!("Reading unitigs..");
    let (seq, ranges) = PackedSeqVec::from_fastx(input);
    let input_bases = ranges.iter().map(Range::len).sum();
    log_file_stats("Read", input, Some((ranges.len(), input_bases))).unwrap();
    let superstring = if k <= 32 {
        masked_superstring::<u64>(k, seq, ranges)
    } else {
        masked_superstring::<u128>(k, seq, ranges)
    };
    let output = output
        .map(Path::to_path_buf)
        .unwrap_or_else(|| input.with_extension("matchtigs.msfa"));
    let mut writer = BufWriter::new(std::fs::File::create(&output).unwrap());
    writeln!(writer, ">matchtigs-masked-superstring").unwrap();
    writer.write_all(&superstring).unwrap();
    writer.write_all(b"\n").unwrap();
    drop(writer);
    log_file_stats("Wrote", &output, Some((1, superstring.len()))).unwrap();
    info!(
        "matchtigs: {}, output {} bases ({})",
        timing.finish(),
        compact(superstring.len()),
        output.display()
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matching_repairs_before_the_next_receiver_at_same_distance() {
        // Each edge has its reverse-complement counterpart. The first match
        // consumes the physical ends represented by sources 0 and 3. Their
        // removal reveals a direct match between source 6 and receiver 4.
        let graph = Graph {
            edges: vec![2, 4, 1, 7, 1, 7, 2, 4],
            offsets: vec![0, 2, 2, 2, 4, 4, 6, 8, 8],
        };
        let k = 5;
        let ranges = vec![0..k; 4];
        let (links, total_bases) = match_ends(k, &ranges, &graph, 4 * k);
        let direct = links
            .iter()
            .filter(|link| link.overlap as usize == k - 1)
            .count();
        assert_eq!(direct, 4);
        assert_eq!(total_bases, 4 * k - 2 * (k - 1));
        let first_target = links[tail_slot(0)].target;
        let second_target = links[tail_slot(6)].target;
        assert_eq!(first_target, 2);
        assert_eq!(second_target, 4);
    }
}
