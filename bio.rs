// Rust implementations of the functions bio specifies: types, coordinate
// mappings, feature trees, gene models, views, and the text formats.
//
// Sourced Rust is included into one pool crate with every other module's, so
// every item here is prefixed `bio_` / `Bio`. Field and tag order are the
// wire contract: a record is positional and a constructor's tag is its
// position, so they follow bio/main.loc exactly.

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct BioGenome;
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct BioTranscript;

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BioStrand { Plus = 0, Minus = 1, Unstranded = 2, Unknown = 3 }

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BioOrient { Fwd = 0, Rev = 1 }

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BioPhase { Phase0 = 0, Phase1 = 1, Phase2 = 2, PhaseNone = 3 }

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct BioSegment {
    pub start: u64,
    pub end: u64,
    pub strand: BioStrand,
}

// A kind of sequence (genome, transcript). Locations and mappings are
// generic over their kinds, as in morloc, so the kinds stay distinct types;
// a sequence id of any kind is an Int.
pub trait BioKind {
    type Id: Copy + PartialEq + std::fmt::Debug;
}
impl BioKind for BioGenome {
    type Id = i64;
}
impl BioKind for BioTranscript {
    type Id = i64;
}

#[derive(Clone, Debug)]
pub struct BioLoc<K: BioKind<Id = i64>> {
    pub seqid: <K as BioKind>::Id,
    pub segment: BioSegment,
}

#[allow(non_snake_case)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct BioBlock {
    pub srcStart: u64,
    pub dstStart: u64,
    pub len: u64,
    pub orient: BioOrient,
}

#[derive(Clone, Debug)]
pub struct BioMapping<A: BioKind<Id = i64>, B: BioKind<Id = i64>> {
    pub src: <A as BioKind>::Id,
    pub dst: <B as BioKind>::Id,
    pub blocks: Vec<BioBlock>,
}

#[derive(Clone, Debug)]
pub enum BioPiece<K: BioKind<Id = i64>> {
    Mapped(std::boxed::Box<(BioLoc<K>,)>),
    Unmapped(std::boxed::Box<(u64, u64)>),
}

#[derive(Clone, Debug)]
pub struct BioSequence<H, S> {
    pub header: H,
    pub bases: S,
}

#[derive(Clone, Debug, Default)]
pub struct BioBaseQual {
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
}

pub type BioFastaEntry = BioSequence<String, Vec<u8>>;
pub type BioFastqEntry = BioSequence<String, BioBaseQual>;

#[derive(Clone, Debug)]
pub struct BioFeature<N> {
    pub node: N,
    pub parts: Vec<BioSegment>,
    pub children: Vec<BioFeature<N>>,
}

#[derive(Clone, Debug)]
pub struct BioRoot<S, N> {
    pub seqid: S,
    pub feature: BioFeature<N>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BioGffNode {
    pub kind: String,
    pub origin: String,
    pub id: String,
    pub attrs: String,
    pub scores: Vec<Option<f64>>,
    pub phases: Vec<BioPhase>,
    pub created: bool,
}

pub type BioGffRoot = BioRoot<String, BioGffNode>;
pub type BioSpan = (u64, u64);

#[derive(Clone, Debug)]
pub struct BioBed3 {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

#[derive(Clone, Debug)]
pub struct BioBed6 {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub score: f64,
    pub strand: BioStrand,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub struct BioBed12 {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub score: f64,
    pub strand: BioStrand,
    pub thickStart: u64,
    pub thickEnd: u64,
    pub rgb: (bool, u8, u8, u8),
    pub blocks: Vec<(u64, u64)>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BioCds {
    pub id: String,
    pub range: BioSpan,
    pub phase: BioPhase,
}

#[derive(Clone, Debug, Default)]
pub struct BioGeneTx {
    pub id: String,
    pub kind: String,
    pub origin: String,
    pub attrs: String,
    pub exons: Vec<BioSpan>,
    pub cds: Option<BioCds>,
}

#[derive(Clone, Debug)]
pub struct BioGene {
    pub seqid: String,
    pub id: String,
    pub kind: String,
    pub origin: String,
    pub attrs: String,
    pub span: BioSpan,
    pub strand: BioStrand,
    pub transcripts: Vec<BioGeneTx>,
}

#[allow(non_snake_case)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BioComplexity {
    pub mixedStrand: bool,
    pub splitFeature: bool,
    pub missingExons: bool,
    pub overlappingExons: bool,
    pub spanMismatch: bool,
    pub multipleCds: bool,
    pub cdsNotInExons: bool,
    pub phaseShift: bool,
    pub underivedLines: bool,
    pub extraChildren: bool,
    pub extraColumns: bool,
    pub anonymousTranscript: bool,
}

#[derive(Clone, Debug)]
pub enum BioGffEntry {
    SimpleCodingGene(std::boxed::Box<(BioGene,)>),
    SimpleNonCodingGene(std::boxed::Box<(BioGene,)>),
    ComplexCodingGene(std::boxed::Box<(BioComplexity, BioGffRoot)>),
    ComplexNonCodingGene(std::boxed::Box<(BioComplexity, BioGffRoot)>),
    OtherFeature(std::boxed::Box<(BioGffRoot,)>),
}

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BioLibStrand { LibUnstranded = 0, LibForward = 1, LibReverse = 2 }

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BioLayout { PairedEnd = 0, SingleEnd = 1 }

#[allow(non_snake_case)]
#[derive(Clone, Debug, PartialEq)]
pub struct BioTxHit {
    pub tx: i64,
    pub pos: i64,
    pub fwd: bool,
    pub fragLen: i64,
    pub mates: i64,
    pub score: i64,
    pub weight: f64,
    pub format: i64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BioFragHits {
    pub hits: Vec<BioTxHit>,
    pub decoy: bool,
}

#[derive(Clone, Debug)]
pub struct BioTxModel {
    pub id: String,
    pub gene: String,
    pub seqid: String,
    pub mapping: BioMapping<BioTranscript, BioGenome>,
    pub cds: Option<(u64, u64)>,
}


// Called by modules whose Rust names bio-rust types, so that this file is
// included in their pools (compiler model/runtime/issues/build-and-text.md).
pub fn bio_rust_anchor(x: i64) -> i64 {
    x
}


// Coordinate mappings
// -------------------

fn bio_flip(s: BioStrand) -> BioStrand {
    match s {
        BioStrand::Plus => BioStrand::Minus,
        BioStrand::Minus => BioStrand::Plus,
        x => x,
    }
}

fn bio_orient_mul(a: BioOrient, b: BioOrient) -> BioOrient {
    if a == b { BioOrient::Fwd } else { BioOrient::Rev }
}

fn bio_sort_blocks(bs: &mut Vec<BioBlock>) {
    bs.sort_by_key(|b| b.srcStart);
}

// The image in B of A's [s, e), which must lie within the block.
fn bio_block_image(b: &BioBlock, s: u64, e: u64) -> (u64, u64) {
    if b.orient == BioOrient::Fwd {
        (b.dstStart + (s - b.srcStart), b.dstStart + (e - b.srcStart))
    } else {
        let src_end = b.srcStart + b.len;
        (b.dstStart + (src_end - e), b.dstStart + (src_end - s))
    }
}

pub fn bio_tiled<A: BioKind<Id = i64>, B: BioKind<Id = i64>>(src: i64, dst: i64, parts: &Vec<(u64, u64, BioOrient)>) -> BioMapping<A, B> {
    let mut blocks = Vec::with_capacity(parts.len());
    let mut at = 0u64;
    for &(dst_start, len, orient) in parts {
        blocks.push(BioBlock { srcStart: at, dstStart: dst_start, len, orient });
        at += len;
    }
    BioMapping { src, dst, blocks }
}

pub fn bio_invert<A: BioKind<Id = i64>, B: BioKind<Id = i64>>(m: &BioMapping<A, B>) -> BioMapping<B, A> {
    let mut blocks: Vec<BioBlock> = m
        .blocks
        .iter()
        .map(|b| BioBlock { srcStart: b.dstStart, dstStart: b.srcStart, len: b.len, orient: b.orient })
        .collect();
    bio_sort_blocks(&mut blocks);
    BioMapping { src: m.dst, dst: m.src, blocks }
}

pub fn bio_compose_blocks<A: BioKind<Id = i64>, B: BioKind<Id = i64>, C: BioKind<Id = i64>>(ab: &BioMapping<A, B>, bc: &BioMapping<B, C>) -> BioMapping<A, C> {
    let mut blocks = Vec::new();
    for x in &ab.blocks {
        let (x0, x1) = (x.dstStart, x.dstStart + x.len);
        let x_inv = BioBlock { srcStart: x.dstStart, dstStart: x.srcStart, len: x.len, orient: x.orient };
        for y in &bc.blocks {
            let is = x0.max(y.srcStart);
            let ie = x1.min(y.srcStart + y.len);
            if is >= ie {
                continue;
            }
            let a = bio_block_image(&x_inv, is, ie);
            let c = bio_block_image(y, is, ie);
            blocks.push(BioBlock { srcStart: a.0, dstStart: c.0, len: ie - is, orient: bio_orient_mul(x.orient, y.orient) });
        }
    }
    bio_sort_blocks(&mut blocks);
    BioMapping { src: ab.src, dst: bc.dst, blocks }
}

fn bio_project_from<A: BioKind<Id = i64>, B: BioKind<Id = i64>>(m: &BioMapping<A, B>, loc: &BioLoc<A>, first: usize) -> Vec<BioPiece<B>> {
    if loc.seqid != m.src {
        panic!("project: the location is not on the mapping's source sequence");
    }
    let (s, e) = (loc.segment.start, loc.segment.end);
    if s > e {
        panic!("project: segment start is past its end");
    }
    let mut out = Vec::new();
    let mapped = |out: &mut Vec<BioPiece<B>>, img: (u64, u64), o: BioOrient| {
        let st = if o == BioOrient::Rev { bio_flip(loc.segment.strand) } else { loc.segment.strand };
        out.push(BioPiece::Mapped(Box::new((BioLoc {
            seqid: m.dst,
            segment: BioSegment { start: img.0, end: img.1, strand: st },
        },))));
    };
    if s == e {
        let mut i = first;
        while i < m.blocks.len() && m.blocks[i].srcStart <= s {
            let b = &m.blocks[i];
            if s <= b.srcStart + b.len {
                mapped(&mut out, bio_block_image(b, s, s), b.orient);
                return out;
            }
            i += 1;
        }
        out.push(BioPiece::Unmapped(Box::new((s, s))));
        return out;
    }
    let mut covered = s;
    let mut i = first;
    while i < m.blocks.len() && m.blocks[i].srcStart < e {
        let b = &m.blocks[i];
        i += 1;
        let bs = s.max(b.srcStart);
        let be = e.min(b.srcStart + b.len);
        if bs >= be {
            continue;
        }
        if bs > covered {
            out.push(BioPiece::Unmapped(Box::new((covered, bs))));
        }
        mapped(&mut out, bio_block_image(b, bs, be), b.orient);
        covered = covered.max(be);
    }
    if covered < e {
        out.push(BioPiece::Unmapped(Box::new((covered, e))));
    }
    out
}

pub fn bio_project<A: BioKind<Id = i64>, B: BioKind<Id = i64>>(m: &BioMapping<A, B>, loc: &BioLoc<A>) -> Vec<BioPiece<B>> {
    bio_project_from(m, loc, 0)
}

pub fn bio_project_many<A: BioKind<Id = i64>, B: BioKind<Id = i64>>(m: &BioMapping<A, B>, locs: &Vec<BioLoc<A>>) -> Vec<Vec<BioPiece<B>>> {
    let mut reach = Vec::with_capacity(m.blocks.len());
    let mut hi = 0u64;
    for (i, b) in m.blocks.iter().enumerate() {
        if i > 0 && b.srcStart < m.blocks[i - 1].srcStart {
            panic!("projectMany: the mapping's blocks are not sorted by srcStart");
        }
        hi = hi.max(b.srcStart + b.len);
        reach.push(hi);
    }
    locs.iter()
        .map(|loc| {
            let s = loc.segment.start;
            let point = s == loc.segment.end;
            let first = reach.partition_point(|&r| if point { r < s } else { r <= s });
            bio_project_from(m, loc, first)
        })
        .collect()
}

// The IUPAC complement of a base, case kept. Other bytes are unchanged.
#[inline]
pub fn bio_complement(c: u8) -> u8 {
    match c {
        b'A' => b'T', b'T' => b'A', b'U' => b'A', b'G' => b'C', b'C' => b'G',
        b'R' => b'Y', b'Y' => b'R', b'K' => b'M', b'M' => b'K',
        b'B' => b'V', b'V' => b'B', b'D' => b'H', b'H' => b'D',
        b'a' => b't', b't' => b'a', b'u' => b'a', b'g' => b'c', b'c' => b'g',
        b'r' => b'y', b'y' => b'r', b'k' => b'm', b'm' => b'k',
        b'b' => b'v', b'v' => b'b', b'd' => b'h', b'h' => b'd',
        x => x,
    }
}

pub fn bio_extract<A: BioKind<Id = i64>, B: BioKind<Id = i64>>(m: &BioMapping<A, B>, dst: &Vec<u8>) -> Vec<u8> {
    let mut out = Vec::new();
    for b in &m.blocks {
        let (s, e) = (b.dstStart as usize, (b.dstStart + b.len) as usize);
        if e > dst.len() {
            panic!("extract: a block runs past the end of the sequence");
        }
        if b.orient == BioOrient::Fwd {
            out.extend_from_slice(&dst[s..e]);
        } else {
            out.extend(dst[s..e].iter().rev().map(|&c| bio_complement(c)));
        }
    }
    out
}

pub fn bio_extract_blocks<A: BioKind<Id = i64>, B: BioKind<Id = i64>>(m: &BioMapping<A, B>, slices: &Vec<Vec<u8>>) -> Vec<u8> {
    if slices.len() != m.blocks.len() {
        panic!("extractBlocks: expected one slice per block");
    }
    let mut out = Vec::new();
    for (b, sl) in m.blocks.iter().zip(slices) {
        if sl.len() as u64 != b.len {
            panic!("extractBlocks: a slice's length differs from its block's");
        }
        if b.orient == BioOrient::Fwd {
            out.extend_from_slice(sl);
        } else {
            out.extend(sl.iter().rev().map(|&c| bio_complement(c)));
        }
    }
    out
}


// Features
// --------

// Parts in transcript order: ascending, or descending when every part is on
// the minus strand. Returns the new order of the old indices.
fn bio_transcript_order(parts: &[BioSegment]) -> Vec<usize> {
    let mut perm: Vec<usize> = (0..parts.len()).collect();
    perm.sort_by_key(|&i| parts[i].start);
    if !parts.is_empty() && parts.iter().all(|p| p.strand == BioStrand::Minus) {
        perm.reverse();
    }
    perm
}

fn bio_permute<T: Clone>(xs: &[T], perm: &[usize]) -> Vec<T> {
    perm.iter().map(|&i| xs[i].clone()).collect()
}

fn bio_spans<N>(f: &BioFeature<N>) -> Vec<BioSpan> {
    f.parts.iter().map(|p| (p.start, p.end)).collect()
}

fn bio_feature<N>(node: N, segs: &[BioSegment]) -> BioFeature<N> {
    BioFeature { node, parts: bio_permute(segs, &bio_transcript_order(segs)), children: Vec::new() }
}


// Views
// -----

fn bio_child_parts(f: &BioFeature<BioGffNode>, kind: &str) -> Vec<BioSegment> {
    let mut parts = Vec::new();
    for c in &f.children {
        if c.node.kind == kind {
            parts.extend_from_slice(&c.parts);
        }
    }
    parts
}

fn bio_tx_model(
    id: String,
    gene: String,
    seqid: String,
    src: i64,
    dst: i64,
    mut exons: Vec<BioSegment>,
    cds_parts: &[BioSegment],
) -> BioTxModel {
    exons.sort_by_key(|e| e.start);
    let st = exons[0].strand;
    for e in &exons {
        if e.strand != st || (st != BioStrand::Plus && st != BioStrand::Minus) {
            panic!("transcriptModels: {} has exons on mixed or unknown strands", id);
        }
    }
    if st == BioStrand::Minus {
        exons.reverse();
    }
    let o = if st == BioStrand::Minus { BioOrient::Rev } else { BioOrient::Fwd };
    let parts: Vec<(u64, u64, BioOrient)> = exons.iter().map(|e| (e.start, e.end - e.start, o)).collect();
    let mapping: BioMapping<BioTranscript, BioGenome> = bio_tiled(src, dst, &parts);
    let mut cds = None;
    if !cds_parts.is_empty() {
        let inv = bio_invert(&mapping);
        let (mut lo, mut hi) = (u64::MAX, 0u64);
        for p in cds_parts {
            let loc: BioLoc<BioGenome> = BioLoc { seqid: dst, segment: *p };
            for piece in bio_project(&inv, &loc) {
                match piece {
                    BioPiece::Mapped(b) => {
                        lo = lo.min(b.0.segment.start);
                        hi = hi.max(b.0.segment.end);
                    }
                    BioPiece::Unmapped(_) => panic!("transcriptModels: a CDS of {} lies outside its exons", id),
                }
            }
        }
        cds = Some((lo, hi));
    }
    BioTxModel { id, gene, seqid, mapping, cds }
}

fn bio_seq_index(index: &std::collections::HashMap<&str, i64>, seqid: &str, what: &str) -> i64 {
    match index.get(seqid) {
        Some(&i) => i,
        None => panic!("transcriptModels: {} is on {}, which is not a genome sequence", what, seqid),
    }
}

fn bio_root_tx_models(index: &std::collections::HashMap<&str, i64>, root: &BioGffRoot, out: &mut Vec<BioTxModel>) {
    let g = &root.feature;
    let gene = &g.node.id;
    let mut txs: Vec<&BioFeature<BioGffNode>> = Vec::new();
    if !bio_child_parts(g, "exon").is_empty() {
        txs.push(g);
    }
    for c in &g.children {
        if !bio_child_parts(c, "exon").is_empty() {
            txs.push(c);
        }
    }
    for t in txs {
        let id = if !t.node.id.is_empty() { t.node.id.clone() } else { format!("{}:{}", gene, out.len()) };
        let dst = bio_seq_index(index, &root.seqid, &id);
        let src = out.len() as i64;
        out.push(bio_tx_model(
            id,
            gene.clone(),
            root.seqid.clone(),
            src,
            dst,
            bio_child_parts(t, "exon"),
            &bio_child_parts(t, "CDS"),
        ));
    }
}

// A standard transcript's CDS parts: its exons clipped to the coding range.
fn bio_cds_spans(exons: &[BioSpan], cds: &Option<BioCds>) -> Vec<BioSpan> {
    let mut out = Vec::new();
    if let Some(c) = cds {
        let (cs, ce) = c.range;
        for &(s, e) in exons {
            let (a, b) = (s.max(cs), e.min(ce));
            if a < b {
                out.push((a, b));
            }
        }
    }
    out
}

pub fn bio_entry_seqid(e: &BioGffEntry) -> String {
    match e {
        BioGffEntry::SimpleCodingGene(g) | BioGffEntry::SimpleNonCodingGene(g) => g.0.seqid.clone(),
        BioGffEntry::ComplexCodingGene(r) | BioGffEntry::ComplexNonCodingGene(r) => r.1.seqid.clone(),
        BioGffEntry::OtherFeature(r) => r.0.seqid.clone(),
    }
}

pub fn bio_transcript_models(names: &Vec<String>, entries: &Vec<BioGffEntry>) -> Vec<BioTxModel> {
    let mut index: std::collections::HashMap<&str, i64> = std::collections::HashMap::new();
    for (i, n) in names.iter().enumerate() {
        index.insert(n.as_str(), i as i64);
    }
    let mut out: Vec<BioTxModel> = Vec::new();
    for e in entries {
        match e {
            BioGffEntry::SimpleCodingGene(b) | BioGffEntry::SimpleNonCodingGene(b) => {
                let g = &b.0;
                for t in &g.transcripts {
                    let id = if !t.id.is_empty() {
                        t.id.clone()
                    } else if !g.id.is_empty() {
                        g.id.clone()
                    } else {
                        format!("{}:{}", g.id, out.len())
                    };
                    let dst = bio_seq_index(&index, &g.seqid, &id);
                    let exons: Vec<BioSegment> =
                        t.exons.iter().map(|&(s, e)| BioSegment { start: s, end: e, strand: g.strand }).collect();
                    let cds: Vec<BioSegment> = bio_cds_spans(&t.exons, &t.cds)
                        .iter()
                        .map(|&(s, e)| BioSegment { start: s, end: e, strand: g.strand })
                        .collect();
                    let src = out.len() as i64;
                    out.push(bio_tx_model(id, g.id.clone(), g.seqid.clone(), src, dst, exons, &cds));
                }
            }
            BioGffEntry::ComplexCodingGene(r) | BioGffEntry::ComplexNonCodingGene(r) => {
                bio_root_tx_models(&index, &r.1, &mut out)
            }
            BioGffEntry::OtherFeature(r) => bio_root_tx_models(&index, &r.0, &mut out),
        }
    }
    out
}


// Gene models
// -----------

fn bio_is_codon(k: &str) -> bool {
    k == "start_codon" || k == "stop_codon"
}

fn bio_sorted_spans(mut xs: Vec<BioSpan>) -> Vec<BioSpan> {
    xs.sort();
    xs
}

// The spans of `n` bases of `spans` (given in transcript order) starting
// `skip` bases in, in transcript order; fewer when they run out.
fn bio_take_bases(spans: &[BioSpan], minus: bool, mut skip: u64, mut n: u64) -> Vec<BioSpan> {
    let mut out = Vec::new();
    for &(s, e) in spans {
        if n == 0 {
            break;
        }
        let len = e - s;
        if skip >= len {
            skip -= len;
            continue;
        }
        let take = n.min(len - skip);
        if minus {
            out.push((e - skip - take, e - skip));
        } else {
            out.push((s + skip, s + skip + take));
        }
        n -= take;
        skip = 0;
    }
    out
}

// The UTR spans of a coding transcript on either side of its CDS, each in
// transcript order: (five prime, three prime).
fn bio_utr_spans(exons_asc: &[BioSpan], cds: BioSpan, minus: bool) -> (Vec<BioSpan>, Vec<BioSpan>) {
    let (cs, ce) = cds;
    let mut lo = Vec::new();
    let mut hi = Vec::new();
    for &(s, e) in exons_asc {
        if s < cs {
            lo.push((s, e.min(cs)));
        }
        if e > ce {
            hi.push((s.max(ce), e));
        }
    }
    if minus {
        lo.reverse();
        hi.reverse();
        (hi, lo)
    } else {
        (lo, hi)
    }
}

fn bio_phase_of(n: u64) -> BioPhase {
    match n {
        0 => BioPhase::Phase0,
        1 => BioPhase::Phase1,
        _ => BioPhase::Phase2,
    }
}

// Phases of CDS parts in transcript order, derived from the first.
fn bio_derived_phases(parts_tx: &[BioSpan], first: BioPhase) -> Vec<BioPhase> {
    let mut out = Vec::new();
    if first == BioPhase::PhaseNone {
        return out;
    }
    let mut done = 0u64;
    let p0 = first as u64;
    for (i, &(s, e)) in parts_tx.iter().enumerate() {
        if i == 0 {
            out.push(first);
        } else {
            let left = (done + 3 - p0 % 3) % 3;
            out.push(bio_phase_of((3 - left) % 3));
        }
        done += e - s;
    }
    out
}

fn bio_any_flag(c: &BioComplexity) -> bool {
    c.mixedStrand || c.splitFeature || c.missingExons || c.overlappingExons || c.spanMismatch
        || c.multipleCds || c.cdsNotInExons || c.phaseShift || c.underivedLines || c.extraChildren
        || c.extraColumns || c.anonymousTranscript
}

fn bio_tree_has(f: &BioFeature<BioGffNode>, kind: &str) -> bool {
    f.node.kind == kind || f.children.iter().any(|c| bio_tree_has(c, kind))
}

fn bio_tree_segments(f: &BioFeature<BioGffNode>, out: &mut Vec<BioSegment>) {
    out.extend_from_slice(&f.parts);
    for c in &f.children {
        bio_tree_segments(c, out);
    }
}

// The feature of a built node: parts in transcript order, with scores and
// phases in the same order, each dropped when every line's is '.'.
fn bio_built_feature(node: &BioGffNode, parts: &[BioSegment]) -> BioFeature<BioGffNode> {
    let perm = bio_transcript_order(parts);
    let mut n = node.clone();
    n.scores = bio_permute(&node.scores, &perm);
    n.phases = bio_permute(&node.phases, &perm);
    if n.scores.iter().all(|x| x.is_none()) {
        n.scores.clear();
    }
    if n.phases.iter().all(|&x| x == BioPhase::PhaseNone) {
        n.phases.clear();
    }
    BioFeature { node: n, parts: bio_permute(parts, &perm), children: Vec::new() }
}

// Checks one transcript -- `t`'s children, with `t` itself absent for a gene
// whose CDS lines sit directly under it -- and builds its model.
fn bio_check_tx(
    t: Option<&BioFeature<BioGffNode>>,
    kids: &[&BioFeature<BioGffNode>],
    gene: &BioGffNode,
    minus: bool,
    c: &mut BioComplexity,
) -> BioGeneTx {
    let mut tx = BioGeneTx::default();
    match t {
        Some(t) => {
            tx.id = t.node.id.clone();
            tx.kind = t.node.kind.clone();
            tx.origin = t.node.origin.clone();
            tx.attrs = t.node.attrs.clone();
            if tx.id.is_empty() {
                c.anonymousTranscript = true;
            }
            if t.parts.len() != 1 {
                c.splitFeature = true;
            }
        }
        None => {
            tx.origin = if kids.is_empty() { gene.origin.clone() } else { kids[0].node.origin.clone() };
        }
    }
    let mut cds_nodes: Vec<&BioFeature<BioGffNode>> = Vec::new();
    let (mut five, mut three, mut generic, mut start, mut stop) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for &k in kids {
        let n = &k.node;
        if !k.children.is_empty() {
            c.extraChildren = true;
        }
        if !n.scores.is_empty() || n.origin != tx.origin {
            c.extraColumns = true;
        }
        if !n.phases.is_empty() && n.kind != "CDS" {
            if !bio_is_codon(&n.kind) || bio_derived_phases(&bio_spans(k), BioPhase::Phase0) != n.phases {
                c.extraColumns = true;
            }
        }
        match n.kind.as_str() {
            "exon" => {
                if k.parts.len() != 1 {
                    c.splitFeature = true;
                }
                tx.exons.extend(bio_spans(k));
            }
            "CDS" => cds_nodes.push(k),
            "five_prime_UTR" => five.extend(bio_spans(k)),
            "three_prime_UTR" => three.extend(bio_spans(k)),
            "UTR" => generic.extend(bio_spans(k)),
            "start_codon" => start.extend(bio_spans(k)),
            "stop_codon" => stop.extend(bio_spans(k)),
            _ => c.extraChildren = true,
        }
    }
    if t.is_none() {
        for k in &cds_nodes {
            tx.exons.extend(bio_spans(k));
        }
    }
    tx.exons.sort();
    if tx.exons.is_empty() {
        c.missingExons = true;
    }
    for i in 1..tx.exons.len() {
        if tx.exons[i].0 < tx.exons[i - 1].1 {
            c.overlappingExons = true;
        }
    }
    if let Some(t) = t {
        if t.parts.len() == 1 && !tx.exons.is_empty()
            && (t.parts[0].start, t.parts[0].end) != (tx.exons[0].0, tx.exons[tx.exons.len() - 1].1)
        {
            c.spanMismatch = true;
        }
    }
    // CDS lines without an ID are one CDS written line by line.
    let merged;
    if cds_nodes.len() > 1 && cds_nodes.iter().all(|k| k.node.id.is_empty()) {
        let mut node = cds_nodes[0].node.clone();
        node.scores.clear();
        node.phases.clear();
        let mut parts = Vec::new();
        for k in &cds_nodes {
            for i in 0..k.parts.len() {
                parts.push(k.parts[i]);
                node.scores.push(if k.node.scores.is_empty() { None } else { k.node.scores[i] });
                node.phases.push(if k.node.phases.is_empty() { BioPhase::PhaseNone } else { k.node.phases[i] });
            }
        }
        merged = bio_built_feature(&node, &parts);
        cds_nodes = vec![&merged];
    }
    if cds_nodes.len() > 1 {
        c.multipleCds = true;
    }
    let mut exons_tx = tx.exons.clone();
    if minus {
        exons_tx.reverse();
    }
    if cds_nodes.len() == 1 && !tx.exons.is_empty() {
        let cn = cds_nodes[0];
        let cn_spans = bio_spans(cn);
        let cds_asc = bio_sorted_spans(cn_spans.clone());
        let range = (cds_asc[0].0, cds_asc[cds_asc.len() - 1].1);
        let probe = Some(BioCds { id: cn.node.id.clone(), range, phase: BioPhase::PhaseNone });
        if bio_cds_spans(&tx.exons, &probe) != cds_asc {
            c.cdsNotInExons = true;
        }
        let first = if cn.node.phases.is_empty() { BioPhase::PhaseNone } else { cn.node.phases[0] };
        if !cn.node.phases.is_empty() && bio_derived_phases(&cn_spans, first) != cn.node.phases {
            c.phaseShift = true;
        }
        tx.cds = Some(BioCds { id: cn.node.id.clone(), range, phase: first });
        let (f5, f3) = bio_utr_spans(&tx.exons, range, minus);
        let mut both = f5.clone();
        both.extend_from_slice(&f3);
        if (!five.is_empty() || !three.is_empty())
            && (bio_sorted_spans(five) != bio_sorted_spans(f5) || bio_sorted_spans(three) != bio_sorted_spans(f3))
        {
            c.underivedLines = true;
        }
        if !generic.is_empty() && bio_sorted_spans(generic) != bio_sorted_spans(both) {
            c.underivedLines = true;
        }
        let cds_tx = &cn_spans;
        if !start.is_empty() && bio_sorted_spans(start) != bio_sorted_spans(bio_take_bases(cds_tx, minus, 0, 3)) {
            c.underivedLines = true;
        }
        if !stop.is_empty() {
            let cds_len: u64 = cds_tx.iter().map(|&(s, e)| e - s).sum();
            let st = bio_sorted_spans(stop);
            let last3 = cds_len >= 3 && st == bio_sorted_spans(bio_take_bases(cds_tx, minus, cds_len - 3, 3));
            let mut before = 0u64;
            for &(s, e) in &exons_tx {
                if minus {
                    before += if e <= range.0 { 0 } else { e - s.max(range.0) };
                } else {
                    before += if s >= range.1 { 0 } else { e.min(range.1) - s };
                }
            }
            let after3 = st == bio_sorted_spans(bio_take_bases(&exons_tx, minus, before, 3));
            if !last3 && !after3 {
                c.underivedLines = true;
            }
        }
    } else if !five.is_empty() || !three.is_empty() || !generic.is_empty() || !start.is_empty() || !stop.is_empty() {
        c.underivedLines = true;
    }
    tx
}

fn bio_gff_entry(r: &BioGffRoot) -> BioGffEntry {
    let g = &r.feature;
    let other = || BioGffEntry::OtherFeature(Box::new((r.clone(),)));
    if g.children.is_empty() {
        return other();
    }
    let mut txs: Vec<&BioFeature<BioGffNode>> = Vec::new();
    let mut leaves: Vec<&BioFeature<BioGffNode>> = Vec::new();
    for ch in &g.children {
        if ch.children.is_empty() { leaves.push(ch) } else { txs.push(ch) }
    }
    let implicit = txs.is_empty();
    if implicit && !leaves.iter().all(|k| k.node.kind == "CDS") {
        return other();
    }
    let mut c = BioComplexity::default();
    let mut all = Vec::new();
    bio_tree_segments(g, &mut all);
    let st = all[0].strand;
    if all.iter().any(|p| p.strand != st) {
        c.mixedStrand = true;
    }
    let minus = st == BioStrand::Minus;
    if g.parts.len() != 1 {
        c.splitFeature = true;
    }
    if !g.node.scores.is_empty() || !g.node.phases.is_empty() {
        c.extraColumns = true;
    }
    if !implicit && !leaves.is_empty() {
        c.extraChildren = true;
    }
    for t in &txs {
        if !t.node.scores.is_empty() || !t.node.phases.is_empty() {
            c.extraColumns = true;
        }
    }
    let mut gene = BioGene {
        seqid: r.seqid.clone(),
        id: g.node.id.clone(),
        kind: g.node.kind.clone(),
        origin: g.node.origin.clone(),
        attrs: g.node.attrs.clone(),
        span: if g.parts.is_empty() { (0, 0) } else { (g.parts[0].start, g.parts[0].end) },
        strand: st,
        transcripts: Vec::new(),
    };
    if implicit {
        gene.transcripts.push(bio_check_tx(None, &leaves, &g.node, minus, &mut c));
    } else {
        for t in &txs {
            let kids: Vec<&BioFeature<BioGffNode>> = t.children.iter().collect();
            gene.transcripts.push(bio_check_tx(Some(t), &kids, &g.node, minus, &mut c));
        }
    }
    let coding = bio_tree_has(g, "CDS");
    if bio_any_flag(&c) {
        if coding {
            BioGffEntry::ComplexCodingGene(Box::new((c, r.clone())))
        } else {
            BioGffEntry::ComplexNonCodingGene(Box::new((c, r.clone())))
        }
    } else if coding {
        BioGffEntry::SimpleCodingGene(Box::new((gene,)))
    } else {
        BioGffEntry::SimpleNonCodingGene(Box::new((gene,)))
    }
}

pub fn bio_gff_entries(roots: &Vec<BioGffRoot>) -> Vec<BioGffEntry> {
    roots.iter().map(bio_gff_entry).collect()
}

fn bio_model_leaf(kind: &str, id: &str, origin: &str, st: BioStrand, parts_tx: &[BioSpan], phases: Vec<BioPhase>) -> BioFeature<BioGffNode> {
    BioFeature {
        node: BioGffNode {
            kind: kind.to_string(),
            origin: origin.to_string(),
            id: id.to_string(),
            attrs: String::new(),
            scores: Vec::new(),
            phases,
            created: false,
        },
        parts: parts_tx.iter().map(|&(a, b)| BioSegment { start: a, end: b, strand: st }).collect(),
        children: Vec::new(),
    }
}

// The tree of a standard gene: the lines a GFF3 has for it.
fn bio_gene_root(g: &BioGene) -> BioGffRoot {
    let minus = g.strand == BioStrand::Minus;
    let tx_order = |mut xs: Vec<BioSpan>| {
        if minus {
            xs.reverse();
        }
        xs
    };
    let mut gf = BioFeature {
        node: BioGffNode {
            kind: g.kind.clone(),
            origin: g.origin.clone(),
            id: g.id.clone(),
            attrs: g.attrs.clone(),
            scores: Vec::new(),
            phases: Vec::new(),
            created: false,
        },
        parts: vec![BioSegment { start: g.span.0, end: g.span.1, strand: g.strand }],
        children: Vec::new(),
    };
    for t in &g.transcripts {
        let mut kids = Vec::new();
        let exons_tx = tx_order(t.exons.clone());
        if !t.id.is_empty() {
            for e in &exons_tx {
                kids.push(bio_model_leaf("exon", "", &t.origin, g.strand, &[*e], Vec::new()));
            }
        }
        if let Some(cds) = &t.cds {
            let cds_tx = tx_order(bio_cds_spans(&t.exons, &t.cds));
            let ph = bio_derived_phases(&cds_tx, cds.phase);
            kids.push(bio_model_leaf("CDS", &cds.id, &t.origin, g.strand, &cds_tx, ph));
            if !t.id.is_empty() {
                let (f5, f3) = bio_utr_spans(&t.exons, cds.range, minus);
                for u in &f5 {
                    kids.push(bio_model_leaf("five_prime_UTR", "", &t.origin, g.strand, &[*u], Vec::new()));
                }
                for u in &f3 {
                    kids.push(bio_model_leaf("three_prime_UTR", "", &t.origin, g.strand, &[*u], Vec::new()));
                }
            }
        }
        if t.id.is_empty() {
            gf.children.extend(kids);
            continue;
        }
        gf.children.push(BioFeature {
            node: BioGffNode {
                kind: t.kind.clone(),
                origin: t.origin.clone(),
                id: t.id.clone(),
                attrs: t.attrs.clone(),
                scores: Vec::new(),
                phases: Vec::new(),
                created: false,
            },
            parts: vec![BioSegment { start: t.exons[0].0, end: t.exons[t.exons.len() - 1].1, strand: g.strand }],
            children: kids,
        });
    }
    BioRoot { seqid: g.seqid.clone(), feature: gf }
}

pub fn bio_entry_roots(entries: &Vec<BioGffEntry>) -> Vec<BioGffRoot> {
    entries
        .iter()
        .map(|e| match e {
            BioGffEntry::SimpleCodingGene(g) | BioGffEntry::SimpleNonCodingGene(g) => bio_gene_root(&g.0),
            BioGffEntry::ComplexCodingGene(r) | BioGffEntry::ComplexNonCodingGene(r) => r.1.clone(),
            BioGffEntry::OtherFeature(r) => r.0.clone(),
        })
        .collect()
}

// Text formats
// ============
//
// Line reading (plain or gzip), atomic writers, FASTA, FASTQ, GFF3, GTF and
// BED, as streaming readers and batch formatters.

// Sequence batches hold about this many bytes; annotation and BED batches
// this many records.
const BIO_BATCH_BYTES: usize = 8 << 20;
const BIO_BATCH_TREES: usize = 4096;

// Reads lines from a plain or gzip-compressed file; concatenated gzip members
// read as one stream.
pub struct BioLineReader {
    path: String,
    inner: Box<dyn std::io::BufRead + Send>,
    lineno: usize,
}

impl BioLineReader {
    pub fn open(path: &str) -> BioLineReader {
        use std::io::BufRead;
        let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("cannot open {}: {}", path, e));
        let mut br = std::io::BufReader::with_capacity(1 << 18, f);
        let gz = {
            let head = br.fill_buf().unwrap_or_else(|e| panic!("read error in {}: {}", path, e));
            head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
        };
        let inner: Box<dyn std::io::BufRead + Send> = if gz {
            Box::new(BioInflated::spawn(br))
        } else {
            Box::new(br)
        };
        BioLineReader { path: path.to_string(), inner, lineno: 0 }
    }

    // The next line, without its line ending, into `line`; false at end of file.
    pub fn next(&mut self, line: &mut Vec<u8>) -> bool {
        use std::io::BufRead;
        line.clear();
        let n = self
            .inner
            .read_until(b'\n', line)
            .unwrap_or_else(|e| panic!("read error in {}: {}", self.path, e));
        if n == 0 {
            return false;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        self.lineno += 1;
        true
    }

    // "path:line: msg", for the line last returned.
    pub fn fail(&self, msg: &str) -> ! {
        panic!("{}:{}: {}", self.path, self.lineno, msg)
    }
}

// The decompressed bytes of a gzip stream (concatenated members read as
// one), inflated in blocks. An empty block marks the end, so a lost
// inflater is an error, not an end.
struct BioInflated {
    blocks: std::sync::mpsc::Receiver<std::io::Result<Vec<u8>>>,
    buf: Vec<u8>,
    pos: usize,
    done: bool,
}

const BIO_INFLATE_BLOCK: usize = 1 << 20;

impl BioInflated {
    fn spawn<R: std::io::BufRead + Send + 'static>(input: R) -> BioInflated {
        use std::io::Read;
        let (tx, rx) = std::sync::mpsc::sync_channel(4);
        std::thread::spawn(move || {
            let mut d = flate2::bufread::MultiGzDecoder::new(input);
            loop {
                let mut b = vec![0u8; BIO_INFLATE_BLOCK];
                let mut n = 0;
                let r = loop {
                    match d.read(&mut b[n..]) {
                        Ok(0) => break Ok(()),
                        Ok(k) => {
                            n += k;
                            if n == b.len() {
                                break Ok(());
                            }
                        }
                        Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                        Err(e) => break Err(e),
                    }
                };
                let end = n < b.len();
                b.truncate(n);
                let sent = match r {
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        return;
                    }
                    Ok(()) if n > 0 => tx.send(Ok(b)),
                    Ok(()) => Ok(()),
                };
                if sent.is_err() {
                    return;
                }
                if end {
                    let _ = tx.send(Ok(Vec::new()));
                    return;
                }
            }
        });
        BioInflated { blocks: rx, buf: Vec::new(), pos: 0, done: false }
    }
}

impl std::io::Read for BioInflated {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        use std::io::BufRead;
        let n = {
            let b = self.fill_buf()?;
            let n = b.len().min(out.len());
            out[..n].copy_from_slice(&b[..n]);
            n
        };
        self.consume(n);
        Ok(n)
    }
}

impl std::io::BufRead for BioInflated {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        while self.pos == self.buf.len() && !self.done {
            match self.blocks.recv() {
                Ok(Ok(b)) if b.is_empty() => self.done = true,
                Ok(Ok(b)) => {
                    self.buf = b;
                    self.pos = 0;
                }
                Ok(Err(e)) => {
                    self.done = true;
                    return Err(e);
                }
                Err(_) => {
                    self.done = true;
                    return Err(std::io::Error::other("the gzip inflater stopped before the end of the stream"));
                }
            }
        }
        Ok(&self.buf[self.pos..])
    }

    fn consume(&mut self, n: usize) {
        self.pos = (self.pos + n).min(self.buf.len());
    }
}

fn bio_utf8(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

fn bio_is_gz(path: &str) -> bool {
    path.ends_with(".gz")
}

// Writes a file, gzip-compressed when the path ends in ".gz". Output goes to
// "<path>.part" and is renamed to <path> by close(); dropping an unclosed
// writer removes the part file, so a write that fails part-way leaves no file
// that looks complete.
enum BioOut {
    Plain(std::io::BufWriter<std::fs::File>),
    Gz(flate2::write::GzEncoder<std::io::BufWriter<std::fs::File>>),
}

pub struct BioWriter {
    path: String,
    part: String,
    out: Option<BioOut>,
}

impl BioWriter {
    pub fn create(path: &str) -> BioWriter {
        let part = format!("{}.part", path);
        let f = std::fs::File::create(&part).unwrap_or_else(|e| panic!("cannot open {} for writing: {}", part, e));
        let buf = std::io::BufWriter::with_capacity(1 << 18, f);
        let out = if bio_is_gz(path) {
            BioOut::Gz(flate2::write::GzEncoder::new(buf, flate2::Compression::default()))
        } else {
            BioOut::Plain(buf)
        };
        BioWriter { path: path.to_string(), part, out: Some(out) }
    }

    pub fn write(&mut self, data: &[u8]) {
        use std::io::Write;
        let r = match self.out.as_mut().unwrap() {
            BioOut::Plain(w) => w.write_all(data),
            BioOut::Gz(w) => w.write_all(data),
        };
        r.unwrap_or_else(|e| panic!("write error in {}: {}", self.path, e));
    }

    pub fn close(mut self) {
        use std::io::Write;
        let r = match self.out.take().unwrap() {
            BioOut::Plain(mut w) => w.flush(),
            BioOut::Gz(w) => w.finish().and_then(|mut b| b.flush()),
        };
        r.unwrap_or_else(|e| panic!("error closing {}: {}", self.path, e));
        std::fs::rename(&self.part, &self.path)
            .unwrap_or_else(|e| panic!("cannot rename {} to {}: {}", self.part, self.path, e));
    }
}

impl Drop for BioWriter {
    fn drop(&mut self) {
        if self.out.take().is_some() {
            let _ = std::fs::remove_file(&self.part);
        }
    }
}

// Appends text to a file, as a new gzip member when the path ends in ".gz".
pub fn bio_append_text(path: &String, text: &String) {
    use std::io::Write;
    let f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .unwrap_or_else(|e| panic!("cannot open {} for appending: {}", path, e));
    if bio_is_gz(path) {
        let mut g = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        g.write_all(text.as_bytes()).unwrap_or_else(|e| panic!("write error in {}: {}", path, e));
        g.finish().unwrap_or_else(|e| panic!("write error in {}: {}", path, e));
    } else {
        let mut f = f;
        f.write_all(text.as_bytes()).unwrap_or_else(|e| panic!("write error in {}: {}", path, e));
    }
}

// Creates an empty file, replacing any file at the path.
pub fn bio_reset_file(path: &String) {
    std::fs::File::create(path).unwrap_or_else(|e| panic!("cannot create {}: {}", path, e));
}

// Collects records and hands them to a sink in batches.
struct BioBatcher<'a, T, S: Fn(&Vec<T>) -> R, R> {
    sink: &'a S,
    limit: usize,
    weight: usize,
    batch: Vec<T>,
    _r: std::marker::PhantomData<R>,
}

impl<'a, T, S: Fn(&Vec<T>) -> R, R> BioBatcher<'a, T, S, R> {
    fn new(sink: &'a S, limit: usize) -> Self {
        BioBatcher { sink, limit, weight: 0, batch: Vec::new(), _r: std::marker::PhantomData }
    }
    fn add(&mut self, x: T, w: usize) {
        self.weight += w;
        self.batch.push(x);
        if self.weight >= self.limit {
            self.flush();
        }
    }
    fn flush(&mut self) {
        if self.batch.is_empty() {
            return;
        }
        (self.sink)(&self.batch);
        self.batch.clear();
        self.weight = 0;
    }
}

fn bio_trim_ws(s: &str) -> &str {
    s.trim_matches(|c| c == ' ' || c == '\t')
}

fn bio_parse_u64(r: &BioLineReader, s: &str, what: &str) -> u64 {
    if s.is_empty() || !s.bytes().all(|c| c.is_ascii_digit()) {
        r.fail(&format!("bad {} '{}'", what, s));
    }
    match s.parse::<u64>() {
        Ok(v) => v,
        Err(_) => r.fail(&format!("{} out of range '{}'", what, s)),
    }
}

// strtod's grammar, minus hex floats: what C++'s reader accepts.
fn bio_parse_score(r: &BioLineReader, s: &str) -> Option<f64> {
    if s == "." {
        return None;
    }
    let t = s.trim_start_matches(|c: char| c.is_ascii_whitespace());
    match t.parse::<f64>() {
        Ok(v) if !s.is_empty() && (v.is_finite() || t.to_ascii_lowercase().contains("inf") || t.to_ascii_lowercase().contains("nan")) => Some(v),
        _ => r.fail(&format!("bad score '{}'", s)),
    }
}

// printf's "%.{prec}g".
fn bio_fmt_g(v: f64, prec: usize) -> String {
    if v == 0.0 {
        return if v.is_sign_negative() { "-0".to_string() } else { "0".to_string() };
    }
    if !v.is_finite() {
        return if v.is_nan() { "nan".to_string() } else if v > 0.0 { "inf".to_string() } else { "-inf".to_string() };
    }
    let p = prec.max(1);
    let e = format!("{:.*e}", p - 1, v);
    let (mant, exp) = e.split_at(e.find('e').unwrap());
    let x: i64 = exp[1..].parse().unwrap();
    let strip = |s: String| -> String {
        if s.contains('.') {
            let t = s.trim_end_matches('0');
            t.trim_end_matches('.').to_string()
        } else {
            s
        }
    };
    if (x as i64) < (p as i64) && x >= -4 {
        let decimals = (p as i64 - 1 - x).max(0) as usize;
        strip(format!("{:.*}", decimals, v))
    } else {
        let m = strip(mant.to_string());
        let sign = if x < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", m, sign, x.abs())
    }
}

// An integral value as an integer, else the shortest decimal text that
// reads back as the same double.
fn bio_format_score(v: Option<f64>) -> String {
    let v = match v {
        None => return ".".to_string(),
        Some(v) => v,
    };
    if v.is_finite() && v.abs() < 1e15 && v == (v as i64) as f64 {
        return (v as i64).to_string();
    }
    let mut s = String::new();
    for prec in 1..=17 {
        s = bio_fmt_g(v, prec);
        if s.parse::<f64>().map(|x| x == v).unwrap_or(false) {
            break;
        }
    }
    s
}

fn bio_parse_strand(r: &BioLineReader, s: &str) -> BioStrand {
    match s {
        "+" => BioStrand::Plus,
        "-" => BioStrand::Minus,
        "." => BioStrand::Unstranded,
        "?" => BioStrand::Unknown,
        _ => r.fail(&format!("bad strand '{}'", s)),
    }
}

fn bio_strand_char(s: BioStrand) -> char {
    match s {
        BioStrand::Plus => '+',
        BioStrand::Minus => '-',
        BioStrand::Unstranded => '.',
        BioStrand::Unknown => '?',
    }
}


// FASTA
// -----

// Whether a line holds a space or a tab. Blocks of 32 bytes are tested
// whole, without an early exit, so the test compiles to vector compares.
fn bio_has_blank(b: &[u8]) -> bool {
    let mut blocks = b.chunks_exact(32);
    for c in &mut blocks {
        let mut hit = 0u8;
        for &x in c {
            hit |= (x == b' ') as u8 | (x == b'\t') as u8;
        }
        if hit != 0 {
            return true;
        }
    }
    blocks.remainder().iter().any(|&x| x == b' ' || x == b'\t')
}

pub fn bio_stream_fasta<S: Fn(&Vec<BioFastaEntry>) -> R, R>(path: &String, sink: S) {
    let mut r = BioLineReader::open(path);
    let mut out = BioBatcher::new(&sink, BIO_BATCH_BYTES);
    let mut line = Vec::new();
    let mut cur: Option<BioFastaEntry> = None;
    while r.next(&mut line) {
        if line.is_empty() || line[0] == b';' {
            continue;
        }
        if line[0] == b'>' {
            if let Some(e) = cur.take() {
                let w = e.header.len() + e.bases.len();
                out.add(e, w);
            }
            cur = Some(BioSequence { header: bio_utf8(&line[1..]), bases: Vec::new() });
            continue;
        }
        match cur.as_mut() {
            None => r.fail("sequence before the first '>' header"),
            Some(e) => {
                if bio_has_blank(&line) {
                    e.bases.extend(line.iter().filter(|&&c| c != b' ' && c != b'\t'));
                } else {
                    e.bases.extend_from_slice(&line);
                }
            }
        }
    }
    if let Some(e) = cur.take() {
        let w = e.header.len() + e.bases.len();
        out.add(e, w);
    }
    out.flush();
}

const BIO_FASTA_WIDTH: usize = 60;

fn bio_fasta_text(s: &mut Vec<u8>, e: &BioFastaEntry) {
    s.push(b'>');
    s.extend_from_slice(e.header.as_bytes());
    s.push(b'\n');
    for chunk in e.bases.chunks(BIO_FASTA_WIDTH) {
        s.extend_from_slice(chunk);
        s.push(b'\n');
    }
}

fn bio_bytes_string(v: Vec<u8>) -> String {
    String::from_utf8(v).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned())
}

pub fn bio_format_fasta(batch: &Vec<BioFastaEntry>) -> String {
    let mut s = Vec::new();
    for e in batch {
        bio_fasta_text(&mut s, e);
    }
    bio_bytes_string(s)
}


// FASTQ
// -----

// The next four-line record, or None at a clean end of file.
fn bio_next_fastq(r: &mut BioLineReader, h: &mut Vec<u8>, s: &mut Vec<u8>, p: &mut Vec<u8>, q: &mut Vec<u8>) -> Option<BioFastqEntry> {
    loop {
        if !r.next(h) {
            return None;
        }
        if !h.is_empty() {
            break;
        }
    }
    if h[0] != b'@' {
        r.fail("expected a FASTQ header starting with '@'");
    }
    if !r.next(s) {
        r.fail("truncated FASTQ record");
    }
    if !r.next(p) || p.is_empty() || p[0] != b'+' {
        r.fail("expected a '+' line");
    }
    if !r.next(q) {
        r.fail("truncated FASTQ record");
    }
    if q.len() != s.len() {
        r.fail("quality length differs from sequence length");
    }
    Some(BioSequence { header: bio_utf8(&h[1..]), bases: BioBaseQual { bases: std::mem::take(s), quals: std::mem::take(q) } })
}

fn bio_fastq_text(s: &mut Vec<u8>, e: &BioFastqEntry) {
    if e.bases.bases.len() != e.bases.quals.len() {
        panic!("FASTQ entry '{}': quality length differs from sequence length", e.header);
    }
    s.push(b'@');
    s.extend_from_slice(e.header.as_bytes());
    s.push(b'\n');
    s.extend_from_slice(&e.bases.bases);
    s.extend_from_slice(b"\n+\n");
    s.extend_from_slice(&e.bases.quals);
    s.push(b'\n');
}

pub fn bio_format_fastq(batch: &Vec<BioFastqEntry>) -> String {
    let mut s = Vec::new();
    for e in batch {
        bio_fastq_text(&mut s, e);
    }
    bio_bytes_string(s)
}

pub fn bio_stream_fastq<S: Fn(&Vec<BioFastqEntry>) -> R, R>(path: &String, sink: S) {
    let mut r = BioLineReader::open(path);
    let mut out = BioBatcher::new(&sink, BIO_BATCH_BYTES);
    let (mut h, mut s, mut p, mut q) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    while let Some(e) = bio_next_fastq(&mut r, &mut h, &mut s, &mut p, &mut q) {
        let w = e.header.len() + 2 * e.bases.bases.len();
        out.add(e, w);
    }
    out.flush();
}

pub fn bio_stream_fastq_pairs<S: Fn(&Vec<(BioFastqEntry, BioFastqEntry)>) -> R, R>(path1: &String, path2: &String, sink: S) {
    let mut r1 = BioLineReader::open(path1);
    let mut r2 = BioLineReader::open(path2);
    let mut out = BioBatcher::new(&sink, BIO_BATCH_BYTES);
    let (mut h, mut s, mut p, mut q) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    loop {
        let a = bio_next_fastq(&mut r1, &mut h, &mut s, &mut p, &mut q);
        let b = bio_next_fastq(&mut r2, &mut h, &mut s, &mut p, &mut q);
        match (a, b) {
            (None, None) => break,
            (Some(a), Some(b)) => {
                let w = a.header.len() + b.header.len() + 2 * (a.bases.bases.len() + b.bases.bases.len());
                out.add((a, b), w);
            }
            _ => panic!("mate files {} and {} hold different numbers of reads", path1, path2),
        }
    }
    out.flush();
}


// GFF3 and GTF
// ------------

type BioAttrs = Vec<(String, Vec<String>)>;

fn bio_attr_add(attrs: &mut BioAttrs, key: String, values: Vec<String>) {
    for (k, vs) in attrs.iter_mut() {
        if *k == key {
            vs.extend(values);
            return;
        }
    }
    attrs.push((key, values));
}

fn bio_attr_get<'a>(attrs: &'a BioAttrs, key: &str) -> Option<&'a Vec<String>> {
    attrs.iter().find(|(k, _)| k == key).map(|(_, v)| v)
}

fn bio_hex(c: u8) -> Option<u8> {
    match c {
        b'0'..=b'9' => Some(c - b'0'),
        b'a'..=b'f' => Some(c - b'a' + 10),
        b'A'..=b'F' => Some(c - b'A' + 10),
        _ => None,
    }
}

fn bio_percent_decode_with<F: Fn(&str) -> String>(s: &str, fail: F) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] != b'%' {
            out.push(b[i]);
            i += 1;
            continue;
        }
        if i + 2 >= b.len() {
            panic!("{}", fail(&format!("truncated percent escape in '{}'", s)));
        }
        match (bio_hex(b[i + 1]), bio_hex(b[i + 2])) {
            (Some(hi), Some(lo)) => out.push(hi * 16 + lo),
            _ => panic!("{}", fail(&format!("bad percent escape in '{}'", s))),
        }
        i += 3;
    }
    bio_bytes_string(out)
}

fn bio_reader_msg(r: &BioLineReader, m: &str) -> String {
    format!("{}:{}: {}", r.path, r.lineno, m)
}

fn bio_percent_encode(s: &str) -> String {
    const HEX: &[u8] = b"0123456789ABCDEF";
    let mut out = Vec::with_capacity(s.len());
    for &c in s.as_bytes() {
        if c < 0x20 || c == 0x7f || c == b'%' || c == b';' || c == b'=' || c == b'&' || c == b',' {
            out.push(b'%');
            out.push(HEX[(c >> 4) as usize]);
            out.push(HEX[(c & 15) as usize]);
        } else {
            out.push(c);
        }
    }
    bio_bytes_string(out)
}

fn bio_attrs_parse(col: &str) -> BioAttrs {
    let fail = |m: &str| format!("attributes: {}", m);
    let mut attrs = BioAttrs::new();
    if col.is_empty() || col == "." {
        return attrs;
    }
    for field in col.split(';') {
        let f = bio_trim_ws(field);
        if f.is_empty() {
            continue;
        }
        let eq = match f.find('=') {
            Some(eq) => eq,
            None => panic!("{}", fail(&format!("attribute without '=': '{}'", f))),
        };
        let values = f[eq + 1..].split(',').map(|v| bio_percent_decode_with(v, fail)).collect();
        bio_attr_add(&mut attrs, bio_percent_decode_with(&f[..eq], fail), values);
    }
    attrs
}

fn bio_attrs_text(attrs: &BioAttrs) -> String {
    let mut s = String::new();
    for (k, vs) in attrs {
        if !s.is_empty() {
            s.push(';');
        }
        s.push_str(&bio_percent_encode(k));
        s.push('=');
        for (i, v) in vs.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&bio_percent_encode(v));
        }
    }
    s
}

pub fn bio_attr_values(key: &String, attrs: &String) -> Vec<String> {
    bio_attr_get(&bio_attrs_parse(attrs), key).cloned().unwrap_or_default()
}

// Splits a GFF3 attribute column into its ID, its Parents and the rest of
// the column, which is kept as written. Escapes are checked everywhere.
fn bio_split_gff3_attrs(r: &BioLineReader, col: &str, id: &mut String, parents: &mut Vec<String>) -> String {
    let mut rest = String::new();
    if col == "." {
        return rest;
    }
    let fail = |m: &str| bio_reader_msg(r, m);
    let mut have_id = false;
    for field in col.split(';') {
        let f = bio_trim_ws(field);
        if f.is_empty() {
            continue;
        }
        let eq = match f.find('=') {
            Some(eq) => eq,
            None => r.fail(&format!("attribute without '=': '{}'", f)),
        };
        let key = bio_percent_decode_with(&f[..eq], fail);
        let values: Vec<String> = f[eq + 1..].split(',').map(|v| bio_percent_decode_with(v, fail)).collect();
        if key == "ID" {
            if have_id || values.len() != 1 {
                r.fail("attribute ID has several values");
            }
            *id = values[0].clone();
            have_id = true;
        } else if key == "Parent" {
            parents.extend(values);
        } else {
            if !rest.is_empty() {
                rest.push(';');
            }
            rest.push_str(f);
        }
    }
    rest
}

// GTF attributes: `key "value"; key value;`, a key repeated for several
// values.
fn bio_parse_gtf_attrs(r: &BioLineReader, col: &str) -> BioAttrs {
    let mut attrs = BioAttrs::new();
    let b = col.as_bytes();
    let n = b.len();
    let mut i = 0;
    let skip = |i: &mut usize| {
        while *i < n && (b[*i] == b' ' || b[*i] == b'\t') {
            *i += 1;
        }
    };
    loop {
        skip(&mut i);
        while i < n && b[i] == b';' {
            i += 1;
            skip(&mut i);
        }
        if i >= n {
            break;
        }
        let k0 = i;
        while i < n && b[i] != b' ' && b[i] != b'\t' && b[i] != b';' {
            i += 1;
        }
        let key = col[k0..i].to_string();
        skip(&mut i);
        let value;
        if i < n && b[i] == b'"' {
            let close = match col[i + 1..].find('"') {
                Some(c) => i + 1 + c,
                None => r.fail("unterminated quote in attributes"),
            };
            value = col[i + 1..close].to_string();
            i = close + 1;
        } else {
            let v0 = i;
            while i < n && b[i] != b';' {
                i += 1;
            }
            value = bio_trim_ws(&col[v0..i]).to_string();
        }
        if key.is_empty() {
            r.fail("attribute without a key");
        }
        bio_attr_add(&mut attrs, key.clone(), vec![value]);
        skip(&mut i);
        if i < n && b[i] != b';' {
            r.fail(&format!("expected ';' after attribute {}", key));
        }
    }
    attrs
}

fn bio_phase_char(p: BioPhase) -> char {
    match p {
        BioPhase::Phase0 => '0',
        BioPhase::Phase1 => '1',
        BioPhase::Phase2 => '2',
        BioPhase::PhaseNone => '.',
    }
}

struct BioGffLine {
    seqid: String,
    origin: String,
    kind: String,
    part: BioSegment,
    score: Option<f64>,
    phase: BioPhase,
}

fn bio_parse_gff_columns(r: &BioLineReader, line: &str) -> (BioGffLine, String) {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() != 9 {
        r.fail(&format!("expected 9 tab-separated columns, found {}", cols.len()));
    }
    let start = bio_parse_u64(r, cols[3], "start");
    let end = bio_parse_u64(r, cols[4], "end");
    if start == 0 {
        r.fail("start is 0; GFF coordinates are 1-based");
    }
    if start > end + 1 {
        r.fail("start is past end");
    }
    let part = BioSegment { start: start - 1, end, strand: bio_parse_strand(r, cols[6]) };
    let score = bio_parse_score(r, cols[5]);
    let phase = match cols[7] {
        "." => BioPhase::PhaseNone,
        "0" => BioPhase::Phase0,
        "1" => BioPhase::Phase1,
        "2" => BioPhase::Phase2,
        x => r.fail(&format!("bad phase '{}'", x)),
    };
    (
        BioGffLine {
            seqid: cols[0].to_string(),
            origin: if cols[1] == "." { String::new() } else { cols[1].to_string() },
            kind: cols[2].to_string(),
            part,
            score,
            phase,
        },
        cols[8].to_string(),
    )
}

// A node under construction: the node, one segment, score and phase per
// line, and where it sits.
#[derive(Clone)]
struct BioGffBuild {
    seqid: String,
    node: BioGffNode,
    parts: Vec<BioSegment>,
    parents: Vec<String>,
    children: Vec<usize>,
}

fn bio_built(b: &BioGffBuild) -> BioFeature<BioGffNode> {
    bio_built_feature(&b.node, &b.parts)
}

// Builds the trees of a set of nodes whose Parent references resolve among
// themselves, in order of first appearance, copying a node under each of its
// parents.
fn bio_gff_trees(
    path: &str,
    nodes: &mut Vec<BioGffBuild>,
    by_id: &std::collections::HashMap<String, usize>,
    roots: &mut Vec<BioGffRoot>,
) {
    let mut tops = Vec::new();
    for i in 0..nodes.len() {
        if nodes[i].parents.is_empty() {
            tops.push(i);
            continue;
        }
        for p in nodes[i].parents.clone() {
            let j = match by_id.get(&p) {
                Some(&j) => j,
                None => panic!("{}: Parent '{}' of a {} feature is not defined", path, p, nodes[i].node.kind),
            };
            if nodes[j].seqid != nodes[i].seqid {
                panic!("{}: feature tree '{}' spans several sequences", path, p);
            }
            nodes[j].children.push(i);
        }
    }
    let mut state = vec![0u8; nodes.len()];
    fn build(path: &str, nodes: &Vec<BioGffBuild>, state: &mut Vec<u8>, i: usize) -> BioFeature<BioGffNode> {
        if state[i] == 1 {
            panic!("{}: Parent references form a cycle", path);
        }
        state[i] = 1;
        let mut f = bio_built(&nodes[i]);
        for &c in &nodes[i].children {
            f.children.push(build(path, nodes, state, c));
        }
        state[i] = 2;
        f
    }
    for &t in &tops {
        let f = build(path, nodes, &mut state, t);
        roots.push(BioRoot { seqid: nodes[t].seqid.clone(), feature: f });
    }
    if state.iter().any(|&s| s == 0) {
        panic!("{}: Parent references form a cycle", path);
    }
}

// Reads GFF3. Trees are complete at a "###" line or at the end of the file;
// a file without "###" lines is held in memory until its end.
pub fn bio_stream_gff3<S: Fn(&Vec<BioGffRoot>) -> R, R>(path: &String, sink: S) {
    let mut r = BioLineReader::open(path);
    let mut out = BioBatcher::new(&sink, BIO_BATCH_TREES);
    let mut nodes: Vec<BioGffBuild> = Vec::new();
    let mut by_id: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let flush = |nodes: &mut Vec<BioGffBuild>, by_id: &mut std::collections::HashMap<String, usize>, out: &mut BioBatcher<BioGffRoot, S, R>| {
        let mut roots = Vec::new();
        bio_gff_trees(path, nodes, by_id, &mut roots);
        for root in roots {
            out.add(root, 1);
        }
        nodes.clear();
        by_id.clear();
    };
    let mut raw = Vec::new();
    while r.next(&mut raw) {
        if raw.is_empty() {
            continue;
        }
        let line = bio_utf8(&raw);
        if line.starts_with('#') {
            if line == "###" {
                flush(&mut nodes, &mut by_id, &mut out);
            }
            if line.starts_with("##FASTA") {
                break;
            }
            continue;
        }
        let (g, attr_col) = bio_parse_gff_columns(&r, &line);
        let mut id = String::new();
        let mut parents = Vec::new();
        let attrs = bio_split_gff3_attrs(&r, &attr_col, &mut id, &mut parents);
        if !id.is_empty() {
            if let Some(&j) = by_id.get(&id) {
                let b = &mut nodes[j];
                if b.seqid != g.seqid || b.node.kind != g.kind {
                    r.fail(&format!("ID '{}' is reused by a different feature", id));
                }
                b.parts.push(g.part);
                b.node.scores.push(g.score);
                b.node.phases.push(g.phase);
                continue;
            }
            by_id.insert(id.clone(), nodes.len());
        }
        nodes.push(BioGffBuild {
            seqid: g.seqid,
            node: BioGffNode {
                kind: g.kind,
                origin: g.origin,
                id,
                attrs,
                scores: vec![g.score],
                phases: vec![g.phase],
                created: false,
            },
            parts: vec![g.part],
            parents,
            children: Vec::new(),
        });
    }
    flush(&mut nodes, &mut by_id, &mut out);
    out.flush();
}

fn bio_check_parallel(n: &BioGffNode, nparts: usize) {
    if (!n.scores.is_empty() && n.scores.len() != nparts) || (!n.phases.is_empty() && n.phases.len() != nparts) {
        panic!("a {} node's scores and phases do not match its parts", n.kind);
    }
}

fn bio_node_score(n: &BioGffNode, i: usize) -> String {
    if n.scores.is_empty() { ".".to_string() } else { bio_format_score(n.scores[i]) }
}

fn bio_node_phase(n: &BioGffNode, i: usize) -> char {
    if n.phases.is_empty() { '.' } else { bio_phase_char(n.phases[i]) }
}

fn bio_feature_line(s: &mut String, seqid: &str, n: &BioGffNode, p: &BioSegment, i: usize, attrs: &str) {
    use std::fmt::Write;
    let origin = if n.origin.is_empty() { "." } else { &n.origin };
    let _ = write!(
        s,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        seqid,
        origin,
        n.kind,
        p.start + 1,
        p.end,
        bio_node_score(n, i),
        bio_strand_char(p.strand),
        bio_node_phase(n, i),
        attrs
    );
}

// One tree as GFF3 lines, followed by "###". A node copied under several
// parents is written once, listing all of them. A node without an ID that
// has children gets `prefix` followed by a count kept in `fresh`.
fn bio_gff3_text(s: &mut String, root: &BioGffRoot, prefix: &str, fresh: &mut usize) {
    let mut parents_of: std::collections::HashMap<String, Vec<String>> = std::collections::HashMap::new();
    // (feature, key): key is the node's id, or "\x01" + its parent for a
    // leaf without an id.
    let mut order: Vec<(&BioFeature<BioGffNode>, String)> = Vec::new();
    fn visit<'a>(
        f: &'a BioFeature<BioGffNode>,
        parent: Option<&str>,
        prefix: &str,
        fresh: &mut usize,
        parents_of: &mut std::collections::HashMap<String, Vec<String>>,
        order: &mut Vec<(&'a BioFeature<BioGffNode>, String)>,
    ) {
        let mut id = f.node.id.clone();
        if id.is_empty() && !f.children.is_empty() {
            *fresh += 1;
            id = format!("{}{}", prefix, fresh);
        }
        let seen = !id.is_empty() && parents_of.contains_key(&id);
        if !id.is_empty() {
            let ps = parents_of.entry(id.clone()).or_default();
            if let Some(p) = parent {
                ps.push(p.to_string());
            }
        }
        if !seen {
            order.push((f, id.clone()));
        }
        if id.is_empty() {
            if let Some(p) = parent {
                order.last_mut().unwrap().1 = format!("\x01{}", p);
            }
        }
        if !seen {
            for c in &f.children {
                visit(c, if id.is_empty() { None } else { Some(&id) }, prefix, fresh, parents_of, order);
            }
        }
    }
    visit(&root.feature, None, prefix, fresh, &mut parents_of, &mut order);
    for (f, key) in &order {
        let n = &f.node;
        bio_check_parallel(n, f.parts.len());
        let mut attrs = String::new();
        let mut append = |a: &str| {
            if a.is_empty() {
                return;
            }
            if !attrs.is_empty() {
                attrs.push(';');
            }
            attrs.push_str(a);
        };
        if let Some(p) = key.strip_prefix('\x01') {
            append(&format!("Parent={}", bio_percent_encode(p)));
        } else if !key.is_empty() {
            append(&format!("ID={}", bio_percent_encode(key)));
            let ps = &parents_of[key];
            if !ps.is_empty() {
                let enc: Vec<String> = ps.iter().map(|p| bio_percent_encode(p)).collect();
                append(&format!("Parent={}", enc.join(",")));
            }
        }
        append(&n.attrs);
        if attrs.is_empty() {
            attrs = ".".to_string();
        }
        for (i, p) in f.parts.iter().enumerate() {
            bio_feature_line(s, &root.seqid, n, p, i, &attrs);
        }
    }
    s.push_str("###\n");
}

pub fn bio_format_gff3(offset: i64, batch: &Vec<BioGffRoot>) -> String {
    let mut s = if offset == 0 { "##gff-version 3\n".to_string() } else { String::new() };
    let prefix = format!("bio_node_{}_", offset);
    let mut fresh = 0usize;
    for root in batch {
        bio_gff3_text(&mut s, root, &prefix, &mut fresh);
    }
    s
}

// A gene assembled from GTF lines.
#[derive(Default)]
struct BioGtfGene {
    seqid: String,
    gene: Option<BioGffBuild>,
    tx_order: Vec<String>,
    tx: std::collections::HashMap<String, BioGffBuild>,
    tx_kids: std::collections::HashMap<String, Vec<BioGffBuild>>,
    gene_kids: Vec<BioGffBuild>,
}

// The span covering a set of parts: strand shared by all, else Unknown.
fn bio_envelope(parts: &[BioSegment]) -> BioSegment {
    let s = parts.iter().map(|p| p.start).min().unwrap_or(0);
    let e = parts.iter().map(|p| p.end).max().unwrap_or(0);
    let mut st = if parts.is_empty() { BioStrand::Unknown } else { parts[0].strand };
    if parts.iter().any(|p| p.strand != st) {
        st = BioStrand::Unknown;
    }
    BioSegment { start: s, end: e, strand: st }
}

fn bio_created(kind: &str, id: &str, attrs: &BioAttrs, origin: &str, parts: &[BioSegment]) -> BioFeature<BioGffNode> {
    let n = BioGffNode {
        kind: kind.to_string(),
        origin: origin.to_string(),
        id: id.to_string(),
        attrs: bio_attrs_text(attrs),
        scores: Vec::new(),
        phases: Vec::new(),
        created: true,
    };
    bio_feature(n, &[bio_envelope(parts)])
}

fn bio_gtf_tree(gene_id: &str, g: &BioGtfGene) -> BioGffRoot {
    let mut gene_parts = Vec::new();
    let mut origin = String::new();
    let mut txs = Vec::new();
    let empty = Vec::new();
    for t in &g.tx_order {
        let kids = g.tx_kids.get(t).unwrap_or(&empty);
        let mut tx_parts = Vec::new();
        for k in kids {
            tx_parts.extend_from_slice(&k.parts);
            if origin.is_empty() {
                origin = k.node.origin.clone();
            }
        }
        let mut tf = match g.tx.get(t) {
            Some(b) => {
                if origin.is_empty() {
                    origin = b.node.origin.clone();
                }
                bio_built(b)
            }
            None => {
                let a: BioAttrs = vec![
                    ("gene_id".to_string(), vec![gene_id.to_string()]),
                    ("transcript_id".to_string(), vec![t.clone()]),
                ];
                let o = if kids.is_empty() { String::new() } else { kids[0].node.origin.clone() };
                bio_created("transcript", t, &a, &o, &tx_parts)
            }
        };
        for k in kids {
            tf.children.push(bio_built(k));
        }
        gene_parts.extend_from_slice(&tf.parts);
        txs.push(tf);
    }
    for k in &g.gene_kids {
        gene_parts.extend_from_slice(&k.parts);
        if origin.is_empty() {
            origin = k.node.origin.clone();
        }
    }
    let mut gf = match &g.gene {
        Some(b) => bio_built(b),
        None => {
            let a: BioAttrs = vec![("gene_id".to_string(), vec![gene_id.to_string()])];
            bio_created("gene", gene_id, &a, &origin, &gene_parts)
        }
    };
    gf.children.extend(txs);
    for k in &g.gene_kids {
        gf.children.push(bio_built(k));
    }
    BioRoot { seqid: g.seqid.clone(), feature: gf }
}


// Reads GTF. A gene's lines must be contiguous; a gene_id that reappears
// after another gene has started is an error.
pub fn bio_stream_gtf<S: Fn(&Vec<BioGffRoot>) -> R, R>(path: &String, sink: S) {
    bio_stream_gtf_with(path, |_: &str, _: &BioAttrs| true, sink)
}

// A GTF stream keeping only lines on the named sequences and, when
// `check_tx` holds, with a non-empty transcript_id (nf-core's gtffilter).
pub fn bio_stream_gtf_filtered<S: Fn(&Vec<BioGffRoot>) -> R, R>(names: &Vec<String>, check_tx: bool, path: &String, sink: S) {
    let keep: std::collections::HashSet<&str> = names.iter().map(|s| s.as_str()).collect();
    bio_stream_gtf_with(
        path,
        |seqid: &str, attrs: &BioAttrs| {
            keep.contains(seqid)
                && (!check_tx || bio_attr_get(attrs, "transcript_id").map_or(false, |v| v.iter().any(|t| !t.is_empty())))
        },
        sink,
    )
}

fn bio_stream_gtf_with<K: Fn(&str, &BioAttrs) -> bool, S: Fn(&Vec<BioGffRoot>) -> R, R>(path: &String, keep: K, sink: S) {
    let r = std::cell::RefCell::new(BioLineReader::open(path));
    let mut raw = Vec::new();
    let lines = std::iter::from_fn(|| loop {
        let mut rr = r.borrow_mut();
        if !rr.next(&mut raw) {
            return None;
        }
        if raw.is_empty() || raw[0] == b'#' {
            continue;
        }
        let line = bio_utf8(&raw);
        let (g, attr_col) = bio_parse_gff_columns(&rr, &line);
        let attrs = bio_parse_gtf_attrs(&rr, &attr_col);
        if keep(&g.seqid, &attrs) {
            return Some((g, attrs));
        }
    });
    bio_gtf_build(lines, |m: &str| bio_reader_msg(&r.borrow(), m), sink)
}

// Builds a GTF's trees from its lines. A gene's lines must be contiguous on
// one sequence; a sequence's genes are held until the sequence ends, so a
// GTF sorted by position, whose genes interleave, reads like one grouped by
// gene. `msg` places an error message.
fn bio_gtf_build<I, M, S, R>(lines: I, msg: M, sink: S)
where
    I: Iterator<Item = (BioGffLine, BioAttrs)>,
    M: Fn(&str) -> String,
    S: Fn(&Vec<BioGffRoot>) -> R,
{
    let fail = |m: &str| -> ! { panic!("{}", msg(m)) };
    let one = |attrs: &BioAttrs, key: &str| -> Option<String> {
        match bio_attr_get(attrs, key) {
            None => None,
            Some(vs) if vs.is_empty() => None,
            Some(vs) if vs.len() > 1 => fail(&format!("attribute {} has several values", key)),
            Some(vs) => Some(vs[0].clone()),
        }
    };
    let mut out = BioBatcher::new(&sink, BIO_BATCH_TREES);
    let mut done: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut pending: Vec<(String, BioGtfGene)> = Vec::new();
    let mut slot: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut cur_seq: Option<String> = None;
    let flush = |pending: &mut Vec<(String, BioGtfGene)>,
                 slot: &mut std::collections::HashMap<String, usize>,
                 done: &mut std::collections::HashSet<String>,
                 out: &mut BioBatcher<BioGffRoot, S, R>| {
        for (id, gene) in pending.drain(..) {
            out.add(bio_gtf_tree(&id, &gene), 1);
            done.insert(id);
        }
        slot.clear();
    };
    for (g, attrs) in lines {
        let gene_id = one(&attrs, "gene_id");
        let tx_id = one(&attrs, "transcript_id");
        let gene_id = match gene_id {
            Some(x) if !x.is_empty() => x,
            _ => fail("GTF line without a gene_id"),
        };
        if cur_seq.as_deref() != Some(g.seqid.as_str()) {
            if slot.contains_key(&gene_id) {
                fail(&format!("gene {} spans several sequences", gene_id));
            }
            flush(&mut pending, &mut slot, &mut done, &mut out);
            cur_seq = Some(g.seqid.clone());
        }
        let k = match slot.get(&gene_id) {
            Some(&k) => k,
            None => {
                if done.contains(&gene_id) {
                    fail(&format!("lines of gene {} are not contiguous", gene_id));
                }
                pending.push((gene_id.clone(), BioGtfGene { seqid: g.seqid.clone(), ..Default::default() }));
                slot.insert(gene_id.clone(), pending.len() - 1);
                pending.len() - 1
            }
        };
        let c = &mut pending[k].1;
        let is_gene = g.kind == "gene";
        let is_tx = g.kind == "transcript";
        let mut id = String::new();
        if is_gene {
            id = gene_id.clone();
        }
        if is_tx {
            if let Some(t) = &tx_id {
                id = t.clone();
            }
        }
        let b = BioGffBuild {
            seqid: g.seqid.clone(),
            node: BioGffNode {
                kind: g.kind.clone(),
                origin: g.origin.clone(),
                id,
                attrs: bio_attrs_text(&attrs),
                scores: vec![g.score],
                phases: vec![g.phase],
                created: false,
            },
            parts: vec![g.part],
            parents: Vec::new(),
            children: Vec::new(),
        };
        match &tx_id {
            _ if is_gene => {
                if c.gene.is_some() {
                    fail(&format!("gene {} has two gene lines", gene_id));
                }
                c.gene = Some(b);
            }
            Some(t) if is_tx && !t.is_empty() => {
                if c.tx.contains_key(t) {
                    fail(&format!("transcript {} has two transcript lines", t));
                }
                if !c.tx_kids.contains_key(t) {
                    c.tx_order.push(t.clone());
                }
                c.tx_kids.entry(t.clone()).or_default();
                c.tx.insert(t.clone(), b);
            }
            Some(t) if !t.is_empty() => {
                if !c.tx_kids.contains_key(t) {
                    c.tx_order.push(t.clone());
                }
                c.tx_kids.entry(t.clone()).or_default().push(b);
            }
            _ => c.gene_kids.push(b),
        }
    }
    flush(&mut pending, &mut slot, &mut done, &mut out);
    out.flush();
}

fn bio_gtf_attrs(attrs: &BioAttrs) -> String {
    let mut s = String::new();
    for (k, vs) in attrs {
        for v in vs {
            if !s.is_empty() {
                s.push(' ');
            }
            s.push_str(k);
            s.push_str(" \"");
            s.push_str(v);
            s.push_str("\";");
        }
    }
    s
}

// Writes GTF. Nodes the reader created are not written. A line's gene_id is
// its own gene_id attribute or its root's, and its transcript_id is its own
// transcript_id attribute or that of the root's child it lies under (for
// either, the node's id stands in for a missing attribute).
fn bio_gtf_text(s: &mut String, root: &BioGffRoot) {
    fn visit(s: &mut String, seqid: &str, f: &BioFeature<BioGffNode>, depth: usize, gene_id: Option<String>, tx_id: Option<String>) {
        let n = &f.node;
        let nattrs = bio_attrs_parse(&n.attrs);
        let nid = if n.id.is_empty() { None } else { Some(n.id.clone()) };
        let mut g = gene_id;
        let mut t = tx_id;
        match bio_attr_get(&nattrs, "gene_id") {
            Some(v) if !v.is_empty() => g = Some(v[0].clone()),
            _ if depth == 0 => g = nid.clone(),
            _ => {}
        }
        match bio_attr_get(&nattrs, "transcript_id") {
            Some(v) if !v.is_empty() => t = Some(v[0].clone()),
            _ if depth == 1 => t = nid.clone(),
            _ => {}
        }
        if !n.created {
            let gid = match &g {
                Some(x) => x.clone(),
                None => panic!("cannot write a {} feature to GTF without a gene_id", n.kind),
            };
            bio_check_parallel(n, f.parts.len());
            let mut a: BioAttrs = vec![("gene_id".to_string(), vec![gid])];
            if let Some(tt) = &t {
                if depth > 0 {
                    a.push(("transcript_id".to_string(), vec![tt.clone()]));
                }
            }
            for kv in &nattrs {
                if kv.0 != "gene_id" && kv.0 != "transcript_id" {
                    a.push(kv.clone());
                }
            }
            let attrs = bio_gtf_attrs(&a);
            for (i, p) in f.parts.iter().enumerate() {
                bio_feature_line(s, seqid, n, p, i, &attrs);
            }
        }
        for c in &f.children {
            visit(s, seqid, c, depth + 1, g.clone(), t.clone());
        }
    }
    visit(s, &root.seqid, &root.feature, 0, None, None);
}

pub fn bio_format_gtf(batch: &Vec<BioGffRoot>) -> String {
    let mut s = String::new();
    for root in batch {
        bio_gtf_text(&mut s, root);
    }
    s
}


// BED
// ---

fn bio_parse_bed_score(r: &BioLineReader, s: &str) -> f64 {
    if s == "." { f64::NAN } else { bio_parse_score(r, s).unwrap() }
}

fn bio_format_bed_score(v: f64) -> String {
    if v.is_nan() { ".".to_string() } else { bio_format_score(Some(v)) }
}

fn bio_parse_rgb(r: &BioLineReader, s: &str) -> (bool, u8, u8, u8) {
    if s == "0" {
        return (false, 0, 0, 0);
    }
    let xs: Vec<&str> = s.split(',').collect();
    if xs.len() != 3 {
        r.fail(&format!("bad itemRgb '{}'", s));
    }
    let mut c = [0u8; 3];
    for i in 0..3 {
        let v = bio_parse_u64(r, xs[i], "itemRgb");
        if v > 255 {
            r.fail(&format!("bad itemRgb '{}'", s));
        }
        c[i] = v as u8;
    }
    (true, c[0], c[1], c[2])
}

fn bio_format_rgb(c: &(bool, u8, u8, u8)) -> String {
    if !c.0 { "0".to_string() } else { format!("{},{},{}", c.1, c.2, c.3) }
}

fn bio_bed_columns<'a>(r: &BioLineReader, line: &'a str, n: usize) -> Vec<&'a str> {
    let cols: Vec<&str> = line.splitn(n + 1, '\t').take(n).collect();
    if cols.len() < n {
        r.fail(&format!("expected at least {} tab-separated columns", n));
    }
    cols
}

fn bio_bed_skip(line: &str) -> bool {
    line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser")
}

fn bio_bed3_of(r: &BioLineReader, cols: &[&str]) -> BioBed3 {
    let x = BioBed3 { chrom: cols[0].to_string(), start: bio_parse_u64(r, cols[1], "start"), end: bio_parse_u64(r, cols[2], "end") };
    if x.start > x.end {
        r.fail("start is past end");
    }
    x
}

fn bio_bed6_of(r: &BioLineReader, cols: &[&str]) -> BioBed6 {
    let b = bio_bed3_of(r, cols);
    BioBed6 {
        chrom: b.chrom,
        start: b.start,
        end: b.end,
        name: cols[3].to_string(),
        score: bio_parse_bed_score(r, cols[4]),
        strand: bio_parse_strand(r, cols[5]),
    }
}

fn bio_bed12_of(r: &BioLineReader, cols: &[&str]) -> BioBed12 {
    let b = bio_bed6_of(r, cols);
    let mut x = BioBed12 {
        chrom: b.chrom,
        start: b.start,
        end: b.end,
        name: b.name,
        score: b.score,
        strand: b.strand,
        thickStart: bio_parse_u64(r, cols[6], "thickStart"),
        thickEnd: bio_parse_u64(r, cols[7], "thickEnd"),
        rgb: bio_parse_rgb(r, cols[8]),
        blocks: Vec::new(),
    };
    if x.thickStart > x.thickEnd || x.thickStart < x.start || x.thickEnd > x.end {
        r.fail("thickStart and thickEnd must lie in order within the record");
    }
    let count = bio_parse_u64(r, cols[9], "blockCount");
    let list = |s: &str, what: &str| -> Vec<u64> {
        s.split(',').filter(|x| !x.is_empty()).map(|x| bio_parse_u64(r, x, what)).collect()
    };
    let sizes = list(cols[10], "blockSize");
    let starts = list(cols[11], "blockStart");
    if count == 0 || sizes.len() as u64 != count || starts.len() as u64 != count {
        r.fail("blockCount does not match blockSizes and blockStarts");
    }
    let mut prev = x.start;
    for i in 0..count as usize {
        let bs = x.start + starts[i];
        let be = bs + sizes[i];
        if if i == 0 { bs != x.start } else { bs < prev } {
            r.fail("blocks must start at chromStart, ascending and apart");
        }
        x.blocks.push((bs, be));
        prev = be;
    }
    if prev != x.end {
        r.fail("the last block must end at chromEnd");
    }
    x
}

fn bio_stream_bed_n<T, F: Fn(&BioLineReader, &[&str]) -> T, S: Fn(&Vec<T>) -> R, R>(path: &str, n: usize, of: F, sink: S) {
    let mut r = BioLineReader::open(path);
    let mut out = BioBatcher::new(&sink, BIO_BATCH_TREES);
    let mut raw = Vec::new();
    while r.next(&mut raw) {
        let line = bio_utf8(&raw);
        if bio_bed_skip(&line) {
            continue;
        }
        let cols = bio_bed_columns(&r, &line, n);
        let x = of(&r, &cols);
        out.add(x, 1);
    }
    out.flush();
}

pub fn bio_stream_bed3<S: Fn(&Vec<BioBed3>) -> R, R>(path: &String, sink: S) {
    bio_stream_bed_n(path, 3, bio_bed3_of, sink)
}

pub fn bio_stream_bed6<S: Fn(&Vec<BioBed6>) -> R, R>(path: &String, sink: S) {
    bio_stream_bed_n(path, 6, bio_bed6_of, sink)
}

pub fn bio_stream_bed12<S: Fn(&Vec<BioBed12>) -> R, R>(path: &String, sink: S) {
    bio_stream_bed_n(path, 12, bio_bed12_of, sink)
}

fn bio_bed6_cols(s: &mut String, chrom: &str, start: u64, end: u64, name: &str, score: f64, strand: BioStrand) {
    use std::fmt::Write;
    let _ = write!(s, "{}\t{}\t{}\t{}\t{}\t{}", chrom, start, end, name, bio_format_bed_score(score), bio_strand_char(strand));
}

pub fn bio_format_bed3(batch: &Vec<BioBed3>) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    for x in batch {
        let _ = writeln!(s, "{}\t{}\t{}", x.chrom, x.start, x.end);
    }
    s
}

pub fn bio_format_bed6(batch: &Vec<BioBed6>) -> String {
    let mut s = String::new();
    for x in batch {
        bio_bed6_cols(&mut s, &x.chrom, x.start, x.end, &x.name, x.score, x.strand);
        s.push('\n');
    }
    s
}

pub fn bio_format_bed12(batch: &Vec<BioBed12>) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    for x in batch {
        if x.blocks.is_empty() {
            panic!("a BED12 record needs at least one block");
        }
        bio_bed6_cols(&mut s, &x.chrom, x.start, x.end, &x.name, x.score, x.strand);
        let mut sizes = String::new();
        let mut starts = String::new();
        for &(a, b) in &x.blocks {
            let _ = write!(sizes, "{},", b - a);
            let _ = write!(starts, "{},", a - x.start);
        }
        let _ = writeln!(s, "\t{}\t{}\t{}\t{}\t{}\t{}", x.thickStart, x.thickEnd, bio_format_rgb(&x.rgb), x.blocks.len(), sizes, starts);
    }
    s
}


// Writers
// -------
//
// A writer receives a producer and hands it a sink that writes each batch.
// The producer gets the sink as a morloc function value; a Rust closure
// becomes one through rustmorloc's blanket impl.

pub type BioSink<A> = rustmorloc::FnValue1<Vec<A>, ()>;

fn bio_write_with<A: 'static, P, R>(path: &str, produce: P, format: fn(&Vec<A>) -> String)
where
    P: Fn(&BioSink<A>) -> R,
{
    let out = std::sync::Arc::new(std::sync::Mutex::new(Some(BioWriter::create(path))));
    let inner = out.clone();
    let sink: BioSink<A> = std::sync::Arc::new(move |xs: &Vec<A>| {
        let text = format(xs);
        inner.lock().unwrap_or_else(|p| p.into_inner()).as_mut().unwrap().write(text.as_bytes());
    });
    produce(&sink);
    drop(sink);
    let w = out.lock().unwrap_or_else(|p| p.into_inner()).take().unwrap();
    w.close();
}

pub fn bio_write_fasta<P: Fn(&BioSink<BioFastaEntry>) -> R, R>(path: &String, produce: P) {
    bio_write_with(path, produce, bio_format_fasta)
}

pub fn bio_write_fastq<P: Fn(&BioSink<BioFastqEntry>) -> R, R>(path: &String, produce: P) {
    bio_write_with(path, produce, bio_format_fastq)
}

pub fn bio_write_fastq_pairs<P: Fn(&BioSink<(BioFastqEntry, BioFastqEntry)>) -> R, R>(path1: &String, path2: &String, produce: P) {
    let out = std::sync::Arc::new(std::sync::Mutex::new(Some((BioWriter::create(path1), BioWriter::create(path2)))));
    let inner = out.clone();
    let sink: BioSink<(BioFastqEntry, BioFastqEntry)> = std::sync::Arc::new(move |xs: &Vec<(BioFastqEntry, BioFastqEntry)>| {
        let (mut s1, mut s2) = (Vec::new(), Vec::new());
        for (a, b) in xs {
            bio_fastq_text(&mut s1, a);
            bio_fastq_text(&mut s2, b);
        }
        let mut ws = inner.lock().unwrap_or_else(|p| p.into_inner());
        let ws = ws.as_mut().unwrap();
        ws.0.write(&s1);
        ws.1.write(&s2);
    });
    produce(&sink);
    drop(sink);
    let (w1, w2) = out.lock().unwrap_or_else(|p| p.into_inner()).take().unwrap();
    w1.close();
    w2.close();
}

pub fn bio_write_gff3<P: Fn(&BioSink<BioGffRoot>) -> R, R>(path: &String, produce: P) {
    let out = std::sync::Arc::new(std::sync::Mutex::new(Some((BioWriter::create(path), 0usize))));
    out.lock().unwrap_or_else(|p| p.into_inner()).as_mut().unwrap().0.write(b"##gff-version 3\n");
    let inner = out.clone();
    let sink: BioSink<BioGffRoot> = std::sync::Arc::new(move |xs: &Vec<BioGffRoot>| {
        let mut st = inner.lock().unwrap_or_else(|p| p.into_inner());
        let st = st.as_mut().unwrap();
        let mut s = String::new();
        for root in xs {
            bio_gff3_text(&mut s, root, "bio_node_", &mut st.1);
        }
        st.0.write(s.as_bytes());
    });
    produce(&sink);
    drop(sink);
    let (w, _) = out.lock().unwrap_or_else(|p| p.into_inner()).take().unwrap();
    w.close();
}

pub fn bio_write_gtf<P: Fn(&BioSink<BioGffRoot>) -> R, R>(path: &String, produce: P) {
    bio_write_with(path, produce, bio_format_gtf)
}

pub fn bio_write_bed3<P: Fn(&BioSink<BioBed3>) -> R, R>(path: &String, produce: P) {
    bio_write_with(path, produce, bio_format_bed3)
}

pub fn bio_write_bed6<P: Fn(&BioSink<BioBed6>) -> R, R>(path: &String, produce: P) {
    bio_write_with(path, produce, bio_format_bed6)
}

pub fn bio_write_bed12<P: Fn(&BioSink<BioBed12>) -> R, R>(path: &String, produce: P) {
    bio_write_with(path, produce, bio_format_bed12)
}


// Reference products
// ------------------

// The first word of each header in a FASTA file, in order.
pub fn bio_fasta_names(path: &String) -> Vec<String> {
    let mut r = BioLineReader::open(path);
    let mut line = Vec::new();
    let mut out = Vec::new();
    while r.next(&mut line) {
        if line.first() == Some(&b'>') {
            out.push(bio_utf8(&line[1..]).split_whitespace().next().unwrap_or("").to_string());
        }
    }
    out
}

// Write the spliced sequence of every transcript of `entries` as FASTA,
// reading the genome once, contig by contig; returns the transcript ids in
// the order written (genome order, then model order within a contig). A
// transcript on a sequence absent from the genome is an error.
pub fn bio_write_transcript_fasta(genome: &String, entries: &Vec<BioGffEntry>, out: &String) -> Vec<String> {
    let w = std::cell::RefCell::new(BioTxWriter::new(entries, out));
    bio_stream_fasta(genome, |batch: &Vec<BioFastaEntry>| w.borrow_mut().add(batch));
    w.into_inner().finish(genome)
}

// `bio_write_transcript_fasta` with the genome from a pull.
pub fn bio_write_transcript_fasta_pull<P>(entries: &Vec<BioGffEntry>, out: &String, pull: P) -> Vec<String>
where
    P: Fn() -> (bool, String, Vec<BioFastaEntry>),
{
    let mut w = BioTxWriter::new(entries, out);
    bio_drain_pull("writeTranscriptFastaPull", pull, |batch| w.add(batch));
    w.finish("the genome")
}

// Calls `f` on each batch of a pull until it is exhausted; a failed pull is
// an error.
fn bio_drain_pull<T, P, F>(what: &str, pull: P, mut f: F)
where
    P: Fn() -> (bool, String, Vec<T>),
    F: FnMut(&Vec<T>),
{
    loop {
        let (ok, msg, xs) = pull();
        if !ok {
            panic!("{}: {}", what, msg);
        }
        if xs.is_empty() {
            return;
        }
        f(&xs);
    }
}

struct BioTxWriter {
    names: Vec<String>,
    by_seq: std::collections::HashMap<String, Vec<BioTxModel>>,
    w: BioWriter,
    written: Vec<String>,
    found: std::collections::HashSet<String>,
}

impl BioTxWriter {
    fn new(entries: &Vec<BioGffEntry>, out: &String) -> BioTxWriter {
        let mut names: Vec<String> = Vec::new();
        let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
        for e in entries {
            let s = bio_entry_seqid(e);
            if seen.insert(s.clone()) {
                names.push(s);
            }
        }
        let models = bio_transcript_models(&names, entries);
        let n = models.len();
        let mut by_seq: std::collections::HashMap<String, Vec<BioTxModel>> = std::collections::HashMap::new();
        for m in models {
            by_seq.entry(m.seqid.clone()).or_default().push(m);
        }
        BioTxWriter { names, by_seq, w: BioWriter::create(out), written: Vec::with_capacity(n), found: std::collections::HashSet::new() }
    }

    fn add(&mut self, batch: &Vec<BioFastaEntry>) {
        for g in batch {
            let name = g.header.split_whitespace().next().unwrap_or("").to_string();
            if let Some(ms) = self.by_seq.get(name.as_str()) {
                let mut text = Vec::new();
                for m in ms {
                    let e = BioSequence { header: m.id.clone(), bases: bio_extract(&m.mapping, &g.bases) };
                    bio_fasta_text(&mut text, &e);
                    self.written.push(m.id.clone());
                }
                self.w.write(&text);
                self.found.insert(name);
            }
        }
    }

    fn finish(self, genome: &str) -> Vec<String> {
        for n in &self.names {
            if !self.found.contains(n) {
                panic!("writeTranscriptFasta: sequence {} of the annotation is not in {}", n, genome);
            }
        }
        self.w.close();
        self.written
    }
}

pub fn bio_write_fasta_pull<P>(path: &String, pull: P) -> Vec<String>
where
    P: Fn() -> (bool, String, Vec<BioFastaEntry>),
{
    let mut w = BioWriter::create(path);
    let mut names = Vec::new();
    bio_drain_pull("writeFastaPull", pull, |batch| {
        let mut text = Vec::new();
        for e in batch {
            names.push(e.header.split_whitespace().next().unwrap_or("").to_string());
            bio_fasta_text(&mut text, e);
        }
        w.write(&text);
    });
    w.close();
    names
}

#[derive(Clone, Debug)]
pub struct BioGtfLine {
    pub seqid: String,
    pub origin: String,
    pub kind: String,
    pub part: BioSegment,
    pub score: Option<f64>,
    pub phase: BioPhase,
    pub attrs: String,
}

pub fn bio_stream_gtf_lines<S: Fn(&Vec<BioGtfLine>) -> R, R>(path: &String, sink: S) {
    let mut r = BioLineReader::open(path);
    let mut out = BioBatcher::new(&sink, BIO_BATCH_TREES);
    let mut raw = Vec::new();
    while r.next(&mut raw) {
        if raw.is_empty() || raw[0] == b'#' {
            continue;
        }
        let line = bio_utf8(&raw);
        let (g, attr_col) = bio_parse_gff_columns(&r, &line);
        let attrs = bio_parse_gtf_attrs(&r, &attr_col);
        out.add(
            BioGtfLine { seqid: g.seqid, origin: g.origin, kind: g.kind, part: g.part, score: g.score, phase: g.phase, attrs: bio_attrs_text(&attrs) },
            1,
        );
    }
    out.flush();
}

pub fn bio_gtf_trees(lines: &Vec<BioGtfLine>) -> Vec<BioGffRoot> {
    let roots = std::cell::RefCell::new(Vec::new());
    let it = lines.iter().enumerate().map(|(i, l)| {
        (
            BioGffLine { seqid: l.seqid.clone(), origin: l.origin.clone(), kind: l.kind.clone(), part: l.part, score: l.score, phase: l.phase },
            bio_attrs_parse(&l.attrs),
            i,
        )
    });
    let at = std::cell::Cell::new(0usize);
    bio_gtf_build(
        it.map(|(g, a, i)| {
            at.set(i);
            (g, a)
        }),
        |m: &str| format!("GTF line {}: {}", at.get() + 1, m),
        |xs: &Vec<BioGffRoot>| roots.borrow_mut().extend(xs.iter().cloned()),
    );
    roots.into_inner()
}

pub fn bio_filter_gtf_lines(names: &Vec<String>, check_tx: bool, lines: &Vec<BioGtfLine>) -> Vec<BioGtfLine> {
    let keep: std::collections::HashSet<&str> = names.iter().map(|s| s.as_str()).collect();
    lines
        .iter()
        .filter(|l| {
            keep.contains(l.seqid.as_str())
                && (!check_tx || bio_attr_get(&bio_attrs_parse(&l.attrs), "transcript_id").map_or(false, |v| v.iter().any(|t| !t.is_empty())))
        })
        .cloned()
        .collect()
}

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default)]
pub struct BioTxGene {
    pub tx: String,
    pub gene: String,
    pub geneName: String,
    pub biotype: String,
}

fn bio_first_attr(attrs: &str, keys: &[&str]) -> String {
    let parsed = bio_attrs_parse(attrs);
    for k in keys {
        if let Some(v) = bio_attr_get(&parsed, k) {
            if let Some(x) = v.first() {
                return x.clone();
            }
        }
    }
    String::new()
}

// Each transcript's gene, gene name and gene biotype, in transcriptModels
// order of the trees' entries. A key is looked up on the gene, then on the
// transcript, then on the transcript's children: a GTF without gene and
// transcript lines states names and biotypes only on its exon lines.
pub fn bio_tx_genes(roots: &Vec<BioGffRoot>) -> Vec<BioTxGene> {
    let entries = bio_gff_entries(roots);
    let mut names: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for e in &entries {
        let s = bio_entry_seqid(e);
        if seen.insert(s.clone()) {
            names.push(s);
        }
    }
    // attribute columns to search, by transcript id: gene, transcript, then
    // the transcript's children
    let mut search: std::collections::HashMap<String, Vec<&str>> = std::collections::HashMap::new();
    for r in roots {
        let g = &r.feature;
        let has_exon = |f: &BioFeature<BioGffNode>| f.children.iter().any(|c| c.node.kind == "exon");
        if has_exon(g) {
            let mut cols = vec![g.node.attrs.as_str()];
            cols.extend(g.children.iter().map(|c| c.node.attrs.as_str()));
            search.insert(g.node.id.clone(), cols);
        }
        for t in &g.children {
            if has_exon(t) {
                let mut cols = vec![g.node.attrs.as_str(), t.node.attrs.as_str()];
                cols.extend(t.children.iter().map(|c| c.node.attrs.as_str()));
                search.insert(t.node.id.clone(), cols);
            }
        }
    }
    bio_transcript_models(&names, &entries)
        .into_iter()
        .map(|m| {
            let cols = search.get(&m.id).cloned().unwrap_or_default();
            let pick = |keys: &[&str]| cols.iter().map(|c| bio_first_attr(c, keys)).find(|x| !x.is_empty()).unwrap_or_default();
            let name = pick(&["gene_name", "Name"]);
            BioTxGene {
                tx: m.id,
                geneName: if name.is_empty() { m.gene.clone() } else { name },
                biotype: pick(&["gene_biotype", "gene_type", "biotype"]),
                gene: m.gene,
            }
        })
        .collect()
}


// Read-pair pulls over FASTQ text
// -------------------------------
//
// A pull is a handle into a process-local table of open readers; `next`
// returns the next batch (empty at the end) and drops the reader then.

struct BioPairReader {
    runs: std::collections::VecDeque<(String, String)>,
    cur: Option<(BioLineReader, BioLineReader, String, String)>,
}

fn bio_pair_readers() -> &'static std::sync::Mutex<std::collections::HashMap<i64, BioPairReader>> {
    static R: std::sync::OnceLock<std::sync::Mutex<std::collections::HashMap<i64, BioPairReader>>> = std::sync::OnceLock::new();
    R.get_or_init(|| std::sync::Mutex::new(std::collections::HashMap::new()))
}

pub fn bio_open_pairs(r1: &Vec<String>, r2: &Vec<String>) -> i64 {
    if r1.len() != r2.len() {
        panic!("openPairs: {} first-mate files but {} second-mate files", r1.len(), r2.len());
    }
    static NEXT: std::sync::atomic::AtomicI64 = std::sync::atomic::AtomicI64::new(1);
    let h = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let runs = r1.iter().cloned().zip(r2.iter().cloned()).collect();
    bio_pair_readers().lock().unwrap_or_else(|p| p.into_inner()).insert(h, BioPairReader { runs, cur: None });
    h
}

fn bio_next_pairs_batch(rd: &mut BioPairReader) -> Vec<(BioFastqEntry, BioFastqEntry)> {
    let mut out = Vec::new();
    let mut weight = 0usize;
    let (mut h, mut s, mut p, mut q) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    while weight < BIO_BATCH_BYTES {
        if rd.cur.is_none() {
            match rd.runs.pop_front() {
                None => break,
                Some((a, b)) => rd.cur = Some((BioLineReader::open(&a), BioLineReader::open(&b), a, b)),
            }
        }
        let (ra, rb, pa, pb) = rd.cur.as_mut().unwrap();
        let a = bio_next_fastq(ra, &mut h, &mut s, &mut p, &mut q);
        let b = bio_next_fastq(rb, &mut h, &mut s, &mut p, &mut q);
        match (a, b) {
            (Some(a), Some(b)) => {
                weight += a.header.len() + b.header.len() + 2 * (a.bases.bases.len() + b.bases.bases.len());
                out.push((a, b));
            }
            (None, None) => rd.cur = None,
            _ => panic!("mate files {} and {} hold different numbers of reads", pa, pb),
        }
    }
    out
}

pub fn bio_next_pairs(h: i64) -> (bool, String, Vec<(BioFastqEntry, BioFastqEntry)>) {
    let mut table = bio_pair_readers().lock().unwrap_or_else(|p| p.into_inner());
    let Some(rd) = table.get_mut(&h) else { return (true, String::new(), Vec::new()) };
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| bio_next_pairs_batch(rd))) {
        Ok(xs) => {
            if xs.is_empty() {
                table.remove(&h);
            }
            (true, String::new(), xs)
        }
        Err(e) => {
            table.remove(&h);
            let msg = e
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "read error".to_string());
            (false, msg, Vec::new())
        }
    }
}

struct BioReadReader {
    runs: std::collections::VecDeque<String>,
    cur: Option<(BioLineReader, String)>,
}

fn bio_read_readers() -> &'static std::sync::Mutex<std::collections::HashMap<i64, BioReadReader>> {
    static R: std::sync::OnceLock<std::sync::Mutex<std::collections::HashMap<i64, BioReadReader>>> = std::sync::OnceLock::new();
    R.get_or_init(|| std::sync::Mutex::new(std::collections::HashMap::new()))
}

pub fn bio_open_reads(runs: &Vec<String>) -> i64 {
    static NEXT: std::sync::atomic::AtomicI64 = std::sync::atomic::AtomicI64::new(1);
    let h = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    bio_read_readers().lock().unwrap_or_else(|p| p.into_inner()).insert(h, BioReadReader { runs: runs.iter().cloned().collect(), cur: None });
    h
}

fn bio_next_reads_batch(rd: &mut BioReadReader) -> Vec<BioFastqEntry> {
    let mut out = Vec::new();
    let mut weight = 0usize;
    let (mut h, mut s, mut p, mut q) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    while weight < BIO_BATCH_BYTES {
        if rd.cur.is_none() {
            match rd.runs.pop_front() {
                None => break,
                Some(a) => rd.cur = Some((BioLineReader::open(&a), a)),
            }
        }
        let (r, _) = rd.cur.as_mut().unwrap();
        match bio_next_fastq(r, &mut h, &mut s, &mut p, &mut q) {
            Some(a) => {
                weight += a.header.len() + 2 * a.bases.bases.len();
                out.push(a);
            }
            None => rd.cur = None,
        }
    }
    out
}

pub fn bio_next_reads(h: i64) -> (bool, String, Vec<BioFastqEntry>) {
    let mut table = bio_read_readers().lock().unwrap_or_else(|p| p.into_inner());
    let Some(rd) = table.get_mut(&h) else { return (true, String::new(), Vec::new()) };
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| bio_next_reads_batch(rd))) {
        Ok(xs) => {
            if xs.is_empty() {
                table.remove(&h);
            }
            (true, String::new(), xs)
        }
        Err(e) => {
            table.remove(&h);
            let msg = e
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "read error".to_string());
            (false, msg, Vec::new())
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct BioRefModels {
    pub maps: Vec<BioMapping<BioTranscript, BioGenome>>,
}

pub fn bio_ref_models(names: &Vec<String>, models: &Vec<BioTxModel>) -> BioRefModels {
    let by_id: std::collections::HashMap<&str, &BioTxModel> = models.iter().map(|m| (m.id.as_str(), m)).collect();
    let none = || BioMapping::<BioTranscript, BioGenome> { src: -1, dst: -1, blocks: Vec::new() };
    BioRefModels { maps: names.iter().map(|n| by_id.get(n.as_str()).map(|m| m.mapping.clone()).unwrap_or_else(none)).collect() }
}

pub fn bio_project_point(m: &BioMapping<BioTranscript, BioGenome>, x: i64) -> Option<(i64, bool)> {
    if x < 0 {
        return None;
    }
    let x = x as u64;
    let i = m.blocks.partition_point(|b| b.srcStart + b.len <= x);
    let b = m.blocks.get(i)?;
    if x < b.srcStart {
        return None;
    }
    let off = x - b.srcStart;
    Some(match b.orient {
        BioOrient::Fwd => ((b.dstStart + off) as i64, true),
        BioOrient::Rev => ((b.dstStart + b.len - 1 - off) as i64, false),
    })
}
