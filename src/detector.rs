use crate::chemistry::Chemistry;
use crate::error::{PrimerError, PrimerResult};
use crate::grammar::{Grammar, GrammarOp};
use crate::model::{Orientation, PrimerMatch};
use crate::rhapsody::{BdCellVersion, RhapsodyWhitelist};

#[derive(Debug, Clone)]
pub struct PrimerDetector {
    grammar: Grammar,
    rhapsody: Option<RhapsodyWhitelist>,
    detect_reverse_complement: bool,
}

impl PrimerDetector {
    pub fn from_chemistry(chemistry: Chemistry) -> PrimerResult<Self> {
        let grammar = chemistry.grammar()?;
        let rhapsody = match chemistry {
            Chemistry::BdV1 => Some(RhapsodyWhitelist::toy_v1()),
            Chemistry::BdV2_96 => Some(RhapsodyWhitelist::toy_v2_96()),
            Chemistry::BdV2_384 => Some(RhapsodyWhitelist::toy_v2_384()),
            _ => None,
        };
        Ok(Self { grammar, rhapsody, detect_reverse_complement: true })
    }

    pub fn from_grammar(grammar: Grammar) -> Self {
        Self { grammar, rhapsody: None, detect_reverse_complement: true }
    }

    pub fn from_grammar_with_rhapsody(grammar: Grammar, rhapsody: RhapsodyWhitelist) -> Self {
        Self { grammar, rhapsody: Some(rhapsody), detect_reverse_complement: true }
    }

    pub fn detect(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Option<PrimerMatch>> {
        if seq.len() != qual.len() {
            return Err(PrimerError::InvalidCoordinates("sequence and quality have different lengths".to_string()));
        }
        if let Some(hit) = self.detect_one_orientation(seq, qual, Orientation::Forward)? {
            return Ok(Some(hit));
        }
        if self.detect_reverse_complement {
            let rc_seq = Self::reverse_complement(seq);
            let rc_qual = Self::reverse(qual);
            if let Some(mut hit) = self.detect_one_orientation(&rc_seq, &rc_qual, Orientation::ReverseComplement)? {
                self.remap_reverse_hit(&mut hit, seq.len());
                return Ok(Some(hit));
            }
        }
        Ok(None)
    }

    pub fn detect_all(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Vec<PrimerMatch>> {
        if seq.len() != qual.len() {
            return Err(PrimerError::InvalidCoordinates("sequence and quality have different lengths".to_string()));
        }
        let mut hits = Vec::new();
        let mut offset = 0usize;
        while offset < seq.len() {
            let sub_seq = &seq[offset..];
            let sub_qual = &qual[offset..];
            let Some(mut hit) = self.detect_one_orientation(sub_seq, sub_qual, Orientation::Forward)? else {
                offset += 1;
                continue;
            };
            hit.shift(offset);
            let next = hit.insert_start.max(hit.primer_end).max(offset + 1);
            hits.push(hit);
            offset = next;
        }
        self.close_insert_ranges(&mut hits, seq.len());
        Ok(hits)
    }

    pub fn grammar(&self) -> &Grammar {
        &self.grammar
    }

    pub fn detect_one_orientation(&self, seq: &[u8], qual: &[u8], orientation: Orientation) -> PrimerResult<Option<PrimerMatch>> {
        self.try_from_start(seq, qual, 0, orientation)
    }

    pub fn try_from_start(&self, seq: &[u8], qual: &[u8], start: usize, orientation: Orientation) -> PrimerResult<Option<PrimerMatch>> {
        let mut hit = PrimerMatch::new(self.grammar.name.clone(), orientation);
        hit.primer_start = start;
        let mut pos = start;
        let mut search = (0usize, 0usize);

        for op in &self.grammar.ops {
            match op {
                GrammarOp::Fixed { seq: fixed, mismatches } => {
                    let chosen = Self::find_fixed(seq, pos, fixed, *mismatches, search)?;
                    let Some(next_pos) = chosen else { return Ok(None); };
                    pos = next_pos + fixed.len();
                    search = (0, 0);
                }
                GrammarOp::Cell { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) { return Ok(None); }
                    hit.add_cell_range(pos, pos + *len);
                    pos += *len;
                }
                GrammarOp::Umi { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) { return Ok(None); }
                    hit.set_umi_range(pos, pos + *len);
                    pos += *len;
                }
                GrammarOp::PolyT { min } => {
                    let count = Self::count_base(seq, pos, b'T');
                    if count < *min { return Ok(None); }
                    pos += count;
                }
                GrammarOp::Insert => {
                    hit.insert_start = pos;
                    hit.insert_end = seq.len();
                    hit.primer_end = pos;
                    return Ok(Some(hit));
                }
                GrammarOp::Skip { len } => {
                    if !Self::has_range(seq, pos, *len) { return Ok(None); }
                    pos += *len;
                }
                GrammarOp::Search { start, end } => {
                    search = (*start, *end);
                }
                GrammarOp::BdCell { version } => {
                    let Some(rhapsody) = &self.rhapsody else { return Ok(None); };
                    if rhapsody.version() != *version { return Ok(None); }
                    let Some(call) = rhapsody.call(seq, qual, pos, search.0, search.1) else { return Ok(None); };
                    hit.set_cell_ranges(call.cell_ranges.clone());
                    hit.set_umi_range(call.umi_range.start, call.umi_range.end);
                    hit.bd_cell_id = Some(call.cell_id);
                    pos = call.consumed;
                    search = (0, 0);
                }
                GrammarOp::Tag { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) { return Ok(None); }
                    hit.add_named_range("TAG", pos, pos + *len);
                    pos += *len;
                }
                GrammarOp::Feature { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) { return Ok(None); }
                    hit.add_named_range("FEATURE", pos, pos + *len);
                    pos += *len;
                }
            }
        }
        hit.primer_end = pos;
        hit.insert_start = pos;
        hit.insert_end = seq.len();
        Ok(Some(hit))
    }

    pub fn remap_reverse_hit(&self, hit: &mut PrimerMatch, len: usize) {
        hit.remap_reverse(len);
    }

    pub fn shift_hit(&self, hit: &mut PrimerMatch, offset: usize) {
        hit.shift(offset);
    }

    pub fn close_insert_ranges(&self, hits: &mut [PrimerMatch], read_len: usize) {
        if hits.is_empty() {
            return;
        }
        for index in 0..hits.len() - 1 {
            hits[index].insert_end = hits[index + 1].primer_start;
        }
        let last = hits.len() - 1;
        hits[last].insert_end = read_len;
    }

    pub fn has_range(seq: &[u8], start: usize, len: usize) -> bool {
        start.checked_add(len).is_some_and(|end| end <= seq.len())
    }

    pub fn count_base(seq: &[u8], start: usize, base: u8) -> usize {
        seq.iter().skip(start).take_while(|b| **b == base).count()
    }

    pub fn find_fixed(seq: &[u8], pos: usize, fixed: &[u8], mismatches: usize, search: (usize, usize)) -> PrimerResult<Option<usize>> {
        let start = pos + search.0;
        let end = pos + search.1;
        for candidate in start..=end {
            if !Self::has_range(seq, candidate, fixed.len()) { continue; }
            if Self::hamming(&seq[candidate..candidate + fixed.len()], fixed)? <= mismatches {
                return Ok(Some(candidate));
            }
        }
        Ok(None)
    }

    pub fn hamming(left: &[u8], right: &[u8]) -> PrimerResult<usize> {
        if left.len() != right.len() {
            return Err(PrimerError::InvalidCoordinates("hamming requires equal lengths".to_string()));
        }
        Ok(left.iter().zip(right.iter()).filter(|(a, b)| a != b).count())
    }

    pub fn reverse(seq: &[u8]) -> Vec<u8> {
        seq.iter().rev().copied().collect()
    }

    pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter().rev().map(|b| Self::complement(*b)).collect()
    }

    pub fn complement(base: u8) -> u8 {
        match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            other => other,
        }
    }

    pub fn version_from_bd_op(version: BdCellVersion) -> BdCellVersion {
        version
    }
}
