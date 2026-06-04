use crate::chemistry::Chemistry;
use crate::error::{PrimerError, PrimerResult};
use crate::grammar::{Grammar, GrammarOp};
use crate::model::{Orientation, PrimerMatch, PrimerAttempt, PrimerSegmentAttempt};
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
            Chemistry::BdV1 => Some(RhapsodyWhitelist::bd_v1()),
            Chemistry::BdV2_96 => Some(RhapsodyWhitelist::bd_v2_96()),
            Chemistry::BdV2_384 => Some(RhapsodyWhitelist::bd_v2_384()),
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

    pub fn detect(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Vec<PrimerMatch>> {
        self.detect_all(seq, qual)
    }

    pub fn detect_first(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Option<PrimerMatch>> {
        if seq.len() != qual.len() {
            return Err(PrimerError::invalid_coordinates(
                "sequence and quality have different lengths",
            ));
        }

        for offset in 0..seq.len() {
            if let Some(mut hit) = self.detect_one_orientation(
                &seq[offset..],
                &qual[offset..],
                Orientation::Forward,
            )? {
                hit.shift_by(offset);
                hit.insert_end = seq.len();
                return Ok(Some(hit));
            }
        }

        if self.detect_reverse_complement {
            let rc_seq = Self::reverse_complement(seq);
            let rc_qual = Self::reverse(qual);

            for offset in 0..rc_seq.len() {
                if let Some(mut hit) = self.detect_one_orientation(
                    &rc_seq[offset..],
                    &rc_qual[offset..],
                    Orientation::ReverseComplement,
                )? {
                    hit.shift_by(offset);
                    hit.insert_end = rc_seq.len();
                    hit.remap_reverse_coordinates(seq.len());
                    return Ok(Some(hit));
                }
            }
        }

        Ok(None)
    }

    pub fn detect_all(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Vec<PrimerMatch>> {
        if seq.len() != qual.len() {
            return Err(PrimerError::invalid_coordinates("sequence and quality have different lengths"));
        }

        let mut hits = Vec::new();
        let mut offset = 0usize;
        while offset < seq.len() {
            let Some(mut hit) = self.try_from_start(seq, qual, offset, Orientation::Forward)? else {
                offset += 1;
                continue;
            };
            hit.shift_by(offset);
            let next_offset = hit.primer_end.max(offset + 1);
            hits.push(hit);
            offset = next_offset;
        }

        self.finish_insert_ends(&mut hits, seq.len());
        Ok(hits)
    }

    pub fn explain_all(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Vec<PrimerAttempt>> {
        if seq.len() != qual.len() {
            return Err(PrimerError::invalid_coordinates(
                "sequence and quality have different lengths",
            ));
        }

        let mut attempts = Vec::new();

        for offset in 0..seq.len() {
            //let sub_seq = &seq[offset..];
            //let sub_qual = &qual[offset..];
            match self.try_from_start(seq, qual, offset, Orientation::Forward)? {
                Some(hit) => {
                    attempts.push(PrimerAttempt {
                        offset,
                        orientation: Orientation::Forward,
                        ok: true,
                        reason: "matched".to_string(),
                        segments: hit
                            .segments
                            .iter()
                            .flat_map(|segment| {
                                segment.ranges.iter().map(|range| PrimerSegmentAttempt {
                                    name: segment.name.clone(),
                                    range: range.clone(),
                                    dna: String::from_utf8_lossy(&seq[range.start..range.end])
                                        .to_string(),
                                    ok: true,
                                    reason: "matched".to_string(),
                                })
                            })
                            .collect(),
                    });

                    break;
                }

                None => {
                    attempts.push(PrimerAttempt {
                        offset,
                        orientation: Orientation::Forward,
                        ok: false,
                        reason: "no complete primer match".to_string(),
                        segments: Vec::new(),
                    });
                }
            }
        }

        Ok(attempts)
    }

    pub fn grammar(&self) -> &Grammar {
        &self.grammar
    }

    pub fn detect_one_orientation(&self, seq: &[u8], qual: &[u8], orientation: Orientation) -> PrimerResult<Option<PrimerMatch>> {
        self.try_from_start(seq, qual, 0, orientation)
    }

    pub fn try_from_start(&self, seq: &[u8], qual: &[u8], start: usize, orientation: Orientation) -> PrimerResult<Option<PrimerMatch>> {
        let mut primer_match = PrimerMatch::new(self.grammar.name.clone(), orientation);
        primer_match.primer_start = start;
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
                    primer_match.add_segment("CELL", pos..pos + *len);
                    pos += *len;
                }
                GrammarOp::Umi { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) { return Ok(None); }
                    primer_match.add_segment("UMI", pos..pos + *len);
                    pos += *len;
                }
                GrammarOp::PolyT { min } => {
                    let count = Self::count_base(seq, pos, b'T');
                    if count < *min { return Ok(None); }
                    pos += count;
                }
                GrammarOp::Insert => {
                    primer_match.insert_start = pos;
                    primer_match.insert_end = seq.len();
                    primer_match.primer_end = pos;
                    return Ok(Some(primer_match));
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
                    primer_match.bd_cell_id = Some(call.cell_id);
                    primer_match.add_segment_ranges(
                        "BD_CELL",
                        vec![call.c1.0..call.c1.1, call.c2.0..call.c2.1, call.c3.0..call.c3.1],
                    );
                    primer_match.add_segment("UMI", call.umi.0..call.umi.1);
                    pos = call.consumed;
                    search = (0, 0);
                }
                GrammarOp::Tag { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) { return Ok(None); }
                    primer_match.add_segment("TAG", pos..pos + *len);
                    pos += *len;
                }
                GrammarOp::Feature { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) { return Ok(None); }
                    primer_match.add_segment("FEATURE", pos..pos + *len);
                    pos += *len;
                }
            }
        }
        primer_match.primer_end = pos;
        primer_match.insert_start = pos;
        primer_match.insert_end = seq.len();
        Ok(Some(primer_match))
    }

    pub fn finish_insert_ends(&self, hits: &mut [PrimerMatch], seq_len: usize) {
        if hits.is_empty() {
            return;
        }
        for idx in 0..hits.len() {
            hits[idx].insert_end = if idx + 1 < hits.len() {
                hits[idx + 1].primer_start
            } else {
                seq_len
            };
        }
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
            return Err(PrimerError::invalid_coordinates("hamming requires equal lengths"));
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
