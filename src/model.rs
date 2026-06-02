use crate::error::{PrimerError, PrimerResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    ReverseComplement,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrimerRange {
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimerPart {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimerSegment {
    pub name: String,
    pub ranges: Vec<PrimerRange>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimerMatch {
    pub chemistry_name: String,
    pub orientation: Orientation,
    pub primer_start: usize,
    pub primer_end: usize,
    pub insert_start: usize,
    pub insert_end: usize,
    pub cell_ranges: Vec<PrimerRange>,
    pub umi_range: Option<PrimerRange>,
    pub bd_cell_id: Option<u64>,
    pub segments: Vec<PrimerSegment>,
}

impl PrimerRange {
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    pub fn len(self) -> usize {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(self) -> bool {
        self.start >= self.end
    }

    pub fn shift(&mut self, offset: usize) {
        self.start += offset;
        self.end += offset;
    }

    pub fn remap_reverse(&mut self, read_len: usize) {
        let old_start = self.start;
        let old_end = self.end;
        self.start = read_len - old_end;
        self.end = read_len - old_start;
    }

    pub fn validate(self, seq: &[u8], qual: &[u8]) -> PrimerResult<()> {
        if seq.len() != qual.len() {
            return Err(PrimerError::InvalidCoordinates(
                "sequence and quality have different lengths".to_string(),
            ));
        }
        if self.start > self.end || self.end > seq.len() {
            return Err(PrimerError::InvalidCoordinates(format!(
                "invalid primer range {}..{} for read length {}",
                self.start,
                self.end,
                seq.len()
            )));
        }
        Ok(())
    }
}

impl PrimerPart {
    pub fn from_reverse_complement_ranges(seq: &[u8], qual: &[u8], ranges: &[PrimerRange]) -> PrimerResult<Self> {
        let mut part = Self::from_ranges(seq, qual, ranges)?;
        part.seq = part.seq.iter().rev().map(|base| Self::complement(*base)).collect();
        part.qual.reverse();
        Ok(part)
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

    pub fn from_ranges(seq: &[u8], qual: &[u8], ranges: &[PrimerRange]) -> PrimerResult<Self> {
        let total_len = ranges.iter().map(|range| range.len()).sum();
        let mut part = Self {
            seq: Vec::with_capacity(total_len),
            qual: Vec::with_capacity(total_len),
        };
        for range in ranges {
            range.validate(seq, qual)?;
            part.seq.extend_from_slice(&seq[range.start..range.end]);
            part.qual.extend_from_slice(&qual[range.start..range.end]);
        }
        Ok(part)
    }

    pub fn from_range(seq: &[u8], qual: &[u8], range: PrimerRange) -> PrimerResult<Self> {
        Self::from_ranges(seq, qual, &[range])
    }

    pub fn seq_string(&self) -> PrimerResult<String> {
        String::from_utf8(self.seq.clone()).map_err(|err| {
            PrimerError::InvalidCoordinates(format!("primer sequence is not UTF-8/ASCII: {err}"))
        })
    }

    pub fn qual_string(&self) -> PrimerResult<String> {
        String::from_utf8(self.qual.clone()).map_err(|err| {
            PrimerError::InvalidCoordinates(format!("primer quality is not UTF-8/ASCII: {err}"))
        })
    }
}

impl PrimerSegment {
    pub fn new(name: &str, ranges: Vec<PrimerRange>) -> Self {
        Self {
            name: name.to_string(),
            ranges,
        }
    }

    pub fn shift(&mut self, offset: usize) {
        for range in &mut self.ranges {
            range.shift(offset);
        }
    }

    pub fn remap_reverse(&mut self, read_len: usize) {
        for range in &mut self.ranges {
            range.remap_reverse(read_len);
        }
        self.ranges.sort_by_key(|range| range.start);
    }
}

impl PrimerMatch {
    pub fn new(chemistry_name: String, orientation: Orientation) -> Self {
        Self {
            chemistry_name,
            orientation,
            primer_start: 0,
            primer_end: 0,
            insert_start: 0,
            insert_end: 0,
            cell_ranges: Vec::new(),
            umi_range: None,
            bd_cell_id: None,
            segments: Vec::new(),
        }
    }

    pub fn add_cell_range(&mut self, start: usize, end: usize) {
        let range = PrimerRange::new(start, end);
        self.cell_ranges.push(range);
        self.segments.push(PrimerSegment::new("CELL", vec![range]));
    }

    pub fn set_cell_ranges(&mut self, ranges: Vec<PrimerRange>) {
        self.cell_ranges = ranges.clone();
        self.segments.push(PrimerSegment::new("BD_CELL", ranges));
    }

    pub fn set_umi_range(&mut self, start: usize, end: usize) {
        let range = PrimerRange::new(start, end);
        self.umi_range = Some(range);
        self.segments.push(PrimerSegment::new("UMI", vec![range]));
    }

    pub fn add_named_range(&mut self, name: &str, start: usize, end: usize) {
        self.segments.push(PrimerSegment::new(name, vec![PrimerRange::new(start, end)]));
    }

    pub fn get_cell(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<PrimerPart> {
        if self.cell_ranges.is_empty() {
            return Err(PrimerError::InvalidCoordinates("primer match does not contain a CELL segment".to_string()));
        }
        self.extract_part(seq, qual, &self.cell_ranges)
    }

    pub fn get_umi(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<PrimerPart> {
        let Some(range) = self.umi_range else {
            return Err(PrimerError::InvalidCoordinates("primer match does not contain a UMI segment".to_string()));
        };
        self.extract_part(seq, qual, &[range])
    }

    pub fn get_insert(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<PrimerPart> {
        self.extract_part(seq, qual, &[PrimerRange::new(self.insert_start, self.insert_end)])
    }

    pub fn extract_part(&self, seq: &[u8], qual: &[u8], ranges: &[PrimerRange]) -> PrimerResult<PrimerPart> {
        match self.orientation {
            Orientation::Forward => PrimerPart::from_ranges(seq, qual, ranges),
            Orientation::ReverseComplement => PrimerPart::from_reverse_complement_ranges(seq, qual, ranges),
        }
    }

    pub fn shift(&mut self, offset: usize) {
        self.primer_start += offset;
        self.primer_end += offset;
        self.insert_start += offset;
        self.insert_end += offset;
        for range in &mut self.cell_ranges {
            range.shift(offset);
        }
        if let Some(range) = &mut self.umi_range {
            range.shift(offset);
        }
        for segment in &mut self.segments {
            segment.shift(offset);
        }
    }

    pub fn remap_reverse(&mut self, read_len: usize) {
        let old_start = self.primer_start;
        let old_end = self.primer_end;
        let old_insert_start = self.insert_start;
        let old_insert_end = self.insert_end;
        self.primer_start = read_len - old_end;
        self.primer_end = read_len - old_start;
        self.insert_start = read_len - old_insert_end;
        self.insert_end = read_len - old_insert_start;
        for range in &mut self.cell_ranges {
            range.remap_reverse(read_len);
        }
        self.cell_ranges.sort_by_key(|range| range.start);
        if let Some(range) = &mut self.umi_range {
            range.remap_reverse(read_len);
        }
        for segment in &mut self.segments {
            segment.remap_reverse(read_len);
        }
    }
}
