use std::collections::HashMap;

use crate::error::{PrimerError, PrimerResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BdCellVersion {
    V1,
    V2_96,
    V2_384,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RhapsodyCellCall {
    pub version: BdCellVersion,
    pub cell_id: u64,
    pub cell_seq: Vec<u8>,
    pub cell_qual: Vec<u8>,
    pub umi_seq: Vec<u8>,
    pub umi_qual: Vec<u8>,
    pub shift: usize,
    pub consumed: usize,
    pub c1: (usize, usize),
    pub c2: (usize, usize),
    pub c3: (usize, usize),
    pub umi: (usize, usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RhapsodyWhitelist {
    version: BdCellVersion,
    block_size: u64,
    map: HashMap<Vec<u8>, u64>,
}

impl BdCellVersion {
    pub fn parse(raw: &str) -> PrimerResult<Self> {
        match raw {
            "v1" => Ok(Self::V1),
            "v2.96" => Ok(Self::V2_96),
            "v2.384" => Ok(Self::V2_384),
            other => Err(PrimerError::rhapsody(format!("unknown BD cell version '{other}'"))),
        }
    }

    pub fn block_size(self) -> u64 {
        match self {
            Self::V1 => 384,
            Self::V2_96 => 96,
            Self::V2_384 => 384,
        }
    }

    pub fn umi_len(self) -> usize {
        match self {
            Self::V1 => 8,
            Self::V2_96 | Self::V2_384 => 6,
        }
    }

    pub fn unshifted_consumed_len(self) -> usize {
        match self {
            Self::V1 => 60,
            Self::V2_96 | Self::V2_384 => 42,
        }
    }
}

impl RhapsodyCellCall {
    pub fn empty(version: BdCellVersion) -> Self {
        Self {
            version,
            cell_id: 0,
            cell_seq: Vec::new(),
            cell_qual: Vec::new(),
            umi_seq: Vec::new(),
            umi_qual: Vec::new(),
            shift: 0,
            consumed: 0,
            c1: (0, 0),
            c2: (0, 0),
            c3: (0, 0),
            umi: (0, 0),
        }
    }
}

impl RhapsodyWhitelist {
    pub fn new(version: BdCellVersion, entries: Vec<(Vec<u8>, u64)>) -> Self {
        Self {
            version,
            block_size: version.block_size(),
            map: entries.into_iter().collect(),
        }
    }

    pub fn toy_v2_384() -> Self {
        Self::new(
            BdCellVersion::V2_384,
            vec![
                (b"ACGTACGTA".to_vec(), 7),
                (b"TGCATGCAT".to_vec(), 11),
                (b"GATTACAGA".to_vec(), 19),
                (b"CCCCGGGGA".to_vec(), 23),
                (b"AAAACCCCG".to_vec(), 29),
                (b"TTTTGGGGA".to_vec(), 31),
            ],
        )
    }

    pub fn toy_v2_96() -> Self {
        Self::new(
            BdCellVersion::V2_96,
            vec![
                (b"ACGTACGTA".to_vec(), 7),
                (b"TGCATGCAT".to_vec(), 11),
                (b"GATTACAGA".to_vec(), 19),
                (b"CCCCGGGGA".to_vec(), 23),
            ],
        )
    }

    pub fn toy_v1() -> Self {
        Self::new(
            BdCellVersion::V1,
            vec![
                (b"ACGTACGTA".to_vec(), 7),
                (b"TGCATGCAT".to_vec(), 11),
                (b"GATTACAGA".to_vec(), 19),
            ],
        )
    }

    pub fn version(&self) -> BdCellVersion {
        self.version
    }

    pub fn call(&self, seq: &[u8], qual: &[u8], offset: usize, shift_start: usize, shift_end: usize) -> Option<RhapsodyCellCall> {
        for shift in shift_start..=shift_end {
            if let Some(call) = self.call_exact_shift(seq, qual, offset, shift) {
                return Some(call);
            }
        }
        None
    }

    pub fn call_exact_shift(&self, seq: &[u8], qual: &[u8], offset: usize, shift: usize) -> Option<RhapsodyCellCall> {
        let base = offset.checked_add(shift)?;
        let (c1, c2, c3, umi, consumed) = self.coords(base)?;
        if seq.len() < umi.1 || qual.len() < umi.1 {
            return None;
        }
        let c1_idx = self.index_of(&seq[c1.0..c1.1])?;
        let c2_idx = self.index_of(&seq[c2.0..c2.1])?;
        let c3_idx = self.index_of(&seq[c3.0..c3.1])?;
        let cell_id = c1_idx * self.block_size * self.block_size + c2_idx * self.block_size + c3_idx + 1;
        let mut cell_seq = Vec::with_capacity(27);
        let mut cell_qual = Vec::with_capacity(27);
        self.extend_part(&mut cell_seq, &mut cell_qual, seq, qual, c1);
        self.extend_part(&mut cell_seq, &mut cell_qual, seq, qual, c2);
        self.extend_part(&mut cell_seq, &mut cell_qual, seq, qual, c3);
        Some(RhapsodyCellCall {
            version: self.version,
            cell_id,
            cell_seq,
            cell_qual,
            umi_seq: seq[umi.0..umi.1].to_vec(),
            umi_qual: qual[umi.0..umi.1].to_vec(),
            shift,
            consumed,
            c1,
            c2,
            c3,
            umi,
        })
    }

    pub fn expected_id(&self, c1: u64, c2: u64, c3: u64) -> u64 {
        c1 * self.block_size * self.block_size + c2 * self.block_size + c3 + 1
    }

    pub fn index_of(&self, seq: &[u8]) -> Option<u64> {
        self.map.get(seq).copied()
    }

    pub fn coords(&self, base: usize) -> Option<((usize, usize), (usize, usize), (usize, usize), (usize, usize), usize)> {
        match self.version {
            BdCellVersion::V1 => Some(((base, base + 9), (base + 21, base + 30), (base + 43, base + 52), (base + 52, base + 60), base + 60)),
            BdCellVersion::V2_96 | BdCellVersion::V2_384 => Some(((base, base + 9), (base + 13, base + 22), (base + 26, base + 35), (base + 36, base + 42), base + 42)),
        }
    }

    pub fn extend_part(&self, cell_seq: &mut Vec<u8>, cell_qual: &mut Vec<u8>, seq: &[u8], qual: &[u8], range: (usize, usize)) {
        cell_seq.extend_from_slice(&seq[range.0..range.1]);
        cell_qual.extend_from_slice(&qual[range.0..range.1]);
    }
}
