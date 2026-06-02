use crate::error::{PrimerError, PrimerResult};
use crate::rhapsody::BdCellVersion;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Grammar {
    pub name: String,
    pub ops: Vec<GrammarOp>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GrammarOp {
    Fixed { seq: Vec<u8>, mismatches: usize },
    Cell { len: usize },
    Umi { len: usize },
    PolyT { min: usize },
    Insert,
    Skip { len: usize },
    Search { start: usize, end: usize },
    BdCell { version: BdCellVersion },
    Tag { len: usize },
    Feature { len: usize },
}

impl Grammar {
    pub fn new(name: impl Into<String>, ops: Vec<GrammarOp>) -> Self {
        Self { name: name.into(), ops }
    }

    pub fn parse(name: impl Into<String>, structure: &str) -> PrimerResult<Self> {
        let mut ops = Vec::new();
        for raw in structure.split('+') {
            let token = raw.trim();
            if token.is_empty() {
                continue;
            }
            ops.push(GrammarOp::parse_token(token)?);
        }
        if ops.is_empty() {
            return Err(PrimerError::InvalidGrammar("empty grammar".to_string()));
        }
        Ok(Self::new(name, ops))
    }
}

impl GrammarOp {
    pub fn parse_token(token: &str) -> PrimerResult<Self> {
        if let Some(rest) = token.strip_prefix("FIXED:") {
            return Self::parse_fixed(rest);
        }
        if let Some(rest) = token.strip_prefix("CELL:") {
            return Self::parse_len(rest).map(|len| Self::Cell { len });
        }
        if let Some(rest) = token.strip_prefix("UMI:") {
            return Self::parse_len(rest).map(|len| Self::Umi { len });
        }
        if let Some(rest) = token.strip_prefix("POLYT:") {
            return Self::parse_poly_t(rest);
        }
        if token == "INSERT" {
            return Ok(Self::Insert);
        }
        if let Some(rest) = token.strip_prefix("SKIP:") {
            return Self::parse_len(rest).map(|len| Self::Skip { len });
        }
        if let Some(rest) = token.strip_prefix("SEARCH:") {
            return Self::parse_search(rest);
        }
        if let Some(rest) = token.strip_prefix("BD_CELL:") {
            return BdCellVersion::parse(rest).map(|version| Self::BdCell { version });
        }
        if let Some(rest) = token.strip_prefix("TAG:") {
            return Self::parse_len(rest).map(|len| Self::Tag { len });
        }
        if let Some(rest) = token.strip_prefix("FEATURE:") {
            return Self::parse_len(rest).map(|len| Self::Feature { len });
        }
        Err(PrimerError::InvalidGrammar(format!("unknown token '{token}'")))
    }

    pub fn parse_fixed(rest: &str) -> PrimerResult<Self> {
        let mut fields = rest.split(':');
        let seq = fields.next().unwrap_or_default().as_bytes().to_vec();
        if seq.is_empty() {
            return Err(PrimerError::InvalidGrammar("FIXED sequence is empty".to_string()));
        }
        let mut mismatches = 0usize;
        for field in fields {
            if let Some(mm) = field.strip_prefix("mm=") {
                mismatches = Self::parse_len(mm)?;
            } else {
                return Err(PrimerError::InvalidGrammar(format!("unknown FIXED option '{field}'")));
            }
        }
        Ok(Self::Fixed { seq, mismatches })
    }

    pub fn parse_poly_t(rest: &str) -> PrimerResult<Self> {
        if let Some(min) = rest.strip_prefix("min=") {
            return Self::parse_len(min).map(|min| Self::PolyT { min });
        }
        Err(PrimerError::InvalidGrammar(format!("bad POLYT token '{rest}'")))
    }

    pub fn parse_search(rest: &str) -> PrimerResult<Self> {
        let Some((start, end)) = rest.split_once("..") else {
            return Err(PrimerError::InvalidGrammar(format!("bad SEARCH range '{rest}'")));
        };
        let start = Self::parse_len(start)?;
        let end = Self::parse_len(end)?;
        if end < start {
            return Err(PrimerError::InvalidGrammar(format!("bad SEARCH range '{rest}'")));
        }
        Ok(Self::Search { start, end })
    }

    pub fn parse_len(raw: &str) -> PrimerResult<usize> {
        raw.parse::<usize>()
            .map_err(|_| PrimerError::InvalidGrammar(format!("bad integer '{raw}'")))
    }
}
