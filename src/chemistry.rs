use crate::error::{PrimerError, PrimerResult};
use crate::grammar::Grammar;
use clap::ValueEnum;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Chemistry {
    /// 10x Genomics Single Cell 3' v2
    TenxV2,

    /// 10x Genomics Single Cell 3' v3
    TenxV3,

    /// 10x Genomics Single Cell 3' v4
    TenxV4,

    /// BD Rhapsody v1 (9bp+9bp+9bp cell blocks, 8bp UMI)
    BdV1,

    /// BD Rhapsody v2 96-well barcode set
    BdV2_96,

    /// BD Rhapsody v2 384-well barcode set
    BdV2_384,
}

impl Chemistry {
    pub fn parse(raw: &str) -> PrimerResult<Self> {
        match raw {
            "tenx-v2" => Ok(Self::TenxV2),
            "tenx-v3" => Ok(Self::TenxV3),
            "tenx-v4" => Ok(Self::TenxV4),
            "bd-v1" => Ok(Self::BdV1),
            "bd-v2-96" => Ok(Self::BdV2_96),
            "bd-v2-384" => Ok(Self::BdV2_384),
            other => Err(PrimerError::UnknownChemistry(other.to_string())),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::TenxV2 => "tenx-v2",
            Self::TenxV3 => "tenx-v3",
            Self::TenxV4 => "tenx-v4",
            Self::BdV1 => "bd-v1",
            Self::BdV2_96 => "bd-v2-96",
            Self::BdV2_384 => "bd-v2-384",
        }
    }

    pub fn grammar(self) -> PrimerResult<Grammar> {
        match self {
            Self::TenxV2 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:10+POLYT:min=10+INSERT",
            ),
            Self::TenxV3 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10+INSERT",
            ),
            Self::TenxV4 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10+INSERT",
            ),
            Self::BdV1 => Grammar::parse(self.name(), "BD_CELL:v1+POLYT:min=10+INSERT"),
            Self::BdV2_96 => Grammar::parse(self.name(), "SEARCH:0..4+BD_CELL:v2.96+POLYT:min=10+INSERT"),
            Self::BdV2_384 => Grammar::parse(self.name(), "SEARCH:0..4+BD_CELL:v2.384+POLYT:min=10+INSERT"),
        }
    }
}
