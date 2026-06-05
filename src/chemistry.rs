use crate::error::{PrimerError, PrimerResult};
use crate::grammar::Grammar;

use clap::ValueEnum;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum)]
pub enum Chemistry {
    /// 10x Genomics Chromium Single Cell 3' v2.
    ///
    /// Layout: adapter + 16 bp cell barcode + 10 bp UMI + polyT + insert.
    TenxV2,

    /// 10x Genomics Chromium Single Cell 3' v3.
    ///
    /// Layout: adapter + 16 bp cell barcode + 12 bp UMI + polyT + insert.
    #[default]
    TenxV3,

    /// 10x Genomics Chromium Single Cell 3' v4 / GEM-X style preset.
    ///
    /// Layout: adapter + 16 bp cell barcode + 12 bp UMI + polyT + insert.
    TenxV4,

    /// BD Rhapsody v1 / older layout.
    ///
    /// Uses three 9 bp barcode blocks, BD whitelist lookup, long linker gaps,
    /// and an 8 bp UMI.
    BdV1,

    /// BD Rhapsody v2 96-cell combinatorial barcode layout.
    ///
    /// Uses three 9 bp barcode blocks, BD whitelist lookup, shifts 0..4,
    /// and the 96-block barcode ID formula.
    BdV2_96,

    /// BD Rhapsody v2 384-cell combinatorial barcode layout.
    ///
    /// Uses three 9 bp barcode blocks, BD whitelist lookup, shifts 0..4,
    /// and the 384-block barcode ID formula.
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
            other => Err(PrimerError::unknown_chemistry(other)),
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
        let bd_fixed= "ATAGGAAACTCATGGT";

        match self {
            Self::TenxV2 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:10+POLYT:min=0+INSERT",
            ),
            Self::TenxV3 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=0+INSERT",
            ),
            Self::TenxV4 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=0+INSERT",
            ),
            Self::BdV1 => Grammar::parse(self.name(), "BD_CELL:v1+POLYT:min=0+INSERT"),
            /*
            Self::BdV2_96 => Grammar::parse(self.name(), &format!("SEARCH:0..4+BD_CELL:v2.96+POLYT:min=0+INSERT")),
            Self::BdV2_384 => Grammar::parse(self.name(), &format!("SEARCH:0..4+BD_CELL:v2.384+POLYT:min=0+INSERT")),
            */
            Self::BdV2_96 => Grammar::parse(self.name(), &format!("SEARCH:0..16+BD_CELL:v2.96+POLYT:min=0+INSERT")),
            Self::BdV2_384 => Grammar::parse(self.name(), &format!("SEARCH:0..16+BD_CELL:v2.384+POLYT:min=0+INSERT")),
            
        }
    }
}
