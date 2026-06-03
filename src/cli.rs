use clap::Args;

use crate::{Chemistry, PrimerDetector, RhapsodyWhitelist, Grammar};

#[derive(Debug, Clone, Args)]
pub struct PrimerCli {
    /// Preset single-cell chemistry.
    ///
    /// Ignored if --primer-structure is supplied.
    #[arg(long, value_enum, default_value_t = Chemistry::default())]
    pub chemistry: Chemistry,

    /// Custom primer/read structure grammar.
    ///
    /// Overrides --chemistry.
    #[arg(long)]
    pub primer_structure: Option<String>,

    /// Also search the reverse-complement orientation.
    #[arg(long, default_value_t = true)]
    pub detect_reverse_complement: bool,
}

impl PrimerCli {
    pub fn detector(&self) -> Result<PrimerDetector, String> {
        let detector = if let Some(structure) = self.primer_structure.as_deref() {

            PrimerDetector::from_grammar_with_rhapsody( Grammar::parse( "custom", structure)?, RhapsodyWhitelist::bd_v2_384())
        } else {
            PrimerDetector::from_chemistry(self.chemistry)?
        };

        //detector.set_detect_reverse_complement(self.detect_reverse_complement);

        Ok(detector)
    }
}