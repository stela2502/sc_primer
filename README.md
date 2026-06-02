# sc_primer

`sc_primer` is a reusable Rust crate for detecting and extracting single-cell
primer/read structures from sequencing reads.

It is intended as shared infrastructure for tools such as `bam_tide`, ONT
normalizers, Illumina normalizers, and other single-cell preprocessing binaries.

The crate supports both:

- simple fixed-position chemistries such as 10x Genomics
- shifted/search-based chemistries such as BD Rhapsody
- custom user-supplied primer grammars

A major design goal is that the same detector can be used for both Illumina and
ONT reads, including ONT reads that contain multiple concatenated
primer+insert molecules.

---

# Features

- Chemistry presets:
  
  - `tenx-v2`
  - `tenx-v3`
  - `tenx-v4`
  - `bd-v1`
  - `bd-v2-96`
  - `bd-v2-384`

- Custom primer grammar with `--primer-structure`

- Reusable `clap` CLI block via `PrimerCli`

- Reverse-complement detection

- `detect_all()` for ONT multimers

- Coordinate-backed `PrimerMatch` values

- `get_cell()`, `get_umi()`, and `get_insert()` helpers

- Normal FASTQ sequence/quality handling with `&[u8]`

- BD Rhapsody whitelist/block logic

---

# Installation

Use from crates.io:

```toml
[dependencies]
sc_primer = "0.1"
```

or from a local checkout:

```toml
[dependencies]
sc_primer = { path = "../sc_primer" }
```

If your binary wants to use the reusable CLI block, also enable `clap` in your
binary crate:

```toml
[dependencies]
clap = { version = "4.5", features = ["derive"] }
sc_primer = { path = "../sc_primer" }
```

---

# Quick Rust API Example

```rust
use sc_primer::{Chemistry, PrimerDetector};

fn main() -> Result<(), String> {
    let seq = b"CTACACGACGCTCTTCCGATCTACGTACGTACGTACGTGATCGATCGATTTTTTTTTTTTTTTTGGCCTTAA";
    let qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3)?;
    let hits = detector.detect_all(seq, qual)?;

    for hit in hits {
        let cell = hit.get_cell(seq, qual)?;
        let umi = hit.get_umi(seq, qual)?;
        let insert = hit.get_insert(seq, qual)?;

        println!("cell={}", String::from_utf8_lossy(cell.seq));
        println!("umi={}", String::from_utf8_lossy(umi.seq));
        println!("insert={}", String::from_utf8_lossy(insert.seq));
    }

    Ok(())
}
```

`PrimerMatch` stores coordinates only. It does not copy sequence or quality
data. The original FASTQ `seq` and `qual` remain normal byte strings.

---

# Reusable CLI Integration

`sc_primer` provides a reusable `clap::Args` block:

```rust
use clap::Args;

use crate::{Chemistry, PrimerDetector};

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
        let mut detector = if let Some(structure) = self.primer_structure.as_deref() {
            PrimerDetector::from_structure("custom", structure)?
        } else {
            PrimerDetector::from_chemistry(self.chemistry)?
        };

        detector.set_detect_reverse_complement(self.detect_reverse_complement);

        Ok(detector)
    }
}
```

In a downstream binary, flatten it into your own CLI:

```rust
use clap::Parser;

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Cli {
    #[command(flatten)]
    pub primer: sc_primer::PrimerCli,

    #[arg(long)]
    pub input: String,

    #[arg(long)]
    pub output: String,
}
```

Then create the detector from the parsed CLI object:

```rust
use clap::Parser;

fn main() -> Result<(), String> {
    let cli = Cli::parse();

    let detector = cli.primer.detector()?;

    // ONT-style usage:
    // let hits = detector.detect_all(seq, qual)?;

    // Illumina-style usage:
    // let hit = detector.detect_first(r1_seq, r1_qual)?;

    Ok(())
}
```

This gives every binary the same primer interface:

```bash
my-normalizer \
  --chemistry tenx-v3 \
  --input reads.fastq.gz \
  --output normalized.fastq.gz
```

or:

```bash
my-normalizer \
  --chemistry bd-v2-384 \
  --input reads.fastq.gz \
  --output normalized.fastq.gz
```

or with a completely custom grammar:

```bash
my-normalizer \
  --primer-structure 'FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10+INSERT' \
  --input reads.fastq.gz \
  --output normalized.fastq.gz
```

The rule is simple:

```text
if --primer-structure is supplied, it overrides --chemistry
otherwise --chemistry is used
```

---

# Fancy Chemistry Help Strings

`Chemistry` should derive `clap::ValueEnum`, and the variant documentation is
used by `clap` in generated help text.

Example:

```rust
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
    /// Uses three 9 bp barcode blocks, BD whitelist lookup, longer linker gaps,
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
```

With the flattened `PrimerCli`, a binary automatically exposes:

```bash
my-normalizer --help
```

with options similar to:

```text
--chemistry <CHEMISTRY>
    Preset single-cell chemistry.

    Possible values:
    tenx-v2
    tenx-v3
    tenx-v4
    bd-v1
    bd-v2-96
    bd-v2-384

--primer-structure <PRIMER_STRUCTURE>
    Custom primer/read structure grammar. Overrides --chemistry.

--detect-reverse-complement
    Also search the reverse-complement orientation.
```

---

# Chemistry Presets

## 10x Genomics

10x-style chemistries are adapter anchored.

Conceptually:

```text
FIXED + CELL + UMI + POLYT + INSERT
```

Example preset grammar for 10x v3:

```text
FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10+INSERT
```

## BD Rhapsody

BD Rhapsody is not just a fixed cell barcode plus UMI.

It uses three barcode blocks and whitelist lookup logic.

The v2 layouts also use the historical Rustody-style shift search:

```text
SEARCH:0..4
```

The presets are represented through:

```text
BD_CELL:v1
BD_CELL:v2.96
BD_CELL:v2.384
```

These primitives hide the BD-specific block positions, whitelist lookup, UMI
positioning, and final cell-id calculation.

---

# Custom Grammar

A custom grammar can be supplied with `--primer-structure` or directly through
the Rust API:

```rust
let detector = PrimerDetector::from_structure(
    "custom",
    "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10+INSERT",
)?;
```

Supported primitives:

```text
FIXED:<SEQ>:mm=<N>
CELL:<LEN>
UMI:<LEN>
POLYT:min=<N>
INSERT
SKIP:<LEN>
SEARCH:<START>..<END>
BD_CELL:v1
BD_CELL:v2.96
BD_CELL:v2.384
```

Planned future primitives:

```text
TAG:<LEN>
FEATURE:<LEN>
```

---

# ONT Multimer Detection

ONT reads may contain multiple concatenated primer+insert molecules:

```text
adapter-cell-umi-polyT-insert
adapter-cell-umi-polyT-insert
adapter-cell-umi-polyT-insert
```

Use:

```rust
let hits = detector.detect_all(seq, qual)?;
```

This returns:

```rust
Vec<PrimerMatch>
```

Each `PrimerMatch` can slice its own cell, UMI, and insert:

```rust
for hit in hits {
    let cell = hit.get_cell(seq, qual)?;
    let umi = hit.get_umi(seq, qual)?;
    let insert = hit.get_insert(seq, qual)?;
}
```

This is the preferred API for ONT normalizers.

---

# Illumina Usage

For Illumina read normalization, usually only the first valid primer structure
in R1 is needed:

```rust
let hit = detector.detect_first(r1_seq, r1_qual)?;
let cell = hit.get_cell(r1_seq, r1_qual)?;
let umi = hit.get_umi(r1_seq, r1_qual)?;
```

Then the normalizer can emit artificial R1/R2 FASTQ records and a side-channel
read tag table.

For BD Rhapsody, R1 should be processed for:

- cell id
- UMI
- quality
- trailing/polyT information
- optional sample tags/features

R2 can stay mapper-facing, for example for STAR.

---

# PrimerMatch and PrimerSlice

`PrimerMatch` is coordinate-backed.

It is designed for cheap reuse with large reads:

```rust
let cell = hit.get_cell(seq, qual)?;
```

The returned slice has the shape:

```rust
pub struct PrimerSlice<'a> {
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}
```

This avoids copying sequence and quality data into every match.

---

# Error Handling

The crate intentionally uses simple errors:

```rust
Result<T, String>
```

This keeps downstream crates simple. Users can convert into `anyhow::Error` if
they want:

```rust
let detector = cli.primer.detector().map_err(anyhow::Error::msg)?;
```

---

# Design Philosophy

- Preserve sound Rustody logic
- Keep BD/Rhapsody whitelist logic available
- Avoid locking everything into fixed 10x-style layouts
- Make presets convenient but not mandatory
- Keep custom grammar as an escape hatch
- Support ONT multimers as a first-class use case
- Use coordinate-backed matches instead of copying read data
- Make CLI flattening simple for downstream binaries
