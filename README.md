# sc_primer

`sc_primer` is a small Rust library for detecting single-cell primer/read structures in normal FASTQ-like sequence and quality strings.

It is intended as reusable preprocessing logic for tools such as `bam_tide`, where both ONT and Illumina reads need to be normalized into single-cell molecules before mapping or quantification.

The crate supports two complementary modes:

- chemistry presets for common layouts such as 10x and BD Rhapsody
- custom primer grammars for new or unusual read structures

The important design point is that detection returns coordinate-based `PrimerMatch` records. The original read sequence and quality stay outside the match object.

## Supported chemistry presets

Current presets:

| preset | intended use |
|---|---|
| `tenx-v2` | 10x-style fixed adapter + 16 bp cell barcode + 10 bp UMI |
| `tenx-v3` | 10x-style fixed adapter + 16 bp cell barcode + 12 bp UMI |
| `tenx-v4` | 10x-style fixed adapter + 16 bp cell barcode + 12 bp UMI |
| `bd-v1` | BD Rhapsody v1 layout |
| `bd-v2-96` | BD Rhapsody v2 with 96-block barcode logic |
| `bd-v2-384` | BD Rhapsody v2 with 384-block barcode logic |

The 10x presets are fixed-position after an adapter anchor.

The BD/Rhapsody presets use the `BD_CELL:*` grammar primitive. This is not a simple fixed-length cell slice. It performs BD cell block parsing and whitelist/index-backed cell-id calling.

## Installation

Add the crate to your `Cargo.toml`:

```toml
[dependencies]
sc_primer = "0.1"
```

For local development inside a workspace:

```toml
[dependencies]
sc_primer = { path = "../sc_primer" }
```

## Basic use

Detect the first primer match in a read:

```rust
use sc_primer::{Chemistry, PrimerDetector};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3)?;

    let seq = b"CTACACGACGCTCTTCCGATCTACGTACGTACGTACGTTTGGAACCTTGGTTTTTTTTTTTTTTGATCGATCGATC";
    let qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

    let primer_match = detector.detect_first(seq, qual)?;

    let cell = primer_match.get_cell(seq, qual)?;
    let umi = primer_match.get_umi(seq, qual)?;
    let insert = primer_match.get_insert(seq, qual)?;

    println!("cell: {}", cell.seq_string()?);
    println!("umi: {}", umi.seq_string()?);
    println!("insert: {}", insert.seq_string()?);

    Ok(())
}
```

## Detecting ONT multimers

ONT reads may contain multiple concatenated primer + insert molecules. Use `detect_all` for this case.

```rust
use sc_primer::{Chemistry, PrimerDetector, PrimerMatch};

fn detect_ont_molecules(seq: &[u8], qual: &[u8]) -> sc_primer::PrimerResult<Vec<PrimerMatch>> {
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3)?;
    detector.detect_all(seq, qual)
}
```

Each returned `PrimerMatch` stores coordinates into the original read. You can extract the cell barcode, UMI, and insert sequence later:

```rust
use sc_primer::PrimerMatch;

fn print_molecule(seq: &[u8], qual: &[u8], primer_match: &PrimerMatch) -> sc_primer::PrimerResult<()> {
    let cell = primer_match.get_cell(seq, qual)?;
    let umi = primer_match.get_umi(seq, qual)?;
    let insert = primer_match.get_insert(seq, qual)?;

    println!("cell={} umi={} insert_len={}", cell.seq_string()?, umi.seq_string()?, insert.seq.len());

    Ok(())
}
```

For `detect_all`, the insert of one molecule ends at the start of the next detected primer. The final insert ends at the end of the read.

## Custom grammar

Use a custom grammar when no preset matches the assay.

```rust
use sc_primer::{Grammar, PrimerDetector};

fn main() -> sc_primer::PrimerResult<()> {
    let grammar = Grammar::parse(
        "custom-10x-like",
        "SEARCH:0..2+FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10+INSERT",
    )?;

    let detector = PrimerDetector::from_grammar(grammar);

    Ok(())
}
```

The intended CLI rule for downstream tools is:

```text
if --primer-structure is supplied:
    use custom grammar
else:
    use --chemistry
```

## Grammar primitives

Current grammar primitives:

| primitive | meaning |
|---|---|
| `FIXED:<SEQ>:mm=<N>` | match a fixed sequence with up to `N` mismatches |
| `CELL:<LEN>` | capture a fixed-length cell barcode |
| `UMI:<LEN>` | capture a fixed-length UMI |
| `POLYT:min=<N>` | require at least `N` consecutive `T` bases |
| `INSERT` | mark the start of the biological insert |
| `SKIP:<LEN>` | skip a fixed number of bases |
| `SEARCH:<START>..<END>` | retry the following structure at offsets in the given range |
| `BD_CELL:v1` | BD Rhapsody v1 cell layout and whitelist/index call |
| `BD_CELL:v2.96` | BD Rhapsody v2 96-block layout and whitelist/index call |
| `BD_CELL:v2.384` | BD Rhapsody v2 384-block layout and whitelist/index call |

Example 10x-like grammar:

```text
FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10+INSERT
```

Example BD/Rhapsody grammar:

```text
SEARCH:0..4+BD_CELL:v2.384+POLYT:min=10+INSERT
```

## BD/Rhapsody behavior

BD/Rhapsody detection is represented by `BD_CELL:<version>`.

The BD primitive is expected to handle the old Rustody-style behavior:

- `v1` layout uses three 9 bp cell blocks with longer gaps and an 8 bp UMI
- `v2.96` and `v2.384` use the shorter v2 block layout
- v2 layouts try shifts `0..4` when combined with `SEARCH:0..4`
- cell ids are whitelist/index-backed
- v2.384 final id follows `c1 * 384 * 384 + c2 * 384 + c3 + 1`
- v2.96 uses the corresponding 96-block formula, not the 384 formula

For BD reads, `get_cell(seq, qual)` returns the concatenated BD cell blocks. The numeric BD cell id is available as `primer_match.bd_cell_id`.

## Reverse-complement detection

The detector can report matches in forward or reverse-complement orientation.

`PrimerMatch` stores the orientation. The accessor methods return cell, UMI, and insert sequence in the logical molecule orientation, not merely the raw genomic slice orientation.

```rust
let cell = primer_match.get_cell(seq, qual)?;
let umi = primer_match.get_umi(seq, qual)?;
let insert = primer_match.get_insert(seq, qual)?;
```

## Quality handling

All detection functions require `seq.len() == qual.len()`.

The sequence and quality inputs are ordinary FASTQ-style byte strings:

```rust
let seq: &[u8] = read.seq.as_bytes();
let qual: &[u8] = read.qual.as_bytes();
```

They are not 2-bit encoded.

`get_cell`, `get_umi`, and `get_insert` return both sequence and quality for the requested segment.

## Error handling

Most public APIs return:

```rust
sc_primer::PrimerResult<T>
```

which is an alias for the crate result type.

Typical errors include:

- unknown chemistry
- invalid grammar
- sequence/quality length mismatch
- missing primer match
- invalid segment coordinates
- requesting a cell or UMI segment that the grammar did not define

## Intended use in bam_tide

`bam_tide` can use this crate from both ONT and Illumina normalizers.

For ONT:

- use `detect_all`
- split one concatenated ONT read into multiple normalized molecules
- emit one FASTQ record per detected insert
- use `get_cell`, `get_umi`, and `get_insert` to construct normalized output and side-channel tag tables

For Illumina:

- use `detect_first` on R1 or the relevant index/read structure
- pass R2 to the mapper unchanged when appropriate
- keep STAR/mapper logic independent from cell barcode and UMI parsing

## Development status

This crate is designed to keep common chemistries convenient without locking the project into one fixed primer model.

Presets should stay small and readable. Unusual assays should use custom grammar via `--primer-structure`.
