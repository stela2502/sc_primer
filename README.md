# sc_primer

Reusable single-cell primer/read-structure grammar and detector for `bam_tide`.

The crate is intended to normalize 10x-like and BD/Rhapsody-like read structures before downstream mapping. Preset chemistries are convenience wrappers around the same grammar engine used by custom primer structures.

## API shape

```rust
use sc_primer::{Chemistry, PrimerDetector};

let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3)?;
let matches = detector.detect_all(seq, qual)?;

for primer_match in matches {
    let cell = primer_match.get_cell(seq, qual)?;
    let umi = primer_match.get_umi(seq, qual)?;
    let insert = primer_match.get_insert(seq, qual)?;
}
# Ok::<(), String>(())
```

`PrimerMatch` stores coordinates. `get_cell`, `get_umi`, and `get_insert` slice from normal FASTQ byte strings (`&[u8]`), not 2-bit encoded sequence.

Errors are plain `Result<T, String>` through the `PrimerResult<T>` alias.
