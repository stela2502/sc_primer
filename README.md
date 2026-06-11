# sc_primer

[![Crates.io](https://img.shields.io/crates/v/sc_primer.svg)](https://crates.io/crates/sc_primer)
[![Docs.rs](https://docs.rs/sc_primer/badge.svg)](https://docs.rs/sc_primer)

A high-performance Rust library for detecting and extracting single-cell sequencing structures from FASTQ reads.

`sc_primer` provides a unified primer detection framework for:

* 10x Genomics Chromium chemistries
* BD Rhapsody chemistries
* Custom single-cell library layouts
* Illumina-based workflows
* ONT-based workflows

The crate is designed as a reusable building block for downstream tools such as normalizers, demultiplexers, quantifiers, and single-cell preprocessing pipelines.

---

# Features

* Built-in support for common 10x and BD chemistries
* Forward and reverse-complement detection
* Optional mismatch-tolerant ("fuzzy") matching
* Cell barcode extraction
* UMI extraction
* Insert sequence extraction
* PolyT detection
* Multiple molecule detection within a single read
* Reusable CLI integration via Clap
* Custom chemistry definitions

---

# Installation

```toml
[dependencies]
sc_primer = "CURRENT_VERSION"
```

---

# Quick Start

```rust
use sc_primer::{Chemistry, PrimerDetector};

let detector = PrimerDetector::from_chemistry(
    Chemistry::TenxV3
)?;

let hits = detector.detect_all(read_sequence);

for hit in hits {
    println!("{:?}", hit);
}
```

---

# Library Architecture

The crate is intentionally divided into two layers.

## End-User API

Most applications only need these types:

```rust
use sc_primer::{
    Chemistry,
    PrimerCli,
    PrimerDetector,
    PrimerMatch,
};
```

| Type             | Purpose                      |
| ---------------- | ---------------------------- |
| `Chemistry`      | Built-in chemistry presets   |
| `PrimerCli`      | Reusable clap integration    |
| `PrimerDetector` | Main detection engine        |
| `PrimerMatch`    | Result returned by detection |

---

## Internal Grammar Layer

Internally all chemistries are converted into a common grammar representation.

Examples of grammar tokens include:

```text
FIXED
SEARCH
CELL
BD_CELL
SAMPLE
UMI
POLYT
INSERT
```

This allows:

* 10x Genomics
* BD Rhapsody
* Custom chemistries

to share the same detection engine.

---

# Supported Chemistries

## 10x Genomics

* TenxV2
* TenxV3
* TenxV4

## BD Rhapsody

* BdV1
* BdV2_96
* BdV2_384

These presets automatically configure barcode layouts, primer sequences, and extraction logic.

---

# Public API

## Chemistry

```rust
pub enum Chemistry {
    TenxV2,
    TenxV3,
    TenxV4,
    BdV1,
    BdV2_96,
    BdV2_384,
}
```

The exact variants may evolve as additional chemistries are added.

---

## Orientation

```rust
pub enum Orientation {
    Forward,
    ReverseComplement,
}
```

Every detected molecule records the orientation in which it was found.

This is particularly useful for ONT workflows where reads may arrive in either direction.

---

# Primer Detection

The primary entry point is:

```rust
let hits = detector.detect_all(sequence);
```

This scans the sequence and returns every detectable molecule.

Each returned `PrimerMatch` contains:

* coordinates
* orientation
* extracted cell barcode
* extracted UMI
* insert sequence
* chemistry-specific metadata

---

# Reusable CLI Integration

One of the main goals of `sc_primer` is eliminating duplicated chemistry handling across tools.

Instead of implementing chemistry selection in every application, simply embed the shared CLI.

## Application CLI

```rust
#[derive(clap::Parser)]
pub struct Cli {
    #[command(flatten)]
    pub primer: sc_primer::PrimerCli,
}
```

---

## Creating a Detector

```rust
let detector = cli.primer.detector()?;
```

The detector is fully configured regardless of whether the user selected:

```text
--chemistry tenx-v3
```

or

```text
--chemistry bd-v2-384
```

or

```text
--primer-structure ...
```

This keeps chemistry handling centralized inside `sc_primer`.

---

# Custom Chemistries

For novel protocols a detector can be constructed directly from a grammar definition.

Example:

```text
FIXED(ACGT...)
CELL(16)
UMI(12)
POLYT
INSERT
```

This enables rapid prototyping of new sequencing chemistries without modifying application code.

---

# Fuzzy Matching

`sc_primer` supports optional mismatch-tolerant matching of fixed primer sequences.

This allows recovery of molecules containing sequencing errors.

Example:

```text
CTACACGACGCTCTTCCGATCT
```

may still be detected even when a small number of mismatches are present.

## Performance Impact

### 10x

For 10x chemistries the overhead is negligible because only a small number of anchored primer sequences must be checked.

### BD

For BD chemistries fuzzy matching is more expensive because:

* multiple barcode blocks are evaluated
* whitelist matching is performed
* shifted candidate positions are tested
* mismatch scoring must be computed

Use fuzzy matching when sensitivity is more important than maximum throughput.

---

# Benchmarks

Benchmarks were performed using multimer detection (`detect_all()`).

---

## 10x Chromium

### Standard Detection

```text
10x detect_all multimer
8.348 µs/op
0.835 µs/hit
```

### Fuzzy Detection

```text
10x fuzzy detect_all multimer
8.503 µs/op
0.850 µs/hit
```

| Mode     |        Time |
| -------- | ----------: |
| Standard | 8.348 µs/op |
| Fuzzy    | 8.503 µs/op |

Fuzzy matching increases runtime by only ~2%.

---

## BD Rhapsody v2-384

### Standard Detection

```text
BD detect_all multimer
99.514 µs/op
0.995 µs/hit
```

### Fuzzy Detection

```text
BD fuzzy detect_all multimer
176.714 µs/op
1.767 µs/hit
```

| Mode     |          Time |
| -------- | ------------: |
| Standard |  99.514 µs/op |
| Fuzzy    | 176.714 µs/op |

Even with fuzzy matching enabled, recovered molecules remain below 2 µs per hit.

---

# Typical Use Cases

* FASTQ normalization
* ONT barcode recovery
* Illumina preprocessing
* Single-cell demultiplexing
* Cell barcode extraction
* UMI extraction
* Sequencing chemistry prototyping
* Long-read single-cell workflows

---

# Design Goals

* Fast enough for production-scale sequencing datasets
* Shared chemistry implementation across projects
* Minimal downstream integration effort
* Extensible grammar-based architecture
* Consistent behavior across sequencing platforms

---

# License

MIT License.
