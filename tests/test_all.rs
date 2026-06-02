use pretty_assertions::assert_eq;
use sc_primer::{BdCellVersion, Chemistry, Grammar, Orientation, PrimerDetector, RhapsodyWhitelist};

struct TestData;

impl TestData {
    fn tenx_adapter() -> &'static [u8] {
        b"CTACACGACGCTCTTCCGATCT"
    }

    fn qual(len: usize) -> Vec<u8> {
        vec![b'I'; len]
    }

    fn cat(parts: &[&[u8]]) -> Vec<u8> {
        let mut out = Vec::new();
        for part in parts {
            out.extend_from_slice(part);
        }
        out
    }

    fn rc(seq: &[u8]) -> Vec<u8> {
        PrimerDetector::reverse_complement(seq)
    }

    fn tenx_v3_read() -> Vec<u8> {
        Self::cat(&[
            Self::tenx_adapter(),
            b"AACCGGTTAACCGGTT",
            b"TTAACCGGTTAA",
            b"TTTTTTTTTTTTTTTTTT",
            b"GACCTGACTGACTGACCTGA",
        ])
    }

    fn tenx_adapter_read_with_one_adapter_error() -> Vec<u8> {
        Self::cat(&[
            b"NN",
            b"CTACACGACGCTCTTCCGATCT",
            b"AACCGGTTAACCGGTT",
            b"TTAACCGGTTAA",
            b"TTTTTTTTTTT",
            b"ATGGCCAATTGG",
        ])
    }

    fn bd_v2_shifted_read() -> Vec<u8> {
        Self::cat(&[
            b"NNN",
            b"ACGTACGTA",
            b"GGGG",
            b"TGCATGCAT",
            b"CCCC",
            b"GATTACAGA",
            b"A",
            b"AACCGG",
            b"TTTTTTTTTTTT",
            b"GATCGATCGATC",
        ])
    }

    fn bd_v2_unshifted_read() -> Vec<u8> {
        Self::cat(&[
            b"ACGTACGTA",
            b"GGGG",
            b"TGCATGCAT",
            b"CCCC",
            b"GATTACAGA",
            b"A",
            b"AACCGG",
            b"TTTTTTTTTTTT",
            b"GATCGATCGATC",
        ])
    }

    fn bd_v2_second_read() -> Vec<u8> {
        Self::cat(&[
            b"NN",
            b"CCCCGGGGA",
            b"TTTT",
            b"AAAACCCCG",
            b"AAAA",
            b"TTTTGGGGA",
            b"C",
            b"GGTTAA",
            b"TTTTTTTTTTT",
            b"CCGGAATT",
        ])
    }

    fn bd_v1_read() -> Vec<u8> {
        Self::cat(&[
            b"ACGTACGTA",
            b"NNNNNNNNNNNN",
            b"TGCATGCAT",
            b"NNNNNNNNNNNNN",
            b"GATTACAGA",
            b"TTCCAAGG",
            b"TTTTTTTTTTTT",
            b"GGAATTCC",
        ])
    }

    fn expected_bd_id_v2_384() -> u64 {
        RhapsodyWhitelist::toy_v2_384().expected_id(7, 11, 19)
    }

    fn expected_bd_id_v2_384_second() -> u64 {
        RhapsodyWhitelist::toy_v2_384().expected_id(23, 29, 31)
    }

    fn expected_bd_id_v1() -> u64 {
        RhapsodyWhitelist::toy_v1().expected_id(7, 11, 19)
    }

    fn assert_cell(hit: &sc_primer::PrimerMatch, seq: &[u8], qual: &[u8], expected: &[u8]) {
        let part = hit.get_cell(seq, qual).unwrap();
        assert_eq!(part.seq.as_slice(), expected);
    }

    fn assert_umi(hit: &sc_primer::PrimerMatch, seq: &[u8], qual: &[u8], expected: &[u8]) {
        let part = hit.get_umi(seq, qual).unwrap();
        assert_eq!(part.seq.as_slice(), expected);
    }

    fn assert_insert(hit: &sc_primer::PrimerMatch, seq: &[u8], qual: &[u8], expected: &[u8]) {
        let part = hit.get_insert(seq, qual).unwrap();
        assert_eq!(part.seq.as_slice(), expected);
    }

    fn tenx_cell(index: usize) -> String {
        format!("AACCGGTTAACC{:04}", index)
    }

    fn tenx_umi(index: usize) -> String {
        format!("TTAACCGG{:04}", index)
    }

    fn tenx_insert(index: usize) -> String {
        format!("GACCTGACTGACTGACCTGA{:04}", index)
    }

    fn tenx_monomer(index: usize) -> Vec<u8> {
        let cell = Self::tenx_cell(index);
        let umi = Self::tenx_umi(index);
        let insert = Self::tenx_insert(index);
        Self::cat(&[
            Self::tenx_adapter(),
            cell.as_bytes(),
            umi.as_bytes(),
            b"TTTTTTTTTTTT",
            insert.as_bytes(),
        ])
    }
}

#[test]
fn test_parse_complicated_custom_grammar() {
    let grammar = Grammar::parse(
        "custom",
        "SEARCH:0..4+FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10+INSERT",
    )
    .unwrap();
    assert_eq!(grammar.ops.len(), 6);
}

#[test]
fn test_tenx_v3_fixed_position_cell_umi_polyt_insert() {
    let seq = TestData::tenx_v3_read();
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    let hit = detector.detect(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.chemistry_name, "tenx-v3");
    assert_eq!(hit.orientation, Orientation::Forward);
    TestData::assert_cell(&hit, &seq, &qual, b"AACCGGTTAACCGGTT");
    TestData::assert_umi(&hit, &seq, &qual, b"TTAACCGGTTAA");
    assert_eq!(hit.insert_start, 68);
}

#[test]
fn test_custom_10x_adapter_allows_one_late_error_and_two_base_search() {
    let mut seq = TestData::tenx_adapter_read_with_one_adapter_error();
    seq[2] = b'A';
    let qual = TestData::qual(seq.len());
    let grammar = Grammar::parse(
        "tenx-v3-with-adapter",
        "SEARCH:0..4+FIXED:CTACACGACGCTCTTCCGATCT:mm=1+CELL:16+UMI:12+POLYT:min=10+INSERT",
    )
    .unwrap();
    let detector = PrimerDetector::from_grammar(grammar);
    let hit = detector.detect(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.primer_start, 0);
    TestData::assert_cell(&hit, &seq, &qual, b"AACCGGTTAACCGGTT");
    TestData::assert_umi(&hit, &seq, &qual, b"TTAACCGGTTAA");
    assert_eq!(hit.insert_start, 63);
}

#[test]
fn test_custom_10x_adapter_rejects_too_many_adapter_errors() {
    let mut seq = TestData::tenx_adapter_read_with_one_adapter_error();
    seq[2] = b'A';
    seq[3] = b'A';
    let qual = TestData::qual(seq.len());
    let grammar = Grammar::parse(
        "tenx-v3-with-adapter",
        "SEARCH:0..4+FIXED:CTACACGACGCTCTTCCGATCT:mm=1+CELL:16+UMI:12+POLYT:min=10+INSERT",
    )
    .unwrap();
    let detector = PrimerDetector::from_grammar(grammar);
    assert!(detector.detect(&seq, &qual).unwrap().is_none());
}

#[test]
fn test_tenx_v3_rejects_short_polyt() {
    let seq = TestData::cat(&[
        TestData::tenx_adapter(),
        b"AACCGGTTAACCGGTT",
        b"TTAACCGGTTAA",
        b"TTTTTTTTT",
        b"GACCTGACTGACTGACCTGA",
    ]);
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    assert!(detector.detect(&seq, &qual).unwrap().is_none());
}

#[test]
fn test_bd_v2_384_shifted_rhapsody_whitelist_call() {
    let seq = TestData::bd_v2_shifted_read();
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
    let hit = detector.detect(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.bd_cell_id, Some(TestData::expected_bd_id_v2_384()));
    TestData::assert_cell(&hit, &seq, &qual, b"ACGTACGTATGCATGCATGATTACAGA");
    TestData::assert_umi(&hit, &seq, &qual, b"AACCGG");
    assert_eq!(hit.insert_start, 57);
}

#[test]
fn test_bd_v2_384_unshifted_still_works() {
    let seq = TestData::bd_v2_unshifted_read();
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
    let hit = detector.detect(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.bd_cell_id, Some(TestData::expected_bd_id_v2_384()));
    assert_eq!(hit.insert_start, 54);
}

#[test]
fn test_bd_v2_384_rejects_shift_outside_zero_to_four() {
    let seq = TestData::cat(&[
        b"NNNNN",
        b"ACGTACGTA",
        b"GGGG",
        b"TGCATGCAT",
        b"CCCC",
        b"GATTACAGA",
        b"A",
        b"AACCGG",
        b"TTTTTTTTTTTT",
        b"GATCGATCGATC",
    ]);
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
    assert!(detector.detect(&seq, &qual).unwrap().is_none());
}

#[test]
fn test_bd_v2_384_rejects_non_whitelist_cell_block() {
    let mut seq = TestData::bd_v2_shifted_read();
    seq[3] = b'T';
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
    assert!(detector.detect(&seq, &qual).unwrap().is_none());
}

#[test]
fn test_bd_v2_96_uses_96_block_formula_not_384_formula() {
    let seq = TestData::bd_v2_shifted_read();
    let qual = TestData::qual(seq.len());
    let grammar = Grammar::parse("bd-v2-96", "SEARCH:0..4+BD_CELL:v2.96+POLYT:min=10+INSERT").unwrap();
    let detector = PrimerDetector::from_grammar_with_rhapsody(grammar, RhapsodyWhitelist::toy_v2_96());
    let hit = detector.detect(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.bd_cell_id, Some(RhapsodyWhitelist::toy_v2_96().expected_id(7, 11, 19)));
    assert_ne!(hit.bd_cell_id, Some(TestData::expected_bd_id_v2_384()));
}

#[test]
fn test_bd_v1_old_layout_with_longer_gaps_and_8bp_umi() {
    let seq = TestData::bd_v1_read();
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV1).unwrap();
    let hit = detector.detect(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.bd_cell_id, Some(TestData::expected_bd_id_v1()));
    TestData::assert_umi(&hit, &seq, &qual, b"TTCCAAGG");
    assert_eq!(hit.insert_start, 72);
}

#[test]
fn test_reverse_complement_detection_for_tenx_like_read() {
    let forward = TestData::tenx_v3_read();
    let reverse = TestData::rc(&forward);
    let qual = TestData::qual(reverse.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    let hit = detector.detect(&reverse, &qual).unwrap().unwrap();
    assert_eq!(hit.orientation, Orientation::ReverseComplement);
    TestData::assert_cell(&hit, &reverse, &qual, b"AACCGGTTAACCGGTT");
    TestData::assert_umi(&hit, &reverse, &qual, b"TTAACCGGTTAA");
}

#[test]
fn test_detect_all_two_bd_monomers_in_one_ont_multimer() {
    let first = TestData::bd_v2_shifted_read();
    let last = TestData::bd_v2_unshifted_read();
    let second = TestData::bd_v2_second_read();
    let seq = TestData::cat(&[&first, b"NNNNNN", &second,b"NNNACGTN",&last]);
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
    let hits = detector.detect_all(&seq, &qual).unwrap();
    assert_eq!(hits.len(), 3);
    assert_eq!(hits[0].bd_cell_id, Some(TestData::expected_bd_id_v2_384()));
    assert_eq!(hits[1].bd_cell_id, Some(TestData::expected_bd_id_v2_384_second()));
    assert!(hits[1].primer_start > hits[0].insert_start);
    assert!(hits[2].primer_start > hits[1].insert_start);
}

#[test]
fn test_detect_all_three_tenx_monomers_with_damaged_junctions() {
    let one = TestData::tenx_v3_read();
    let two = TestData::tenx_v3_read();
    let three = TestData::tenx_v3_read();
    let seq = TestData::cat(&[&one, b"NNNCCT", &two, b"GGGGNN", &three]);
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    let hits = detector.detect_all(&seq, &qual).unwrap();
    assert_eq!(hits.len(), 3);
    TestData::assert_cell(&hits[0], &seq, &qual, b"AACCGGTTAACCGGTT");
    TestData::assert_cell(&hits[1], &seq, &qual, b"AACCGGTTAACCGGTT");
    TestData::assert_cell(&hits[2], &seq, &qual, b"AACCGGTTAACCGGTT");
}

#[test]
fn test_detect_all_tenx_ont_read_with_ten_primer_insert_monomers() {
    let mut seq = Vec::new();
    let mut expected_cells = Vec::new();
    let mut expected_umis = Vec::new();
    let mut expected_inserts = Vec::new();

    for index in 0..10 {
        seq.extend_from_slice(&TestData::tenx_monomer(index));
        expected_cells.push(TestData::tenx_cell(index));
        expected_umis.push(TestData::tenx_umi(index));
        expected_inserts.push(TestData::tenx_insert(index));
    }

    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    let matches = detector.detect_all(&seq, &qual).unwrap();

    assert_eq!(matches.len(), 10);

    for index in 0..10 {
        let primer_match = &matches[index];
        assert_eq!(primer_match.chemistry_name, "tenx-v3");
        assert_eq!(primer_match.orientation, Orientation::Forward);
        assert_eq!(primer_match.get_cell(&seq, &qual).unwrap().seq.as_slice(), expected_cells[index].as_bytes());
        assert_eq!(primer_match.get_umi(&seq, &qual).unwrap().seq.as_slice(), expected_umis[index].as_bytes());
        assert_eq!(primer_match.get_insert(&seq, &qual).unwrap().seq.as_slice(), expected_inserts[index].as_bytes());
        assert_eq!(primer_match.get_cell(&seq, &qual).unwrap().qual, vec![b'I'; 16]);
        assert_eq!(primer_match.get_umi(&seq, &qual).unwrap().qual, vec![b'I'; 12]);
        assert_eq!(primer_match.get_insert(&seq, &qual).unwrap().qual, vec![b'I'; expected_inserts[index].len()]);
    }
}

#[test]
fn test_quality_is_sliced_with_cell_and_umi() {
    let seq = TestData::tenx_v3_read();
    let mut qual = TestData::qual(seq.len());
    let cell_start = TestData::tenx_adapter().len();
    let umi_start = cell_start + 16;
    for i in cell_start..cell_start + 16 {
        qual[i] = b'!';
    }
    for i in umi_start..umi_start + 12 {
        qual[i] = b'#';
    }
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    let hit = detector.detect(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.get_cell(&seq, &qual).unwrap().qual, vec![b'!'; 16]);
    assert_eq!(hit.get_umi(&seq, &qual).unwrap().qual, vec![b'#'; 12]);
}

#[test]
fn test_sequence_quality_length_mismatch_is_error() {
    let seq = TestData::tenx_v3_read();
    let qual = TestData::qual(seq.len() - 1);
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    assert!(detector.detect(&seq, &qual).is_err());
}

#[test]
fn test_bd_cell_version_parser() {
    assert_eq!(BdCellVersion::parse("v1").unwrap(), BdCellVersion::V1);
    assert_eq!(BdCellVersion::parse("v2.96").unwrap(), BdCellVersion::V2_96);
    assert_eq!(BdCellVersion::parse("v2.384").unwrap(), BdCellVersion::V2_384);
    assert!(BdCellVersion::parse("v3").is_err());
}
