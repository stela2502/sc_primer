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
        RhapsodyWhitelist::bd_v2_384().expected_id(7, 11, 19)
    }

    fn expected_bd_id_v2_386_second() -> u64 {
        RhapsodyWhitelist::bd_v2_384().expected_id(23, 29, 31)
    }

    fn expected_bd_id_v1() -> u64 {
        RhapsodyWhitelist::bd_v1().expected_id(7, 11, 19)
    }
}

struct TenxOntMultimerTest;

impl TenxOntMultimerTest {
    fn cell(index: usize) -> Vec<u8> {
        format!("ACGTACGTACGT{:04}", index).into_bytes()
    }

    fn umi(index: usize) -> Vec<u8> {
        format!("TTGGAACC{:04}", index).into_bytes()
    }

    fn insert(index: usize) -> Vec<u8> {
        format!("GATCGATCGATCGATCGATCGATC{:04}", index).into_bytes()
    }

    fn build_read() -> (Vec<u8>, Vec<u8>, Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>) {
        let mut seq = Vec::new();
        let mut qual = Vec::new();
        let mut cells = Vec::new();
        let mut umis = Vec::new();
        let mut inserts = Vec::new();

        for index in 0..10 {
            let cell = Self::cell(index);
            let umi = Self::umi(index);
            let insert = Self::insert(index);
            let monomer_quality = b'!'.saturating_add(index as u8);

            seq.extend_from_slice(TestData::tenx_adapter());
            seq.extend_from_slice(&cell);
            seq.extend_from_slice(&umi);
            seq.extend_from_slice(b"TTTTTTTTTTTTTT");
            seq.extend_from_slice(&insert);

            let added = TestData::tenx_adapter().len() + cell.len() + umi.len() + 14 + insert.len();
            qual.extend(std::iter::repeat(monomer_quality).take(added));

            cells.push(cell);
            umis.push(umi);
            inserts.push(insert);
        }

        (seq, qual, cells, umis, inserts)
    }

    fn run() {
        let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
        let (seq, qual, expected_cells, expected_umis, expected_inserts) = Self::build_read();
        let matches = detector.detect_all(&seq, &qual).unwrap();

        assert_eq!(matches.len(), 10);

        for index in 0..10 {
            let primer_match = &matches[index];
            let cell = primer_match.get_cell(&seq, &qual).unwrap();
            let umi = primer_match.get_umi(&seq, &qual).unwrap();
            let insert = primer_match.get_insert(&seq, &qual).unwrap();
            let expected_quality = vec![b'!'.saturating_add(index as u8); cell.seq.len()];

            assert_eq!(primer_match.chemistry_name, "tenx-v3");
            assert_eq!(primer_match.orientation, Orientation::Forward);
            assert_eq!(cell.seq, expected_cells[index]);
            assert_eq!(umi.seq, expected_umis[index]);
            assert_eq!(insert.seq, expected_inserts[index]);
            assert_eq!(cell.qual, expected_quality);
        }
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
    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.chemistry_name, "tenx-v3");
    assert_eq!(hit.orientation, Orientation::Forward);
    assert_eq!(hit.get_cell(&seq, &qual).unwrap().seq.as_slice(), b"AACCGGTTAACCGGTT");
    assert_eq!(hit.get_umi(&seq, &qual).unwrap().seq.as_slice(), b"TTAACCGGTTAA");
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
    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.primer_start, 0);
    assert_eq!(hit.get_cell(&seq, &qual).unwrap().seq.as_slice(), b"AACCGGTTAACCGGTT");
    assert_eq!(hit.get_umi(&seq, &qual).unwrap().seq.as_slice(), b"TTAACCGGTTAA");
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
    assert!(detector.detect_first(&seq, &qual).unwrap().is_none());
}


#[test]
fn bd_v2_386_index_c1_accepts_one_n_via_fuzzy_rescue() {
    let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

    let exact = wl
        .index_c1(b"CGGAGAGAT")
        .expect("exact C1 should exist");

    let rescued = wl
        .index_c1(b"CGGNGAGAT")
        .expect("C1 with one N should be rescued");

    let rescued = wl
        .index_c1(b"CGGTGAGAT")
        .expect("C1 with one T should be rescued");

    assert_eq!(rescued, exact);
}

#[test]
fn test_bd_v2_386_rejects_shift_outside_zero_to_four() {
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
    assert!(detector.detect_first(&seq, &qual).unwrap().is_none());
}

#[test]
fn test_bd_v2_386_rejects_non_whitelist_cell_block() {
    let mut seq = TestData::bd_v2_shifted_read();
    seq[3] = b'T';
    let qual = TestData::qual(seq.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
    assert!(detector.detect_first(&seq, &qual).unwrap().is_none());
}



#[test]
fn test_reverse_complement_detection_for_tenx_like_read() {
    let forward = TestData::tenx_v3_read();
    let reverse = TestData::rc(&forward);
    let qual = TestData::qual(reverse.len());
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    let hit = detector.detect_first(&reverse, &qual).unwrap().unwrap();
    assert_eq!(hit.orientation, Orientation::ReverseComplement);
    assert!(hit.get_cell(&reverse, &qual).is_ok());
    assert!(hit.get_umi(&reverse, &qual).is_ok());
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
    assert_eq!(hits[0].get_cell(&seq, &qual).unwrap().seq.as_slice(), b"AACCGGTTAACCGGTT");
    assert_eq!(hits[1].get_cell(&seq, &qual).unwrap().seq.as_slice(), b"AACCGGTTAACCGGTT");
    assert_eq!(hits[2].get_cell(&seq, &qual).unwrap().seq.as_slice(), b"AACCGGTTAACCGGTT");
}

#[test]
fn test_detect_all_tenx_ont_read_with_ten_primer_insert_monomers() {
    TenxOntMultimerTest::run();
}

#[test]
fn test_quality_is_sliced_with_cell_and_umi() {
    let seq = TestData::tenx_v3_read();
    let mut qual = TestData::qual(seq.len());
    let cell_start = TestData::tenx_adapter().len();
    let umi_start = cell_start + 16;
    for q in qual.iter_mut().skip(cell_start).take(16) {
        *q = b'!';
    }
    for q in qual.iter_mut().skip(umi_start).take(12) {
        *q = b'#';
    }
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.get_cell(&seq, &qual).unwrap().qual, vec![b'!'; 16]);
    assert_eq!(hit.get_umi(&seq, &qual).unwrap().qual, vec![b'#'; 12]);
}

#[test]
fn test_sequence_quality_length_mismatch_is_error() {
    let seq = TestData::tenx_v3_read();
    let qual = TestData::qual(seq.len() - 1);
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxV3).unwrap();
    assert!(detector.detect_all(&seq, &qual).is_err());
}

#[test]
fn test_bd_cell_version_parser() {
    assert_eq!(BdCellVersion::parse("v1").unwrap(), BdCellVersion::V1);
    assert_eq!(BdCellVersion::parse("v2.96").unwrap(), BdCellVersion::V2_96);
    assert_eq!(BdCellVersion::parse("v2.384").unwrap(), BdCellVersion::V2_384);
    assert!(BdCellVersion::parse("v3").is_err());
}
