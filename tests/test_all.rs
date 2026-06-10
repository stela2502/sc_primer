use pretty_assertions::assert_eq;
use read_tag_table::ReadTagRecord;

use sc_primer::{
    BdCellVersion, Grammar, Orientation, PrimerDetector, RhapsodyWhitelist,
};

struct TestData;

impl TestData {
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

    fn tenx_3p_v3_grammar(name: &str) -> Grammar {
        Grammar::parse(
            name,
            "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10+INSERT",
        )
        .unwrap()
    }

    fn tenx_3p_v3_no_polyt_grammar(name: &str) -> Grammar {
        Grammar::parse(
            name,
            "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+INSERT",
        )
        .unwrap()
    }

    fn tenx_cell() -> &'static [u8] {
        b"AACCGGTTAACCGGTT"
    }

    fn tenx_umi() -> &'static [u8] {
        b"TTAACCGGTTAA"
    }

    fn tenx_insert() -> &'static [u8] {
        b"GACCTGACTGACTGACCTGA"
    }

    fn synthesized_tenx_read(grammar: &Grammar) -> Vec<u8> {
        let mut read = grammar
            .synthesize(Self::tenx_cell(), Self::tenx_umi())
            .unwrap();
        read.extend_from_slice(Self::tenx_insert());
        read
    }

    fn bd_v2_384_grammar(name: &str) -> Grammar {
        Grammar::parse(
            name,
            "FIXED:ATAGGAAACTCATGGT:mm=2+BD_CELL:v2.384+POLYT:min=0+INSERT",
        )
        .unwrap()
    }

    fn bd_v2_384_cell(index: u64) -> Vec<u8> {
        let wl = RhapsodyWhitelist::bd_v2_384();
        let (c1, c2, c3) = wl
            .cell_id_to_parts_ids(index + 1)
            .expect("test requested invalid BD cell id");
        wl.create_cell_cassette(c1, c2, c3)
    }
}

struct TenxOntMultimerTest;

impl TenxOntMultimerTest {
    fn dna_tag(index: usize, len: usize) -> Vec<u8> {
        let alphabet = b"ACGT";
        let mut out = vec![b'A'; len];
        let mut x = index;

        for pos in (0..len).rev() {
            out[pos] = alphabet[x & 3];
            x >>= 2;
        }

        out
    }

    fn cell(index: usize) -> Vec<u8> {
        let mut v = b"ACGTACGTACGT".to_vec();
        v.extend(Self::dna_tag(index, 4));
        v
    }

    fn umi(index: usize) -> Vec<u8> {
        let mut v = b"TTGGAACC".to_vec();
        v.extend(Self::dna_tag(index, 4));
        v
    }

    fn insert(index: usize) -> Vec<u8> {
        let mut v = b"GATCGATCGATCGATCGATCGATC".to_vec();
        v.extend(Self::dna_tag(index, 4));
        v
    }

    fn build_read(
        grammar: &Grammar,
    ) -> (Vec<u8>, Vec<u8>, Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>) {
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

            let mut primer = grammar.synthesize(&cell, &umi).unwrap();
            primer.extend_from_slice(&insert);

            qual.extend(std::iter::repeat(monomer_quality).take(primer.len()));
            seq.extend_from_slice(&primer);

            cells.push(cell);
            umis.push(umi);
            inserts.push(insert);
        }

        (seq, qual, cells, umis, inserts)
    }

    fn run() {
        let grammar = TestData::tenx_3p_v3_no_polyt_grammar("tenx-ont-stress");
        let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();
        let (seq, qual, expected_cells, expected_umis, expected_inserts) =
            Self::build_read(&grammar);

        let matches = detector.detect_all(&seq, &qual).unwrap();

        assert_eq!(matches.len(), 10, "expected 10 matches - got {}", matches.len());

        for index in 0..matches.len() {
            let primer_match = &matches[index];
            let cell = primer_match.get_cell(&seq, &qual).unwrap();
            let umi = primer_match.get_umi(&seq, &qual).unwrap();
            let insert = primer_match.get_insert(&seq, &qual).unwrap();
            let expected_quality = vec![b'!'.saturating_add(index as u8); cell.seq.len()];

            assert_eq!(primer_match.orientation, Orientation::Forward);
            assert_eq!(
                cell.seq, expected_cells[index],
                "{index}: cell seq got {:?} - expected {:?}",
                std::str::from_utf8(&cell.seq),
                std::str::from_utf8(&expected_cells[index])
            );
            assert_eq!(
                umi.seq, expected_umis[index],
                "{index}: umi seq got {:?} - expected {:?}",
                std::str::from_utf8(&umi.seq),
                std::str::from_utf8(&expected_umis[index])
            );
            assert_eq!(
                insert.seq, expected_inserts[index],
                "{index}: insert seq got {:?} - expected {:?}",
                std::str::from_utf8(&insert.seq),
                std::str::from_utf8(&expected_inserts[index])
            );
            assert_eq!(
                cell.qual, expected_quality,
                "{index}: cell qual got {:?} - expected {:?}",
                std::str::from_utf8(&cell.qual),
                std::str::from_utf8(&expected_quality)
            );
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
    assert_eq!(grammar.cell_len(), 16);
    assert_eq!(grammar.umi_len(), 12);
}

#[test]
fn test_tenx_roundtrip_synthesize_detect_polyt_insert() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-roundtrip");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let seq = TestData::synthesized_tenx_read(&grammar);
    let qual = TestData::qual(seq.len());

    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();

    assert_eq!(hit.orientation, Orientation::Forward);
    assert_eq!(
        hit.get_cell(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_cell()
    );
    assert_eq!(
        hit.get_umi(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_umi()
    );
    assert_eq!(
        hit.get_insert(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_insert()
    );
}

#[test]
fn test_custom_10x_adapter_allows_one_fixed_error_after_synthesis() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-fixed-error");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let mut seq = TestData::cat(&[b"NN", &TestData::synthesized_tenx_read(&grammar)]);
    seq[2 + 5] = b'A';

    let qual = TestData::qual(seq.len());
    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();

    assert_eq!(hit.primer_start, 2);
    assert_eq!(
        hit.get_cell(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_cell()
    );
    assert_eq!(
        hit.get_umi(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_umi()
    );
}

#[test]
fn test_custom_10x_adapter_rejects_too_many_fixed_errors() {
    let grammar = Grammar::parse(
        "tenx-fixed-error-reject",
        "FIXED:CTACACGACGCTCTTCCGATCT:mm=1+CELL:16+UMI:12+POLYT:min=10+INSERT",
    )
    .unwrap();
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let mut seq = TestData::synthesized_tenx_read(&grammar);
    seq[5] = b'A';
    seq[6] = b'A';

    let qual = TestData::qual(seq.len());

    assert!(detector.detect_first(&seq, &qual).unwrap().is_none());
}

#[test]
fn bd_v2_384_index_c1_accepts_one_base_fuzzy_rescue() {
    let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

    let exact = wl
        .index_c1(b"CGGAGAGAT")
        .expect("exact C1 should exist");

    let rescued = wl
        .index_c1(b"CGGTGAGAT")
        .expect("C1 with one mismatch should be rescued");

    assert_eq!(rescued, exact);
}

#[test]
fn test_bd_v2_384_roundtrip_synthesize_detect() {
    let grammar = TestData::bd_v2_384_grammar("bd-roundtrip");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let cell = TestData::bd_v2_384_cell(0);
    let umi = b"AACCGG";
    let insert = b"GATCGATCGATC";

    let mut seq = grammar.synthesize(&cell, umi).unwrap();
    seq.extend_from_slice(insert);

    let qual = TestData::qual(seq.len());
    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();

    assert_eq!(hit.bd_cell_id, Some(1));
    assert_eq!(hit.get_umi(&seq, &qual).unwrap().seq.as_slice(), umi);
    assert_eq!(hit.get_insert(&seq, &qual).unwrap().seq.as_slice(), insert);
}

#[test]
fn test_bd_v2_384_invalid_cell_is_translated_to_valid_cell() {
    let grammar = TestData::bd_v2_384_grammar("bd-repair");
    let mut detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let mut source_cell = TestData::bd_v2_384_cell(0);

    // Damage C1 enough that this exact cassette is not valid anymore.
    source_cell[0] = b'T';
    source_cell[1] = b'T';
    source_cell[2] = b'T';

    let record = ReadTagRecord {
        read_id: "read1".to_string(),
        original_read_id: None,
        cell_seq: source_cell.clone(),
        cell_qual: vec![b'I'; source_cell.len()],
        umi_seq: b"AACCGG".to_vec(),
        umi_qual: vec![b'I'; 6],
    };

    let (target_cell, mut primer) = detector.generate(&record).unwrap();
    primer.extend_from_slice(b"GATCGATCGATC");

    let qual = TestData::qual(primer.len());

    let hit = detector.detect_first(&primer, &qual).unwrap().unwrap();

    assert_ne!(target_cell, source_cell);
    assert_eq!(target_cell, TestData::bd_v2_384_cell(0));
    assert_eq!(hit.bd_cell_id, Some(1));

    let translation = detector
        .primer_translation();

    assert_eq!(translation.len(), 1);
}

#[test]
fn test_reverse_complement_detection_for_tenx_like_read() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-rc");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let forward = TestData::synthesized_tenx_read(&grammar);
    let reverse = TestData::rc(&forward);
    let qual = TestData::qual(reverse.len());

    let hit = detector.detect_first(&reverse, &qual).unwrap().unwrap();

    assert_eq!(hit.orientation, Orientation::ReverseComplement);
    assert!(hit.get_cell(&reverse, &qual).is_ok());
    assert!(hit.get_umi(&reverse, &qual).is_ok());
}

#[test]
fn test_detect_all_three_tenx_monomers_with_damaged_junctions() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-three");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let one = TestData::synthesized_tenx_read(&grammar);
    let two = TestData::synthesized_tenx_read(&grammar);
    let three = TestData::synthesized_tenx_read(&grammar);
    let seq = TestData::cat(&[&one, b"NNNCCT", &two, b"GGGGNN", &three]);
    let qual = TestData::qual(seq.len());

    let hits = detector.detect_all(&seq, &qual).unwrap();

    assert_eq!(hits.len(), 3);
    for hit in hits {
        assert_eq!(
            hit.get_cell(&seq, &qual).unwrap().seq.as_slice(),
            TestData::tenx_cell()
        );
    }
}

#[test]
fn test_detect_all_tenx_ont_read_with_ten_primer_insert_monomers() {
    TenxOntMultimerTest::run();
}

#[test]
fn test_quality_is_sliced_with_cell_and_umi() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-quality");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let seq = TestData::synthesized_tenx_read(&grammar);
    let mut qual = TestData::qual(seq.len());

    let fixed_len = b"CTACACGACGCTCTTCCGATCT".len();
    let cell_start = fixed_len;
    let umi_start = cell_start + 16;

    for q in qual.iter_mut().skip(cell_start).take(16) {
        *q = b'!';
    }
    for q in qual.iter_mut().skip(umi_start).take(12) {
        *q = b'#';
    }

    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();

    assert_eq!(hit.get_cell(&seq, &qual).unwrap().qual, vec![b'!'; 16]);
    assert_eq!(hit.get_umi(&seq, &qual).unwrap().qual, vec![b'#'; 12]);
}

#[test]
fn test_sequence_quality_length_mismatch_is_error() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-length-error");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let seq = TestData::synthesized_tenx_read(&grammar);
    let qual = TestData::qual(seq.len() - 1);

    assert!(detector.detect_all(&seq, &qual).is_err());
}

#[test]
fn test_bd_cell_version_parser() {
    assert_eq!(BdCellVersion::parse("v1").unwrap(), BdCellVersion::V1);
    assert_eq!(BdCellVersion::parse("v2.96").unwrap(), BdCellVersion::V2_96);
    assert_eq!(BdCellVersion::parse("v2.384").unwrap(), BdCellVersion::V2_384);
    assert!(BdCellVersion::parse("v3").is_err());
}