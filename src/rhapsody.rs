use std::collections::HashMap;

use crate::error::{PrimerError, PrimerResult};

use onehot_dna::{encode_candidates, OneHot9};

const BD_V2_96_C1: &[&[u8; 9]] = &[ 
    b"GTCGCTATA", b"CTTGTACTA", b"CTTCACATA",
    b"ACACGCCGG", b"CGGTCCAGG", b"AATCGAATG", b"CCTAGTATA", b"ATTGGCTAA", b"AAGACATGC",
    b"AAGGCGATC", b"GTGTCCTTA", b"GGATTAGGA", b"ATGGATCCA", b"ACATAAGCG", b"AACTGTATT",
    b"ACCTTGCGG", b"CAGGTGTAG", b"AGGAGATTA", b"GCGATTACA", b"ACCGGATAG", b"CCACTTGGA",
    b"AGAGAAGTT", b"TAAGTTCGA", b"ACGGATATT", b"TGGCTCAGA", b"GAATCTGTA", b"ACCAAGGAC",
    b"AGTATCTGT", b"CACACACTA", b"ATTAAGTGC", b"AAGTAACCC", b"AAATCCTGT", b"CACATTGCA",
    b"GCACTGTCA", b"ATACTTAGG", b"GCAATCCGA", b"ACGCAATCA", b"GAGTATTAG", b"GACGGATTA",
    b"CAGCTGACA", b"CAACATATT", b"AACTTCTCC", b"CTATGAAAT", b"ATTATTACC", b"TACCGAGCA",
    b"TCTCTTCAA", b"TAAGCGTTA", b"GCCTTACAA", b"AGCACACAG", b"ACAGTTCCG", b"AGTAAAGCC",
    b"CAGTTTCAC", b"CGTTACTAA", b"TTGTTCCAA", b"AGAAGCACT", b"CAGCAAGAT", b"CAAACCGCC",
    b"CTAACTCGC", b"AATATTGGG", b"AGAACTTCC", b"CAAAGGCAC", b"AAGCTCAAC", b"TCCAGTCGA",
    b"AGCCATCAC", b"AACGAGAAG", b"CTACAGAAC", b"AGAGCTATG", b"GAGGATGGA", b"TGTACCTTA",
    b"ACACACAAA", b"TCAGGAGGA", b"GAGGTGCTA", b"ACCCTGACC", b"ACAAGGATC", b"ATCCCGGAG",
    b"TATGTGGCA", b"GCTGCCAAT", b"ATCAGAGCT", b"TCGAAGTGA", b"ATAGACGAG", b"AGCCCAATC",
    b"CAGAATCGT", b"ATCTCCACA", b"ACGAAAGGT", b"TAGCTTGTA", b"ACACGAGAT", b"AACCGCCTC",
    b"ATTTAGATG", b"CAAGCAAGC", b"CAAAGTGTG", b"GGCAAGCAA", b"GAGCCAATA", b"ATGTAATGG",
    b"CCTGAGCAA", b"GAGTACATT", b"TGCGATCTA" ];
const BD_V2_96_C2: &[&[u8; 9]] = &[ 
    b"TACAGGATA", b"CACCAGGTA", b"TGTGAAGAA", b"GATTCATCA", b"CACCCAAAG",
    b"CACAAAGGC", b"GTGTGTCGA", b"CTAGGTCCT", b"ACAGTGGTA", b"TCGTTAGCA", b"AGCGACACC",
    b"AAGCTACTT", b"TGTTCTCCA", b"ACGCGAAGC", b"CAGAAATCG", b"ACCAAAATG", b"AGTGTTGTC",
    b"TAGGGATAC", b"AGGGCTGGT", b"TCATCCTAA", b"AATCCTGAA", b"ATCCTAGGA", b"ACGACCACC",
    b"TTCCATTGA", b"TAGTCTTGA", b"ACTGTTAGA", b"ATTCATCGT", b"ACTTCGAGC", b"TTGCGTACA",
    b"CAGTGCCCG", b"GACACTTAA", b"AGGAGGCGC", b"GCCTGTTCA", b"GTACATCTA", b"AATCAGTTT",
    b"ACGATGAAT", b"TGACAGACA", b"ATTAGGCAT", b"GGAGTCTAA", b"TAGAACACA", b"AAATAAATA",
    b"CCGACAAGA", b"CACCTACCC", b"AAGAGTAGA", b"TCATTGAGA", b"GACCTTAGA", b"CAAGACCTA",
    b"GGAATGATA", b"AAACGTACC", b"ACTATCCTC", b"CCGTATCTA", b"ACACATGTC", b"TTGGTATGA",
    b"GTGCAGTAA", b"AGGATTCAA", b"AGAATGGAG", b"CTCTCTCAA", b"GCTAACTCA", b"ATCAACCGA",
    b"ATGAGTTAC", b"ACTTGATGA", b"ACTTTAACT", b"TTGGAGGTA", b"GCCAATGTA", b"ATCCAACCG",
    b"GATGAACTG", b"CCATGCACA", b"TAGTGACTA", b"AAACTGCGC", b"ATTACCAAG", b"CACTCGAGA",
    b"AACTCATTG", b"CTTGCTTCA", b"ACCTGAGTC", b"AGGTTCGCT", b"AAGGACTAT", b"CGTTCGGTA",
    b"AGATAGTTC", b"CAATTGATC", b"GCATGGCTA", b"ACCAGGTGT", b"AGCTGCCGT", b"TATAGCCCT",
    b"AGAGGACCA", b"ACAATATGG", b"CAGCACTTC", b"CACTTATGT", b"AGTGAAAGG", b"AACCCTCGG",
    b"AGGCAGCTA", b"AACCAAAGT", b"GAGTGCGAA", b"CGCTAAGCA", b"AATTATAAC", b"TACTAGTCA",
    b"CAACAACGG" ];
const BD_V2_96_C3: &[&[u8; 9]] = &[ 
    b"AAGCCTTCT", b"ATCATTCTG",
    b"CACAAGTAT", b"ACACCTTAG", b"GAACGACAA", b"AGTCTGTAC", b"AAATTACAG", b"GGCTACAGA",
    b"AATGTATCG", b"CAAGTAGAA", b"GATCTCTTA", b"AACAACGCG", b"GGTGAGTTA", b"CAGGGAGGG",
    b"TCCGTCTTA", b"TGCATAGTA", b"ACTTACGAT", b"TGTATGCGA", b"GCTCCTTGA", b"GGCACAACA",
    b"CTCAAGACA", b"ACGCTGTTG", b"ATATTGTAA", b"AAGTTTACG", b"CAGCCTGGC", b"CTATTAGCC",
    b"CAAACGTGG", b"AAAGTCATT", b"GTCTTGGCA", b"GATCAGCGA", b"ACATTCGGC", b"AGTAATTAG",
    b"TGAAGCCAA", b"TCTACGACA", b"CATAACGTT", b"ATGGGACTC", b"GATAGAGGA", b"CTACATGCG",
    b"CAACGATCT", b"GTTAGCCTA", b"AGTTGCATC", b"AAGGGAACT", b"ACTACATAT", b"CTAAGCTTC",
    b"ACGAACCAG", b"TACTTCGGA", b"AACATCCAT", b"AGCCTGGTT", b"CAAGTTTCC", b"CAGGCATTT",
    b"ACGTGGGAG", b"TCTCACGGA", b"GCAACATTA", b"ATGGTCCGT", b"CTATCATGA", b"CAATACAAG",
    b"AAAGAGGCC", b"GTAGAAGCA", b"GCTATGGAA", b"ACTCCAGGG", b"ACAAGTGCA", b"GATGGTCCA",
    b"TCCTCAATA", b"AATAAACAA", b"CTGTACGGA", b"CTAGATAGA", b"AGCTATGTG", b"AAATGGAGG",
    b"AGCCGCAAG", b"ACAGTAAAC", b"AACGTGTGA", b"ACTGAATTC", b"AAGGGTCAG", b"TGTCTATCA",
    b"TCAGATTCA", b"CACGATCCG", b"AACAGAAAC", b"CATGAATGA", b"CGTACTACG", b"TTCAGCTCA",
    b"AAGGCCGCA", b"GGTTGGACA", b"CGTCTAGGT", b"AATTCGGCG", b"CAACCTCCA", b"CAATAGGGT",
    b"ACAGGCTCC", b"ACAACTAGT", b"AGTTGTTCT", b"AATTACCGG", b"ACAAACTTT", b"TCTCGGTTA",
    b"ACTAGACCG", b"ACTCATACG", b"ATCGAGTCT", b"CATAGGTCA" ];

const BD_V2_386_C1: &[&[u8; 9]] = &[ 
    b"TGTGTTCGC", b"TGTGGCGCC", b"TGTCTAGCG", b"TGGTTGTCC", b"TGGTTCCTC",
    b"TGGTGTGCT", b"TGGCGACCG", b"TGCTGTGGC", b"TGCTGGCAC", b"TGCTCTTCC", b"TGCCTCACC",
    b"TGCCATTAT", b"TGATGTCTC", b"TGATGGCCT", b"TGATGCTTG", b"TGAAGGACC", b"TCTGTCTCC",
    b"TCTGATTAT", b"TCTGAGGTT", b"TCTCGTTCT", b"TCTCATCCG", b"TCCTGGATT", b"TCAGCATTC",
    b"TCACGCCTT", b"TATGTGCAC", b"TATGCGGCC", b"TATGACGAG", b"TATCTCGTG", b"TATATGACC",
    b"TAGGCTGTG", b"TACTGCGTT", b"TACGTGTCC", b"TAATCACAT", b"GTTGTGTTG", b"GTTGTGGCT",
    b"GTTGTCTGT", b"GTTGTCGAG", b"GTTGTCCTC", b"GTTGTATCC", b"GTTGGTTCT", b"GTTGGCGTT",
    b"GTTGGAGCG", b"GTTGCTGCC", b"GTTGCGCAT", b"GTTGCAGGT", b"GTTGCACTG", b"GTTGATGAT",
    b"GTTGATACG", b"GTTGAAGTC", b"GTTCTGTGC", b"GTTCTCTCG", b"GTTCTATAT", b"GTTCGTATG",
    b"GTTCGGCCT", b"GTTCGCGGC", b"GTTCGATTC", b"GTTCCGGTT", b"GTTCCGACG", b"GTTCACGCT",
    b"GTTATCACC", b"GTTAGTCCG", b"GTTAGGTGT", b"GTTAGAGAC", b"GTTAGACTT", b"GTTACCTCT",
    b"GTTAATTCC", b"GTTAAGCGC", b"GTGTTGCTT", b"GTGTTCGGT", b"GTGTTCCAG", b"GTGTTCATC",
    b"GTGTCACAC", b"GTGTCAAGT", b"GTGTACTGC", b"GTGGTTAGT", b"GTGGTACCG", b"GTGGCGATC",
    b"GTGCTTCTG", b"GTGCGTTCC", b"GTGCGGTAT", b"GTGCGCCTT", b"GTGCGAACT", b"GTGCAGCCG",
    b"GTGCAATTG", b"GTGCAAGGC", b"GTCTTGCGC", b"GTCTGGCCG", b"GTCTGAGGC", b"GTCTCAGAT",
    b"GTCTCAACC", b"GTCTATCGT", b"GTCGGTGTG", b"GTCGGAATC", b"GTCGCTCCG", b"GTCCTCGCC",
    b"GTCCTACCT", b"GTCCGCTTG", b"GTCCATTCT", b"GTCCAATAC", b"GTCATGTAT", b"GTCAGTGGT",
    b"GTCAGATAG", b"GTATTAACT", b"GTATCAGTC", b"GTATAGCCT", b"GTATACTTG", b"GTATAAGGT",
    b"GTAGCATCG", b"GTACCGTCC", b"GTACACCTC", b"GTAAGTGCC", b"GTAACAGAG", b"GGTTGTGTC",
    b"GGTTGGCTG", b"GGTTGACGC", b"GGTTCGTCG", b"GGTTCAGTT", b"GGTTATATT", b"GGTTAATAC",
    b"GGTGTACGT", b"GGTGCCGCT", b"GGTGCATGC", b"GGTCGTTGC", b"GGTCGAGGT", b"GGTAGGCAC",
    b"GGTAGCTTG", b"GGTACATAG", b"GGTAATCTG", b"GGCTTGGCC", b"GGCTTCACG", b"GGCTTATGT",
    b"GGCTTACTC", b"GGCTGTCTT", b"GGCTCTGTG", b"GGCTCCGGT", b"GGCTCACCT", b"GGCGTTGAG",
    b"GGCGTGTAC", b"GGCGTGCTG", b"GGCGTATCG", b"GGCGCTCGT", b"GGCGCTACC", b"GGCGAGCCT",
    b"GGCGAGATC", b"GGCGACTTG", b"GGCCTCTTC", b"GGCCTACAG", b"GGCCAGCGC", b"GGCCAACTT",
    b"GGCATTCCT", b"GGCATCCGC", b"GGCATAACC", b"GGCAACGAT", b"GGATGTCCG", b"GGATGAGAG",
    b"GGATCTGGC", b"GGATCCATG", b"GGATAGGTT", b"GGAGTCGTG", b"GGAGAAGGC", b"GGACTCCTT",
    b"GGACTAGTC", b"GGACCGTTG", b"GGAATTAGT", b"GGAATCTCT", b"GGAATCGAC", b"GGAAGCCTC",
    b"GCTTGTAGC", b"GCTTGACCG", b"GCTTCGGAC", b"GCTTCACAT", b"GCTTAGTCT", b"GCTGGATAT",
    b"GCTGGAACC", b"GCTGCGATG", b"GCTGATCAG", b"GCTGAGCGT", b"GCTCTTGTC", b"GCTCTCCTG",
    b"GCTCGGTCC", b"GCTCCAATT", b"GCTATTCGC", b"GCTATGAGT", b"GCTAGTGTT", b"GCTAGGATC",
    b"GCTAGCACT", b"GCTACGTAT", b"GCTAACCTT", b"GCGTTCCGC", b"GCGTGTGCC", b"GCGTGCATT",
    b"GCGTCGGTT", b"GCGTATGTG", b"GCGTATACT", b"GCGGTTCAC", b"GCGGTCTTG", b"GCGGCGTCG",
    b"GCGGCACCT", b"GCGCTGGAC", b"GCGCTCTCC", b"GCGCGGCAG", b"GCGCGATAC", b"GCGCCGACC",
    b"GCGAGCGAG", b"GCGAGAGGT", b"GCGAATTAC", b"GCCTTGCAT", b"GCCTGCGCT", b"GCCTAACTG",
    b"GCCGTCCGT", b"GCCGCTGTC", b"GCCATGCCG", b"GCCAGCTAT", b"GCCAACCAG", b"GCATGGTTG",
    b"GCATCGACG", b"GCAGGCTAG", b"GCAGGACGC", b"GCAGCCATC", b"GCAGATACC", b"GCAGACGTT",
    b"GCACTATGT", b"GCACACGAG", b"GATTGTCAT", b"GATTGGTAG", b"GATTGCACC", b"GATTCTACT",
    b"GATTCGCTT", b"GATTAGGCC", b"GATTACGGT", b"GATGTTGGC", b"GATGTTATG", b"GATGGCCAG",
    b"GATCGTTCG", b"GATCGGAGC", b"GATCGCCTC", b"GATCCTCTG", b"GATCCAGCG", b"GATACACGC",
    b"GAGTTACCT", b"GAGTCGTAT", b"GAGTCGCCG", b"GAGGTGTAG", b"GAGGCATTG", b"GAGCGGACG",
    b"GAGCCTGAG", b"GAGATCTGT", b"GAGATAATT", b"GAGACGGCT", b"GACTTCGTG", b"GACTGTTCT",
    b"GACTCTTAG", b"GACCGCATT", b"GAATTGAGC", b"GAATATTGC", b"GAAGGCTCT", b"GAAGAGACT",
    b"GAACTGCCG", b"GAACGCGTG", b"CTTGTGTAT", b"CTTGTGCGC", b"CTTGTCATG", b"CTTGGTCTT",
    b"CTTGGTACC", b"CTTGGATGT", b"CTTGCTCAC", b"CTTGCAATC", b"CTTGAGGCC", b"CTTGACGGT",
    b"CTTCTGATC", b"CTTCTCGTT", b"CTTCTAGGC", b"CTTCGTTAG", b"CTTATGTCC", b"CTTATGCTT",
    b"CTTATATAG", b"CTTAGGTTG", b"CTTAGGAGC", b"CTTACTTAT", b"CTGTTCTCG", b"CTGTGCCTC",
    b"CTGTCGCAT", b"CTGTCGAGC", b"CTGTAGCTG", b"CTGTACGTT", b"CTGCTTGCC", b"CTGCGTAGT",
    b"CTGCACACC", b"CTGATGGAT", b"CTGAGTCAT", b"CTGACGCCG", b"CTGAACGAG", b"CTCTTGTAG",
    b"CTCTTAGTT", b"CTCTTACCG", b"CTCTGCACC", b"CTCTCGTCC", b"CTCGTATTG", b"CTCGACTAT",
    b"CTCCTGACG", b"CTCACTAGC", b"CTATACGGC", b"CGTTCGCTC", b"CGTTCACCG", b"CGTATAGTT",
    b"CGGTGTTCC", b"CGGTGTCAG", b"CGGTCCTGC", b"CGGCGACTC", b"CGGCACGGT", b"CGGATAGCC",
    b"CGGAGAGAT", b"CGCTAATAG", b"CGCGTTGGC", b"CGCGCAGAG", b"CGCACTGCC", b"CCTTGTCTC",
    b"CCTTGGCGT", b"CCTTCTGAG", b"CCTTCTCCT", b"CCTTCGACC", b"CCTTACTTG", b"CCTGTTCGT",
    b"CCTGTATGC", b"CCTCGGCCG", b"CCGTTAATT", b"CCATGTGCG", b"CCAGTGGTT", b"CCAGGCATT",
    b"CCAGGATCC", b"CCAGCGTTG", b"CATTCCGAT", b"CATTATACC", b"CATGTTGAG", b"ATTGCGTGT",b"ATTGCGGAC", b"ATTGCGCCG", b"ATTGACTTG", b"ATTCGGCTG", b"ATTCGCGAG", b"ATTCCAAGT",
    b"ATTATCTTC", b"ATTACTGTT", b"ATTACACTC", b"ATGTTCTAT", b"ATGTTACGC", b"ATGTGTATC",
    b"ATGTGGCAG", b"ATGTCTGTG", b"ATGGTGCAT", b"ATGCTTACT", b"ATGCTGTCC", b"ATGCTCGGC",
    b"ATGAGGTTC", b"ATGAGAGTG", b"ATCTTGGCT", b"ATCTGTGCG", b"ATCGGTTCC", b"ATCATGCTC",
    b"ATCATCACT", b"ATATCTTAT", b"ATAGGCGCC", b"AGTTGGTAT", b"AGTTGAGCC", b"AGTGCGACC",
    b"AGGTGCTAC", b"AGGCTTGCG", b"AGGCCTTCC", b"AGGCACCTT", b"AGGAATATG", b"AGCGGCCAG",
    b"AGCCTGGTC", b"AGCCTGACT", b"AGCAATCCG", b"AGAGATGTT", b"AGAGAATTC", b"ACTCGCTTG",
    b"ACTCGACCT", b"ACGTACACC", b"ACGGATGGT", b"ACCAGTCTG", b"ACATTCGGC", b"ACATGAGGT",
    b"ACACTAATT" ];

const BD_V2_386_C2: &[&[u8; 9]] = &[ 
            b"TTGTGTTGT", b"TTGTGGTAG",
            b"TTGTGCGGA", b"TTGTCTGTT", b"TTGTCTAAG", b"TTGTCATAT", b"TTGTCACGA", b"TTGTATGAA",
            b"TTGTACAGT", b"TTGGTTAAT", b"TTGGTGCAA", b"TTGGTCGAG", b"TTGGTATTA", b"TTGGCACAG",
            b"TTGGATACA", b"TTGGAAGTG", b"TTGCGGTTA", b"TTGCCATTG", b"TTGCACGCG", b"TTGCAAGGT",
            b"TTGATGTAT", b"TTGATAATT", b"TTGAGACGT", b"TTGACTACT", b"TTGACCGAA", b"TTCTGGTCT",
            b"TTCTGCACA", b"TTCTCCTTA", b"TTCTCCGCT", b"TTCTAGGTA", b"TTCTAATCG", b"TTCGTCGTA",
            b"TTCGTAGAT", b"TTCGGCTTG", b"TTCGGAATA", b"TTCGCCAGA", b"TTCGATTGT", b"TTCGATCAG",
            b"TTCCTCGGT", b"TTCCGGCAG", b"TTCCGCATT", b"TTCCAATTA", b"TTCATTGAA", b"TTCATGCTG",
            b"TTCAGGAGT", b"TTCACTATA", b"TTCAACTCT", b"TTCAACGTG", b"TTATGCGTT", b"TTATGATTG",
            b"TTATCCTGT", b"TTATCCGAG", b"TTATATTAT", b"TTAGGCGCG", b"TTACTGGAA", b"TTACTAGTT",
            b"TTACGTGGT", b"TTACGATAT", b"TTACCTAGA", b"TTACATGAG", b"TTACAGCGT", b"TTACACGGA",
            b"TTACACACT", b"TTAATCAGT", b"TTAATAGGA", b"TTAAGTGTG", b"TTAACCTTG", b"TTAACACAA",
            b"TGTTCACTT", b"TGTTCAAGA", b"TGTTAAGTG", b"TGTGTTATG", b"TGTGTCCAA", b"TGTGGAGCG",
            b"TGTCAGTTA", b"TGTCAGAAG", b"TGGTTAGTT", b"TGGTTACAA", b"TGGCGTTAT", b"TGGCGCCAA",
            b"TGGAGTCTT", b"TGCGTATTG", b"TGATAGAGA", b"TGAGGTATT", b"TGAGAATCT", b"TCTTGGTAA",
            b"TCTTCATAG", b"TCTGTCCTT", b"TCTGGAATT", b"TCTACCGCG", b"TCGTTCGAA", b"TCGTCAGTG",
            b"TCGACGAGA", b"TCATGGCTT", b"TCACACTTA", b"TATTCCGAA", b"TATTATGGT", b"TATGCTATT",
            b"TATCAAGGA", b"TAGTTCAAT", b"TAGCTGCTT", b"TAGAGGAAG", b"TACCTGTTA", b"TACACCTGT",
            b"GTTGTGCGT", b"GTTGGCTAT", b"GTTGCCAAG", b"GTTGACCTT", b"GTTCTGCTA", b"GTTCTGAAT",
            b"GTTCTATCA", b"GTTCGCGTG", b"GTTCCTTAT", b"GTTAGCAGT", b"GTTACTGTG", b"GTTACTCAA",
            b"GTTAAGAGA", b"GTTAACTTA", b"GTGTCGGCA", b"GTGTCCATT", b"GTGCTTGAG", b"GTGCTCGTT",
            b"GTGCTCACA", b"GTGCCTGGA", b"GTCTTGTCG", b"GTCTTGATT", b"GTCTTCCGT", b"GTCTTAAGA",
            b"GTCTCATCT", b"GTCTACGAG", b"GTCGTTGCT", b"GTCGTGTTA", b"GTCGGTAAT", b"GTCGGATGT",
            b"GTCGAGCTG", b"GTCCGGACT", b"GTCCAACAT", b"GTCAGACGA", b"GTCAGAATT", b"GTCACTCTT",
            b"GTCAAGGAA", b"GTATGTCTT", b"GTATGTACA", b"GTATCGGTT", b"GTATATGTA", b"GTATACAAT",
            b"GTAGTTAAG", b"GTAGTCGAT", b"GTAGCCTTA", b"GTAGATACT", b"GTACGATTA", b"GTACAGTCT",
            b"GTAATTCGT", b"GCTTGGCAG", b"GCTTGCTTG", b"GCTTGAGGA", b"GCTTCATTA", b"GCTTATGCG",
            b"GCTGTGTAG", b"GCTGTCATG", b"GCTGGTTGT", b"GCTGGACTG", b"GCTGCCTAA", b"GCTGATATT",
            b"GCTCTTAGT", b"GCTCTATTG", b"GCTCGCCGT", b"GCTCCGCTG", b"GCTATTCTG", b"GCTATACGA",
            b"GCTACTAAG", b"GCTACATGT", b"GCTAACTCT", b"GCGTTGTAA", b"GCGTTCTCT", b"GCGTGCGTA",
            b"GCGTCTTGA", b"GCGTCCGAT", b"GCGTAAGAG", b"GCGCTTACG", b"GCGCGGATT", b"GCGCCATAT",
            b"GCGCATGAA", b"GCGATCAAT", b"GCGAGCCTT", b"GCGAGATTG", b"GCGAGAACA", b"GCCTTGGTA",
            b"GCCTTCTAG", b"GCCTTCACA", b"GCCTGAGTG", b"GCCTCACGT", b"GCCGGCGAA", b"GCCGCACAA",
            b"GCCATGCTT", b"GCCATATAT", b"GCCAATTCG", b"GCATTCGTT", b"GCATGATGT", b"GCAGTTGGA",
            b"GCAGTGTCT", b"GCACTTGTG", b"GCAATCTGT", b"GCAACACTT", b"GATTGTATT", b"GATTGCGAG",
            b"GATTCCAGT", b"GATTCATAT", b"GATTATCAG", b"GATTAGGTT", b"GATGTTGCG", b"GATGGATCT",
            b"GATGCTGAT", b"GATGCCTTG", b"GATCTCCTT", b"GATCGCTTA", b"GATATTGAA", b"GATATTACT",
            b"GAGTGTTAT", b"GAGCTCAGT", b"GAGCGTGCT", b"GAGCGTCGA", b"GAGCGGTTG", b"GAGCGACTT",
            b"GAGCCGAAT", b"GAGATAGAT", b"GAGACCTAT", b"GACGGTCGT", b"GACGCAGGT", b"GACGATATG",
            b"GACCTATCT", b"GAATTAGGA", b"GAATCAGCT", b"GAAGTTCAT", b"GAAGTGGTT", b"GAAGTATTG",
            b"GAAGGCATT", b"GAACGCTGT", b"CTTGTCCAG", b"CTTGGATTG", b"CTTGCTGAA", b"CTTGCCGTG",
            b"CTTGATTCT", b"CTTCTGTCG", b"CTTCGGCGT", b"CTTATGAGT", b"CTTACCGAT", b"CTGTTAGGT",
            b"CTGTCGTCT", b"CTGTATAAT", b"CTGGCTCAT", b"CTGGATGCG", b"CTGCGTGTG", b"CTGCGCGGT",
            b"CTGCCGATT", b"CTGCATTGT", b"CTGATTAAG", b"CTGAGATAT", b"CTGACCTGT", b"CTCGTATCT",
            b"CTCGGCAAG", b"CTCGCAATT", b"CTCCTGCTT", b"CTCCTAAGT", b"CTCCGGATG", b"CTCCGAGCG",
            b"CTCACAGGT", b"CTATTCTAT", b"CTATTAGTG", b"CTATGAATT", b"CTACATATT", b"CGTGGCATT",
            b"CGTCTTAAT", b"CGTCTGGTT", b"CGTCACTGT", b"CGTAGGTCT", b"CGGTTCGAG", b"CGGTTCATT",
            b"CGGTGCTCT", b"CGGTAATTG", b"CGGCCTGAT", b"CGGATATAG", b"CGGAATATT", b"CGCTCCAAT",
            b"CGCGTTCGT", b"CGCAGGTTG", b"CGAGGATGT", b"CGAGCTGTT", b"CGACGGCTT", b"CCTTGTGTG",
            b"CCTGTCTCA", b"CCTGACTAT", b"CCTACCTTG", b"CCGTAGATT", b"CCGGCTGGT", b"CATCGGACG",
            b"CATCGATAA", b"CATCCTTCT", b"CAGTTCTGT", b"CAGTGCCAG", b"CAGGCACTG", b"CAGCCTCTT",
            b"CACTTATAT", b"CACTGGTCG", b"CACTGCATG", b"CACGCGTTG", b"CACGATGTT", b"CACCATCTG",
            b"CACAGGCGT", b"ATTGTACAA", b"ATTGGTATG", b"ATTGCTAAT", b"ATTGCATAG", b"ATTGCAGTT",
            b"ATTCTGCAG", b"ATTCTACGT", b"ATTCGGATT", b"ATTCCGTTG", b"ATTCATCAA", b"ATTCAAGAG",
            b"ATTAGCCTT", b"ATTAATATT", b"ATGTTAGAG", b"ATGTTAACT", b"ATGTAGTCG", b"ATGGTGTAG",
            b"ATGGATTAT", b"ATCTTGAAG", b"ATCTGATAT", b"ATCTCAGAA", b"ATCGCTCAA", b"ATCGCGTCG",
            b"ATCCATGGT", b"ATCATGAGA", b"ATCATAGTT", b"ATCAGCGAG", b"ATCACCATT", b"ATAGTAATT",
            b"ATAGCTGTG", b"ATACTCTCG", b"ATACCTCAT", b"AGTTGCGCG", b"AGTTGAATT", b"AGTTATGAT",
            b"AGTGTCCGT", b"AGTGGCTTG", b"AGTGCTTCT", b"AGTATCATT", b"AGTACACAA", b"AGGTATGCG",
            b"AGGTATAGT", b"AGGCTACTT", b"AGGCCAGGT", b"AGGAGCGAT", b"AGCTTATAG", b"AGCTCTAGA",
            b"AGCGTGTAT", b"AGCGTCACA", b"AGCCTTCAT", b"AGCCTGTCG", b"AGCCTCGAG", b"AGCACTGAA",
            b"AGATGTACG", b"AGAGTTAAT", b"AGACCTCTG", b"ACTTCTATA", b"ACTGTCGAG", b"ACTGTATGT",
            b"ACTCTGTAA", b"ACTCGCGAA", b"ACTAGATCT", b"ACTAACGTT", b"ACGTTACTG", b"ACGTGGAAT",
            b"ACGGACTCT", b"ACGCCTAAT", b"ACGCCGTTA", b"ACGACGTGT", b"ACCTCGCAT", b"ACCATCATA",
            b"ACATATATT", b"ACAGGCACA", b"ACACCTGAG", b"ACACATTCT" ] 
;
const BD_V2_386_C3: &[&[u8; 9]] = &[ 
b"TTGTGGCTG", b"TTGTGGAGT", b"TTGTGCGAC", b"TTGTCTTCA", b"TTGTAAGAT",
            b"TTGGTTCTG", b"TTGGTGCGT", b"TTGGTCTAC", b"TTGGTAACT", b"TTGGCGTGC", b"TTGGATTAG",
            b"TTGGAGACG", b"TTGGAATCA", b"TTGCGGCGA", b"TTGCGCTCG", b"TTGCCTTAC", b"TTGCCGGAT",
            b"TTGCATGCT", b"TTGCACGTC", b"TTGCACCAT", b"TTGAACCTG", b"TTCTCGCGT", b"TTCTCAACT",
            b"TTCTACTCA", b"TTCGTCCAT", b"TTCGGATAC", b"TTCGGACGT", b"TTCGCAATC", b"TTCCGGTGC",
            b"TTCCGACTG", b"TTCATTATG", b"TTCATGGAT", b"TTCAGCGCA", b"TTCACCTCG", b"TTCAAGCAG",
            b"TTCAACTAC", b"TTATGCCAG", b"TTATGCATC", b"TTATCGTAC", b"TTATACCTA", b"TTATAATAG",
            b"TTATAAGTC", b"TTAGTTAGC", b"TTAGCTCAT", b"TTAGCACTA", b"TTAGATATG", b"TTACTACGA",
            b"TTACCGTCA", b"TTACAGAGC", b"TTAATTGCA", b"TTAACAGAT", b"TGTTGGCTA", b"TGTTGATGA",
            b"TGTTAAGCT", b"TGTGGCCGA", b"TGTGCTAGC", b"TGTGCGTCA", b"TGTCGCAGT", b"TGTCGAGCA",
            b"TGTACAACG", b"TGGTTCCGA", b"TGGTTCACT", b"TGGTCAAGT", b"TGGCTTGTA", b"TGGCTGTCG",
            b"TGGCGTATG", b"TGGCGCGCT", b"TGGATGTAC", b"TGGACTTGC", b"TGGAATACT", b"TGCTAGCGA",
            b"TGCGTTGCT", b"TGCGGTCTG", b"TGCGCTTAG", b"TGCGCGACG", b"TGCCTGCAT", b"TGCCTAGAC",
            b"TGCACGAGT", b"TGAGTGTGC", b"TGAGGCTCG", b"TCTTCCGTC", b"TCTTATAGT", b"TCTTACCAT",
            b"TCTGTTGTC", b"TCTGTTACT", b"TCTGGCTAG", b"TCTCAGATC", b"TCTAGTTGA", b"TCTAGTACG",
            b"TCGTACTAC", b"TCGGTGTAG", b"TCGGCTGCT", b"TCGCTACTG", b"TCGATCACG", b"TCGAGGCAT",
            b"TCCGGCGTC", b"TCCGGAGCT", b"TCCGCTCGT", b"TCCGAGTAC", b"TCCATTCAT", b"TCCATGGTC",
            b"TCCAAGTCG", b"TCATTACGT", b"TCATGCACT", b"TCAGGTTGC", b"TCAGACCGT", b"TCACTCAGT",
            b"TCAAGCTCA", b"TATTGCGCA", b"TATTCGGCT", b"TATTCCAGC", b"TATTCATCA", b"TATGTTCAG",
            b"TATGGTATG", b"TATGCAAGT", b"TATCTGGTC", b"TATCTGACT", b"TATCCAGAT", b"TATCAGTCG",
            b"TATCACGCT", b"TAGGCGCGA", b"TAGGCACAT", b"TAGGATCGT", b"TAGCATTGC", b"TAGAGTTAC",
            b"TAGACTGAT", b"TACTTGTCG", b"TACGTCCGA", b"TACCGTACT", b"TACCGCGAT", b"TACCAGGAC",
            b"TACAGAAGT", b"TAAGTGCAT", b"TAAGCTACT", b"GTTGACCGA", b"GTTCTCGAC", b"GTTCCTGCT",
            b"GTTATGATG", b"GTGCTTGCA", b"GTGCCGCGT", b"GTATTGCTG", b"GTATTCCGA", b"GTATTAAGC",
            b"GTATGACGT", b"GTAGTTGTC", b"GTAGTACAT", b"GTAGCTCGA", b"GGTTGCTCA", b"GGTTGAGTA",
            b"GGTTAACGT", b"GGTGTGGCA", b"GGTCTTCAG", b"GGTCGTCTA", b"GGTCGGCGT", b"GGTCCGACT",
            b"GGTCATGTC", b"GGTCACATG", b"GGTAGTGCT", b"GGTAGCGTC", b"GGTACCAGT", b"GGTAAGGAT",
            b"GGCTTGTGC", b"GGCTTGACT", b"GGCTTACGA", b"GGCTGTAGT", b"GGCTGGCAG", b"GGCTCCATC",
            b"GGCGTGGAT", b"GGCGTAATC", b"GGCGCAAGT", b"GGCGAGTAG", b"GGCGACCGT", b"GGCCTGTCA",
            b"GGCCATTGC", b"GGCACTCTG", b"GGATGTCAT", b"GGAGTAACT", b"GGAGAACGA", b"GGACTGGCT",
            b"GGACGTTCA", b"GGAACGTGC", b"GCTGTCCAT", b"GCTGGTTCA", b"GCTGCAACT", b"GCTCGTTAC",
            b"GCTATAGAT", b"GCTAGTCGT", b"GCTACCATG", b"GCGTTCTGA", b"GCGTGTTAG", b"GCGGTATCG",
            b"GCGGAGCAT", b"GCGCGGTGC", b"GCGCCTAGT", b"GCGCCGGCT", b"GCCTTCATG", b"GCCATACTG",
            b"GCATGTTGA", b"GCATGCTAC", b"GCAGTATAC", b"GCAGGTACT", b"GCAGCGCGT", b"GCACCTCAT",
            b"GCAATTCGA", b"GATTGCCGT", b"GATGAACAT", b"GATCTTCGA", b"GATCTGCAT", b"GAGTGGCAT",
            b"GAGTCGGAC", b"GAGTATGAT", b"GAGGCGAGT", b"GAGGCAACG", b"GAGCGCACT", b"GAATAGGCT",
            b"ATTGTCACT", b"ATTGTATCA", b"ATTGGTCAG", b"ATTGGCGAT", b"ATTGATCGT", b"ATTCGTAGT",
            b"ATTCATACG", b"ATTCAGGAC", b"ATTACTTCA", b"ATTAATTAG", b"ATTAAGCAT", b"ATGTCTCTA",
            b"ATGTAGCGT", b"ATGGCATAC", b"ATGGAGATC", b"ATGGACTCG", b"ATGGAACGA", b"ATGCTTCAT",
            b"ATGCTCGCT", b"ATGCGACGT", b"ATGCCGTAG", b"ATGAGTTCG", b"ATGACTATC", b"ATGACCGAC",
            b"ATCTTATGC", b"ATCTTACTA", b"ATCTATCAG", b"ATCGTGTAC", b"ATCGTCTGA", b"ATCGGCATG",
            b"ATCGCGAGC", b"ATCGCAACG", b"ATCGATGCT", b"ATCGAATAG", b"ATCCTTCTG", b"ATCCTGCGT",
            b"ATCCGCACT", b"ATCCATTAC", b"ATCCAAGCA", b"ATCAGATCA", b"ATCACACAT", b"ATCAACGTC",
            b"ATCAACCGA", b"ATATTGAGT", b"ATATTCGTC", b"ATATTACAG", b"ATATCTTGA", b"ATATCGCAT",
            b"ATATCAATC", b"ATAGTCCTG", b"ATAGGTCTA", b"ATAGCTGAC", b"ATAGCGGTA", b"AGTTCGCTG",
            b"AGTTACAGC", b"AGTTAACTA", b"AGTGCAATC", b"AGTCTGGTA", b"AGTCTGAGC", b"AGTCTACAT",
            b"AGTCGAACT", b"AGTCCATCG", b"AGTCATTCA", b"AGTATCCAG", b"AGTAGACTG", b"AGTAATCGA",
            b"AGTAAGTGC", b"AGGTTGGCT", b"AGGTTCTAG", b"AGGTGTTCA", b"AGGTGCCAT", b"AGGTCTGAT",
            b"AGGTCGTAC", b"AGGTCAGCA", b"AGGCTTATC", b"AGGCTATGA", b"AGGCCGACG", b"AGGCCAAGC",
            b"AGGCAGGTC", b"AGGCAAGAT", b"AGGAGCAGT", b"AGGACCGCT", b"AGGAATTAC", b"AGCTTGGAC",
            b"AGCTTAAGT", b"AGCTACACG", b"AGCGTTACG", b"AGCGGTGCA", b"AGCGGAGTC", b"AGCGGACGA",
            b"AGCGCGCTA", b"AGCGATAGC", b"AGCGACTCA", b"AGCCTCTAC", b"AGCCGTCGT", b"AGCATGATC",
            b"AGCACTTCG", b"AGCACGGCA", b"AGATTCTGA", b"AGATTAGAT", b"AGATGATAG", b"AGATATGTA",
            b"AGATACCGT", b"AGAGTGCGT", b"AGAGCCGAT", b"AGACTCACT", b"ACTTGCCTA", b"ACTTGAGCA",
            b"ACTTCTAGC", b"ACTTCGACT", b"ACTTAGTAC", b"ACTGTTGAT", b"ACTGTAACG", b"ACTGGTATC",
            b"ACTGACGTC", b"ACTGAAGCT", b"ACTCTGATG", b"ACTCCTGAC", b"ACTCCGCTA", b"ACTCAACTG",
            b"ACTATTGCA", b"ACTAGGCAG", b"ACTACGCGT", b"ACTAATACT", b"ACGTTCGTA", b"ACGTGTGCT",
            b"ACGTGTATG", b"ACGTGGAGC", b"ACGTCTTCG", b"ACGTCAGTC", b"ACGGTCTCA", b"ACGGTCCGT",
            b"ACGGTACAG", b"ACGGCGCTG", b"ACGCTGCGA", b"ACGCGTGTA", b"ACGCGCCAG", b"ACGATGTCG",
            b"ACGATGGAT", b"ACGATCTAC", b"ACGAGCTGA", b"ACGAGCATC", b"ACGAATCGT", b"ACGAACGCA",
            b"ACCTTGTAG", b"ACCTGTTGC", b"ACCTGTCAT", b"ACCTCGATC", b"ACCTAGGTA", b"ACCTACTGA",
            b"ACCTAATCG", b"ACCGTAGCA", b"ACCGGTAGT", b"ACCGGCTAC", b"ACCGCTTCA", b"ACATTGTGC",
            b"ACATTCTCG", b"ACATGGCTG", b"ACATGACGA", b"ACATATGAT", b"ACATATACG", b"ACAGCGTAC",
            b"ACACTTGCT", b"ACACTATCA", b"ACACGCATG", b"ACACCAGTA", b"ACACCAACT", b"ACACATAGT",
            b"ACACACCTA" ] 
;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BdCellVersion {
    V1,
    V2_96,
    V2_384,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RhapsodyCellCall {
    pub version: BdCellVersion,
    pub cell_id: u64,
    pub cell_seq: Vec<u8>,
    pub cell_qual: Vec<u8>,
    pub umi_seq: Vec<u8>,
    pub umi_qual: Vec<u8>,
    pub shift: usize,
    pub consumed: usize,
    pub c1: (usize, usize),
    pub c2: (usize, usize),
    pub c3: (usize, usize),
    pub umi: (usize, usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RhapsodyWhitelist {
    version: BdCellVersion,
    block_size: u64,

    c1: &'static [ &'static[u8; 9]],
    c2: &'static [ &'static[u8; 9]],
    c3: &'static [ &'static[u8; 9]],

    c1_exact: HashMap<Vec<u8>, u64>,
    c2_exact: HashMap<Vec<u8>, u64>,
    c3_exact: HashMap<Vec<u8>, u64>,

    c1_fuzzy: Vec<OneHot9>,
    c2_fuzzy: Vec<OneHot9>,
    c3_fuzzy: Vec<OneHot9>,
}

impl BdCellVersion {
    pub fn parse(raw: &str) -> PrimerResult<Self> {
        match raw {
            "v1" => Ok(Self::V1),
            "v2.96" => Ok(Self::V2_96),
            "v2.384" => Ok(Self::V2_384),
            other => Err(PrimerError::rhapsody(format!("unknown BD cell version '{other}'"))),
        }
    }

    pub fn block_size(self) -> u64 {
        match self {
            Self::V1 => 96,
            Self::V2_96 => 96,
            Self::V2_384 => 384,
        }
    }

    pub fn umi_len(self) -> usize {
        match self {
            Self::V1 => 8,
            Self::V2_96 | Self::V2_384 => 6,
        }
    }

    pub fn unshifted_consumed_len(self) -> usize {
        match self {
            Self::V1 => 60,
            Self::V2_96 | Self::V2_384 => 42,
        }
    }
}

impl RhapsodyCellCall {
    pub fn empty(version: BdCellVersion) -> Self {
        Self {
            version,
            cell_id: 0,
            cell_seq: Vec::new(),
            cell_qual: Vec::new(),
            umi_seq: Vec::new(),
            umi_qual: Vec::new(),
            shift: 0,
            consumed: 0,
            c1: (0, 0),
            c2: (0, 0),
            c3: (0, 0),
            umi: (0, 0),
        }
    }
}

impl RhapsodyWhitelist {
    pub fn new(
        version: BdCellVersion,
        c1s: &'static [&'static [u8; 9]],
        c2s: &'static [&'static [u8; 9]],
        c3s: &'static [&'static [u8; 9]],
    ) -> Self {
        Self {
            version,
            block_size: version.block_size(),

            c1: c1s,
            c2: c2s,
            c3: c3s,

            c1_exact: Self::make_map(c1s),
            c2_exact: Self::make_map(c2s),
            c3_exact: Self::make_map(c3s),

            c1_fuzzy: encode_candidates::<9, _>(c1s).expect("builtin C1 whitelist must encode"),
            c2_fuzzy: encode_candidates::<9, _>(c2s).expect("builtin C2 whitelist must encode"),
            c3_fuzzy: encode_candidates::<9, _>(c3s).expect("builtin C3 whitelist must encode"),
        }
    }

    fn make_map(entries: &[&[u8; 9]] ) -> HashMap<Vec<u8>, u64> {
        entries
            .into_iter()
            .enumerate()
            .map(|(idx, seq)| (seq.to_vec(), idx as u64))
            .collect()
    }

    pub fn builtin(version: BdCellVersion) -> Self {
        match version {
            BdCellVersion::V1 => Self::bd_v1(),
            BdCellVersion::V2_96 => Self::bd_v2_96(),
            BdCellVersion::V2_384 => Self::bd_v2_384(),
        }
    }

    pub fn bd_v1() -> Self {
        Self::new(
            BdCellVersion::V1,
            BD_V2_96_C1,
            BD_V2_96_C2,
            BD_V2_96_C3,
        )
    }

    pub fn bd_v2_96() -> Self {
        Self::new(
            BdCellVersion::V2_96,
            BD_V2_96_C1,
            BD_V2_96_C2,
            BD_V2_96_C3,
        )
    }

    pub fn bd_v2_384() -> Self {
        Self::new(
            BdCellVersion::V2_384,
            BD_V2_386_C1,
            BD_V2_386_C2,
            BD_V2_386_C3,
        )
    }

    pub fn version(&self) -> BdCellVersion {
        self.version
    }

    pub fn call(
        &self,
        seq: &[u8],
        qual: &[u8],
        offset: usize,
        shift_start: usize,
        shift_end: usize,
    ) -> Option<RhapsodyCellCall> {
        for shift in shift_start..=shift_end {
            if let Some(call) = self.call_exact_shift(seq, qual, offset, shift) {
                return Some(call);
            }
        }
        None
    }

    pub fn cell_id_to_seq(&self, cell_id:u64 ) -> Option<Vec<u8>> {
        if cell_id == 0 {
            return None;
        }

        let id = cell_id - 1;
        let bs = self.block_size as u64;

        let c1_idx = (id / (bs * bs)) as usize;

        let rem = id % (bs * bs);

        let c2_idx = (rem / bs) as usize;
        let c3_idx = (rem % bs) as usize;

        let c1 = self.c1.get(c1_idx)?;
        let c2 = self.c2.get(c2_idx)?;
        let c3 = self.c3.get(c3_idx)?;

        let mut seq = Vec::with_capacity(27);
        seq.extend_from_slice(*c1);
        seq.extend_from_slice(*c2);
        seq.extend_from_slice(*c3);

        Some(seq)
    }

    pub fn call_exact_shift(
        &self,
        seq: &[u8],
        qual: &[u8],
        offset: usize,
        shift: usize,
    ) -> Option<RhapsodyCellCall> {
        let base = offset.checked_add(shift)?;
        let (c1, c2, c3, umi, consumed) = self.coords(base)?;

        if seq.len() < umi.1 || qual.len() < umi.1 {
            return None;
        }

        let c1_idx = self.index_c1(&seq[c1.0..c1.1])?;
        let c2_idx = self.index_c2(&seq[c2.0..c2.1])?;
        let c3_idx = self.index_c3(&seq[c3.0..c3.1])?;

        let cell_id =
            c1_idx * self.block_size * self.block_size +
            c2_idx * self.block_size +
            c3_idx +
            1;

        let mut cell_seq = Vec::with_capacity(27);
        let mut cell_qual = Vec::with_capacity(27);

        self.extend_part(&mut cell_seq, &mut cell_qual, seq, qual, c1);
        self.extend_part(&mut cell_seq, &mut cell_qual, seq, qual, c2);
        self.extend_part(&mut cell_seq, &mut cell_qual, seq, qual, c3);

        Some(RhapsodyCellCall {
            version: self.version,
            cell_id,
            cell_seq,
            cell_qual,
            umi_seq: seq[umi.0..umi.1].to_vec(),
            umi_qual: qual[umi.0..umi.1].to_vec(),
            shift,
            consumed,
            c1,
            c2,
            c3,
            umi,
        })
    }

    pub fn expected_id(&self, c1: u64, c2: u64, c3: u64) -> u64 {
        c1 * self.block_size * self.block_size + c2 * self.block_size + c3 + 1
    }

    #[inline]
    fn index_block(
        seq: &[u8],
        exact: &HashMap<Vec<u8>, u64>,
        fuzzy: &[OneHot9],
        max_mismatches: u32,
    ) -> Option<u64> {
        if let Some(idx) = exact.get(seq) {
            return Some(*idx);
        }

        let obs = OneHot9::from_bytes(seq).ok()?;
        let (idx, _dist) = obs.best_match(fuzzy, max_mismatches)?;

        Some(idx as u64)
    }

    pub fn index_c1(&self, seq: &[u8]) -> Option<u64> {
        Self::index_block(seq, &self.c1_exact, &self.c1_fuzzy, 1)
    }

    pub fn index_c2(&self, seq: &[u8]) -> Option<u64> {
        Self::index_block(seq, &self.c2_exact, &self.c2_fuzzy, 1)
    }

    pub fn index_c3(&self, seq: &[u8]) -> Option<u64> {
        Self::index_block(seq, &self.c3_exact, &self.c3_fuzzy, 1)
    }

    pub fn create_primer(
        &self,
        c1_idx: usize,
        c2_idx: usize,
        c3_idx: usize,
        umi: &[u8],
    ) -> Vec<u8> {

        let (c1s, c2s, c3s) = match self.version {
            BdCellVersion::V1 => (
                BD_V2_96_C1,
                BD_V2_96_C2,
                BD_V2_96_C3,
            ),
            BdCellVersion::V2_96 => (
                BD_V2_96_C1,
                BD_V2_96_C2,
                BD_V2_96_C3,
            ),
            BdCellVersion::V2_384 => (
                BD_V2_386_C1,
                BD_V2_386_C2,
                BD_V2_386_C3,
            ),
        };

        let mut seq = Vec::new();

        match self.version {
            BdCellVersion::V1 => {
                seq.extend_from_slice(c1s[c1_idx]);
                seq.extend_from_slice(b"AAAAAAAAAAAA");
                seq.extend_from_slice(c2s[c2_idx]);
                seq.extend_from_slice(b"AAAAAAAAAAAAA");
                seq.extend_from_slice(c3s[c3_idx]);
                seq.extend_from_slice(umi);
            }

            BdCellVersion::V2_96 | BdCellVersion::V2_384 => {
                seq.extend_from_slice(c1s[c1_idx]);
                seq.extend_from_slice(b"AAAA");
                seq.extend_from_slice(c2s[c2_idx]);
                seq.extend_from_slice(b"AAAA");
                seq.extend_from_slice(c3s[c3_idx]);
                seq.extend_from_slice(b"A");
                seq.extend_from_slice(umi);
            }
        }

        seq
    }

    pub fn coords(
        &self,
        base: usize,
    ) -> Option<(
        (usize, usize),
        (usize, usize),
        (usize, usize),
        (usize, usize),
        usize,
    )> {
        match self.version {
            BdCellVersion::V1 => Some((
                (base, base + 9),
                (base + 21, base + 30),
                (base + 43, base + 52),
                (base + 52, base + 60),
                base + 60,
            )),
            BdCellVersion::V2_96 | BdCellVersion::V2_384 => Some((
                (base, base + 9),
                (base + 13, base + 22),
                (base + 26, base + 35),
                (base + 36, base + 42),
                base + 42,
            )),
        }
    }

    pub fn extend_part(
        &self,
        cell_seq: &mut Vec<u8>,
        cell_qual: &mut Vec<u8>,
        seq: &[u8],
        qual: &[u8],
        range: (usize, usize),
    ) {
        cell_seq.extend_from_slice(&seq[range.0..range.1]);
        cell_qual.extend_from_slice(&qual[range.0..range.1]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn qual(len: usize) -> Vec<u8> {
        vec![40; len]
    }

    #[test]
    fn bd_v2_386_detects_real_r1_like_read_at_shift_1() {
        let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

        let seq = b"TNACGGAGAGATGTGAGCGCCATATGACAGCGGAGCATTGAACCTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let qual = qual(seq.len());

        let call = wl
            .call(seq, &qual, 0, 0, 4)
            .expect("BD v2.384 call should be detected");

        assert_eq!(call.version, BdCellVersion::V2_384);
        assert_eq!(call.shift, 3);
        assert_eq!(&seq[call.c1.0..call.c1.1], b"CGGAGAGAT","expected CGGAGAGAT from {} to {}", call.c1.0, call.c1.1);
        assert_eq!(&seq[call.c2.0..call.c2.1], b"GCGCCATAT","expected GCGCCATAT from {} to {}", call.c2.0, call.c2.1);
        assert_eq!(&seq[call.c3.0..call.c3.1], b"GCGGAGCAT","expected GCGGAGCAT from {} to {}", call.c3.0, call.c3.1);


        assert_eq!(
            call.cell_id,
            45928512,
        );

        assert_eq!(
            call.cell_seq,
            b"CGGAGAGATGCGCCATATGCGGAGCAT".to_vec(),
        );
    }

    #[test]
    fn bd_v2_386_fails_without_required_shift() {
        let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

        let seq = b"TNACGGAGAGATGTGAGCGCCATATGACAGCGGAGCATTGAACCTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let qual = qual(seq.len());

        assert!(
            wl.call(seq, &qual, 0, 0, 0).is_none(),
            "shift 0 should not match this read"
        );
    }

    #[test]
    fn bd_v2_386_detects_real_r1_against_builtin_whitelist() {
        let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

        let seq = b"TNACGGAGAGATGTGAGCGCCATATGACAGCGGAGCATTGAACCTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let qual = vec![40; seq.len()];

        let call = wl
            .call(seq, &qual, 0, 0, 4)
            .expect("BD v2.384 builtin whitelist should detect this read");

        eprintln!("shift: {}", call.shift);
        eprintln!("cell_id: {}", call.cell_id);
        eprintln!("cell_seq: {}", String::from_utf8_lossy(&call.cell_seq));
        eprintln!("umi: {}", String::from_utf8_lossy(&call.umi_seq));

        assert_eq!(call.version, BdCellVersion::V2_384);
        assert_eq!(call.umi_seq.len(), 6);
        assert_eq!(call.cell_seq.len(), 27);
    }
}