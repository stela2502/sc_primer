# sc_primer

Reusable single-cell primer/read-structure grammar and detector.

This crate is intended to hold the shared primer logic used by `bam_tide` ONT and Illumina normalizers.
It supports fixed 10x-like structures, BD/Rhapsody whitelist-backed cell IDs, shifted searches, and custom grammar strings.
