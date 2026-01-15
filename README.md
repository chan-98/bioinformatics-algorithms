# Bioinformatics Algorithm Implementations

This repository contains Python implementations of five fundamental bioinformatics algorithms in task specific formats, developed for BeSpoke Labs as part of an AI model training contractual project. These scripts demonstrate core computational biology techniques used in sequence analysis, evolutionary studies, and genomics research.

## 1. Sequence Alignment Algorithms

> `alignment.py`

This Python script implements multiple biological sequence alignment algorithms using dynamic programming, including Needleman-Wunsch (global alignment), Smith-Waterman (local alignment), and semi-global alignment with affine gap penalties. It parses FASTA-formatted input files containing protein or DNA sequences along with alignment parameters, applies the specified algorithm with the chosen scoring scheme (simple match/mismatch, BLOSUM62, or PAM250 substitution matrices), and generates comprehensive outputs including aligned sequences in FASTA format, optimal alignment scores, and detailed JSON statistics (percent identity, gap counts, match/mismatch counts). The script handles edge cases like empty sequences and uses deterministic traceback rules to ensure reproducible results, making it suitable for comparative genomics, protein homology analysis, and sequence database searching tasks.

## 2. Multiple Sequence Alignment Statistics

> `msa-statistics-calculation.py`

This Python script analyzes multiple sequence alignments (MSAs) of protein sequences to compute comprehensive conservation statistics and sequence similarity metrics. It parses FASTA-formatted aligned sequences, calculates per-position conservation scores including Shannon entropy and gap fractions, determines conservation categories (conserved, variable, or gap-only) for each alignment column, and generates pairwise identity matrices comparing all sequence pairs while excluding gap-gap positions from calculations. The script produces three output files per alignment: a summary statistics file with overall metrics (number of conserved/variable positions, average pairwise identity), a tab-separated per-position conservation file with entropy and most common residue at each position, and a symmetric identity matrix showing percent identity between all sequence pairs, making it useful for evolutionary analysis, identifying conserved functional regions, and assessing alignment quality.

## 3. Phylogenetic Tree Construction

> `phylogenetic-tree-from-protein-alignments.py`

This Python script constructs phylogenetic trees from multiple sequence alignments using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering algorithm with distance-based methods. It calculates pairwise evolutionary distances between aligned protein sequences using either simple p-distance or Jukes-Cantor corrected distance (which accounts for multiple substitutions at the same site and handles saturation), builds rooted ultrametric trees by iteratively merging the closest sequence clusters, and optionally performs bootstrap resampling to assess topological confidence. The script produces three output files per alignment set: a Newick-formatted phylogenetic tree with branch lengths and optional bootstrap support values, a symmetric distance matrix in TSV format showing all pairwise evolutionary distances, and a JSON metadata file containing alignment statistics, saturated/identical sequence pairs, and average distances, making it useful for phylogenetic analysis, evolutionary studies, and assessing sequence relationships.

## 4. Population Genomics Pipeline

> `bam2newick.py`

This Python script implements a complete population genomics pipeline that processes aligned sequencing reads (BAM files) to call genetic variants, filter them by quality criteria, calculate population genetics statistics, and construct phylogenetic relationships. It performs variant calling using a pileup-based approach with configurable quality thresholds, filters variants based on depth, allele frequency, genotype quality, and missing data rates, then computes key population genetics metrics including nucleotide diversity (π), Tajima's D for detecting selection, observed and expected heterozygosity, and Fst values to quantify genetic differentiation between populations. The pipeline generates six output files: raw and filtered VCF files with called variants, a JSON file with comprehensive population statistics, a symmetric genetic distance matrix showing pairwise differences between all samples, a Newick-formatted phylogenetic tree built using UPGMA clustering, and a human-readable summary report with quality metrics and interpretation guidelines, making it suitable for population structure analysis, evolutionary studies, and conservation genetics applications.

## 5. Bacterial Gene Prediction

> `bacterial-gene-prediction.py`

This Python script implements a prokaryotic gene prediction pipeline that identifies protein-coding genes in bacterial genomic DNA sequences by detecting open reading frames (ORFs) using multiple biological criteria. It searches for ORFs starting with ATG codons and ending with valid stop codons (TAA, TAG, TGA), applies stringent filters including minimum length (≥100 amino acids), Shine-Dalgarno ribosome binding site detection (scoring upstream sequences for the AGGAGGU consensus with up to 2 mismatches allowed), and Codon Adaptation Index (CAI) calculation using E. coli reference codon usage tables to assess expression likelihood. The pipeline analyzes both DNA strands independently to predict the coding strand based on total ORF coverage, resolves overlapping genes by keeping longer ORFs when overlap exceeds 50% of the shorter gene, and generates four output files per sequence: a FASTA file with predicted protein sequences formatted at 60 characters per line, a detailed JSON metadata file with gene coordinates and statistics, a tab-separated summary table comparing all sequences, and a codon usage analysis JSON, making it suitable for bacterial genome annotation, comparative genomics, and identification of protein-coding regions in newly sequenced prokaryotic genomes.
