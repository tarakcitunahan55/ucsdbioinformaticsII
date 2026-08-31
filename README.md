# UCSD Bioinformatics II — Genome Sequencing

Python implementations of core algorithms from the Genome Sequencing course, covering genome assembly and antibiotic peptide sequencing.

## Contents

### `genomeassembly_week_1/` — Graph-Based Genome Assembly Basics
Building De Bruijn graphs from sequencing reads.

| File | Description |
|---|---|
| `debrujin_from_kmers.txt` | Sample data for De Bruijn graph construction |
| `debrujin_graph_from_kmers.py` | Constructs a De Bruijn graph from a set of k-mers |

### `euler_week_2/` — Eulerian Path/Cycle Genome Reconstruction
Reconstructing genome sequences via Eulerian graph traversal.

| File | Description |
|---|---|
| `eulerian_cycle.py` / `eulerian_path.py` | Finds Eulerian cycles/paths in a graph |
| `euler_path_string_reconstruct.py` | Reconstructs a genome string from an Eulerian path |
| `k_universal_circular_string.py` | Constructs a k-universal circular string |
| `maximal_nonbranching_paths.py` | Extracts contigs via maximal non-branching paths |
| `contigs_from_imperfect_coverage.py` | Assembles contigs under imperfect sequencing coverage |
| `paired_composition.py` | Builds paired k-mer composition of a genome |
| `string_reconstruct_from_ordered_read_...py` | Reconstructs genome from ordered paired reads |
| `string_reconstruct_shuffled_read_pairs.py` | Reconstructs genome from shuffled (unordered) paired reads |

### `antibiotics_week_3/` — Peptide Sequencing & Antibiotics
Identifying peptides from mass spectrometry data.

| File | Description |
|---|---|
| `rna_to_protein.py` / `peptide_encoding.py` | Translates RNA to protein; finds DNA substrings encoding a peptide |
| `theoretical_spectrum_cyclic_peptide.py` / `linear_spectrum.py` | Generates theoretical mass spectra (cyclic/linear peptides) |
| `cyclopeptide_scoring.py` / `linear_peptide_scoring.py` | Scores candidate peptides against an experimental spectrum |
| `cyclopeptide_seq_branch_bound.py` | Branch-and-bound cyclopeptide sequencing |
| `leaderboard_cyclopeptide_seq_score.py` | Leaderboard-based cyclopeptide sequencing for noisy spectra |
| `trim_leaderboard.py` | Trims the leaderboard to top-scoring candidates |
| `spectral_convolution.py` | Computes spectral convolution to infer amino acid masses |
| `count_peptides_given_mass.py` | Counts peptides of a given total mass |
| `convolution_cyclopeptide_seq_score.py` | Combines convolution with sequencing/scoring |

## Topics covered
De Bruijn graphs · Eulerian path/cycle assembly · paired-read reconstruction · cyclopeptide sequencing · mass spectrometry–based antibiotic peptide identification

## Notes
Coursework implementations; shared for portfolio purposes only, DO NOT COPY.
