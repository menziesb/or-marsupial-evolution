# Comprehensive analysis of marsupial olfactory receptor (OR) subgenomes

This repository contains the bioinformatic pipeline used to extract, filter, classify, map, and phylogenetically analyze functional olfactory receptor (OR) repertoires across 15 marsupial genomes.

## Pipeline Components

### 1. Genome Indexing & Subgenome Extraction
- `scripts/download_genomes.py`: Retrieves marsupial reference genome assemblies via NCBI Datasets.
- `scripts/generate_fai.py`: Generates FASTA sequence index (`.fai`) files for chromosomal coordinate padding.
- `run_pipeline.sh`: Executes `tblastn` ($E \le 10^{-4}$) against multi-family mammalian OR bait consensus sequences and pads candidate hits by $\pm 4,500$ bp using `bedtools slop` and `bedtools getfasta`.
- `scripts/translate_filter_orfs.py` & `scripts/filter_intact_genes.py`: Performs strand-aware six-frame translation and isolates intact single-exon ORFs ($250\text{--}450\text{ aa}$) lacking internal stop codons.

### 2. Automated JSON BLAST Classification & Family Binning
- `scripts/allocate_families.py`, `scripts/classify_all_marsupials.py`, & `scripts/profile_all_families.py`: Parses NCBI Web-BLAST JSON outputs, filters hits at $\ge 70\%$ sequence identity and $\ge 70\%$ query coverage, deduplicates overlapping coordinates, and assigns sequences to canonical OR families.
- `scripts/generate_summary.py`: Compiles species-by-family functional OR repertoire counts (Table 1).

### 3. OR14 Subfamily Isolation & Chromosomal Mapping
- `harvest_all_marsupial_or14s.py`, `true_bilby_coords.py`, `map_true_bilby_clusters.py`, `precise_id_thylacine_mapper.py`, `run_mole_final.py`, & `parse_functional_mole_mappings.py`: Isolates OR14 loci and maps deduplicated physical start/end coordinates (Mb) across chromosomes and super-scaffolds.
- `scripts/plot_chromosome_maps.py` & `scripts/plot_global_chromosomes.py`: Plots chromosomal tandem duplication arrays (Figs. 2 and 4).

### 4. Phylogenetics & Pairwise Identity Validation
- `scripts/run_global_ml_tree.py`: Aligns intact OR14 amino acid sequences with human Class I/II reference ORs using `MAFFT` (`--globalpair --maxiterate 1000`) and infers Maximum-Likelihood trees in `IQ-TREE 2` (`-m MFP -alrt 1000 -bb 1000`) (Figs. 1 and 3).
- `scripts/run_true_alignment_matrix.py`: Computes pairwise global sequence identities (`BLOSUM62`, Needleman-Wunsch) between orthologous and divergent Peramelemorphia OR14 sequences and human Class II ORs (Table 2).
