#!/bin/bash

# Explicitly mapped to your exact filenames on disk
SPECIES_LIST=(
    "Monodelphis_domestica"
    "Antechinus_flavipes"
    "Meug_1.1"
    "Phascolarctos_cinereus"
    "Sarcophilus_harrisii"
    "Sminthopsis_crassicaudata"
    "Thylacinus_cynocephalus"
    "Vombatus_ursinus"
)

BAIT_MATRIX="data/OR_known_queries.fasta"
THREADS=6 

echo "=== INITIALISING COMPARATIVE MARSUPIAL OR PIPELINE ==="

for SPECIES in "${SPECIES_LIST[@]}"; do
    GENOME_FASTA="data/genomes/${SPECIES}.fna"
    FAI_INDEX="data/genomes/${SPECIES}.fna.fai"
    BLAST_OUT="results/${SPECIES}_blast_hits.txt"
    RAW_CANDIDATES="results/${SPECIES}_raw_candidates.fasta"
    FINAL_FUNCTIONAL="results/${SPECIES}_functional_ORs.fasta"
    
    echo "--------------------------------------------------------"
    echo "Processing Target Species Assembly: ${SPECIES}"
    echo "--------------------------------------------------------"

    if [ ! -f "$GENOME_FASTA" ]; then
        echo "[WARN] Genome file missing: $GENOME_FASTA. Skipping species."
        continue
    fi

    # Step 1: Structural Indexing Verification
    if [ ! -f "$FAI_INDEX" ]; then
        echo "[INFO] Chromosome widths index (.fai) missing. Indexing via samtools..."
        samtools faidx "$GENOME_FASTA"
    fi

    # Step 2: Build local Search Database Index
    echo "[STEP 4] Constructing local BLAST database..."
    makeblastdb -in "$GENOME_FASTA" -dbtype nucl -out "data/genomes/${SPECIES}_db"

    # Step 3: High-Sensitivity Local Alignment Tracking (TBLASTN)
    echo "[STEP 4] Running high-sensitivity tblastn (E <= 10^-4)..."
    tblastn \
        -query "$BAIT_MATRIX" \
        -db "data/genomes/${SPECIES}_db" \
        -out "$BLAST_OUT" \
        -evalue 1e-4 \
        -outfmt "6 qseqid sseqid pident length evalue sstart send" \
        -num_threads $THREADS

    # Step 4: Flanking Region Coordinate Padding via Bedtools
    echo "[STEP 5] Parsing hit coordinates and padding intervals (+/- 4500 bp)..."
    awk -v OFS="\t" '{
        start = ($6 < $7) ? $6 : $7;
        end = ($6 < $7) ? $7 : $6;
        print $2, start, end, $1, ".", ($6 < $7) ? "+" : "-"
    }' "$BLAST_OUT" > "results/${SPECIES}_interim.bed"

    # Apply 4,500 bp structural interval padding bound to chromosome boundaries (.fai)
    bedtools slop -i "results/${SPECIES}_interim.bed" -g "$FAI_INDEX" -b 4500 > "results/${SPECIES}_padded.bed"

    # Extract nucleotide streams from primary core assembly
    echo "[STEP 5] Extracting padded candidate nucleotide profiles..."
    bedtools getfasta -fi "$GENOME_FASTA" -bed "results/${SPECIES}_padded.bed" -s -fo "$RAW_CANDIDATES"

    # Step 5: Six-Frame Translation and Biological Filtering
    echo "[STEP 6] Executing adaptive strand-aware translation parser..."
    if [ -f "scripts/translate_filter_orfs.py" ]; then
        python scripts/translate_filter_orfs.py "$RAW_CANDIDATES" "$FINAL_FUNCTIONAL"
    else
        echo "[ERROR] Custom Biopython parsing script missing from scripts/ directory."
    fi

    # Clean up interim footprint files
    rm -f "results/${SPECIES}_interim.bed" "results/${SPECIES}_padded.bed"
    echo "Finished extraction sequence for ${SPECIES}."
done

echo "=== PIPELINE COMPLETE FOR ALL TARGET SPECIES ==="
