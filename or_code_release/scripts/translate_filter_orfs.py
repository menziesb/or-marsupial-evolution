import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    if len(sys.argv) < 3:
        print("Usage: python translate_filter_orfs.py <input.fasta> <output.fasta>")
        sys.exit(1)
        
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    
    print(f"Reading padded genomic chunks from: {input_fasta}")
    functional_records = []
    seen_sequences = set()
    gene_counter = 0

    if not os.path.exists(input_fasta):
        print(f"[ERROR] Target file not found: {input_fasta}")
        return

    for record in SeqIO.parse(input_fasta, "fasta"):
        # Test both original forward strand and reverse complement strand
        strands = [record.seq, record.seq.reverse_complement()]
        
        for strand in strands:
            for frame in range(3):
                length = len(strand) - frame
                trimmed_seq = strand[frame : frame + (length // 3) * 3]
                amino_acids = str(trimmed_seq.translate(to_stop=False))
                
                # Split segments by stop codon boundaries (*)
                orfs = amino_acids.split('*')
                for orf in orfs:
                    # Locate the start codon (M)
                    start_idx = orf.find('M')
                    if start_idx != -1:
                        valid_orf = orf[start_idx:]
                        orf_len = len(valid_orf)
                        
                        # Apply standard mammalian OR length parameters (250-450 aa)
                        if 250 <= orf_len <= 450:
                            if valid_orf not in seen_sequences:
                                seen_sequences.add(valid_orf)
                                gene_counter += 1
                                
                                meta_id = f"gene_id_{gene_counter}_len_{orf_len}"
                                new_rec = SeqRecord(
                                    Seq(valid_orf),
                                    id=meta_id,
                                    description="functional_olfactory_receptor"
                                )
                                functional_records.append(new_rec)

    SeqIO.write(functional_records, output_fasta, "fasta")
    print(f"[SUCCESS] Isolated {len(functional_records)} unique functional OR genes into: {output_fasta}")

if __name__ == "__main__":
    main()
