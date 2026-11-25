import os
from Bio import SeqIO

def index_genome(species_name):
    fasta_path = f"data/genomes/{species_name}.fna"
    fai_path = f"{fasta_path}.fai"
    
    if not os.path.exists(fasta_path):
        print(f"File missing: {fasta_path}")
        return
        
    print(f"Building high-precision index for {species_name}...")
    with open(fai_path, "w") as fai_out:
        for record in SeqIO.parse(fasta_path, "fasta"):
            # Write exact Bedtools columns: Name \t Length
            fai_out.write(f"{record.id}\t{len(record.seq)}\n")
    print(f"Success! Saved {fai_path}")

if __name__ == "__main__":
    species_to_index = [
        "Perameles_gunnii",
        "Phascolarctos_cinereus",
        "Vombatus_ursinus",
        "Trichosurus_vulpecula",
        "Macropus_eugenii"
    ]
    for sp in species_to_index:
        index_genome(sp)
