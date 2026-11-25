import os
from Bio import SeqIO

# Master list of all 12 species organized by taxonomic orders
SPECIES_LIST = [
    "Monodelphis_domestica",
    "Perameles_gunnii",
    "Macrotis_lagotis",
    "Thylacinus_cynocephalus",
    "Sarcophilus_harrisii",
    "Sminthopsis_crassicaudata",
    "Antechinus_flavipes",
    "Myrmecobius_fasciatus",
    "Phascolarctos_cinereus",
    "Vombatus_ursinus",
    "Trichosurus_vulpecula",
    "Macropus_eugenii"
]

print("\n" + "="*70)
print(f"{'Species Name':<30} | {'Raw BLAST Hits':<16} | {'Functional ORs':<16}")
print("="*70)

for species in SPECIES_LIST:
    functional_file = f"results/{species}_functional_ORs.fasta"
    raw_file = f"results/{species}_raw_candidates.fasta"
    
    if os.path.exists(functional_file):
        functional_count = len(list(SeqIO.parse(functional_file, "fasta")))
    else:
        functional_count = "Missing file"
        
    if os.path.exists(raw_file):
        raw_count = len(list(SeqIO.parse(raw_file, "fasta")))
    else:
        raw_count = "Missing file"
        
    print(f"{species:<30} | {raw_count:<16} | {functional_count:<16}")

print("="*70 + "\n")
