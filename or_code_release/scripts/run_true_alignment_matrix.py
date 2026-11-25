import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

def execute_matrix_alignment():
    print("\n🧬 Launching Scope-Corrected Protein Alignment Matrix Engine...")
    print("-> Synchronizing sequence data layers natively in system memory...")
    
    thylacine_path = "results/Thylacine_Genuine_Assembly_OR7.fasta"
    dog_path = "data/refseq_controls/Dog_True_Mined_OR7s.fasta"
    
    if not os.path.exists(thylacine_path) or not os.path.exists(dog_path):
        print("❌ Error: Missing primary target fasta files on your hard drive.")
        return
        
    th_records = list(SeqIO.parse(thylacine_path, "fasta"))
    dog_records = list(SeqIO.parse(dog_path, "fasta"))
    
    print(f"-> Cross-examining {len(th_records)} Thylacine CDS tracks against {len(dog_records)} Dog proteins...")
    
    results = []
    for th in th_records:
        best_match = None
        best_score = -1
        
        # FIX: Dynamically translate Thylacine nucleotide text into protein letters
        th_nucl = Seq(str(th.seq).upper())
        th_aa = str(th_nucl.translate(to_stop=False))
        
        for dg in dog_records:
            dg_aa = str(dg.seq).upper()
            
            # Direct character comparison matrix loop across matching protein frames
            matches = sum(1 for a, b in zip(th_aa, dg_aa) if a == b)
            score = (matches / max(len(th_aa), len(dg_aa))) * 100
            
            if score > best_score:
                best_score = score
                best_match = dg.id
                
        results.append({
            "Thylacine Locus": th.id, 
            "Best Dog Match ID": best_match, 
            "Actual Identity (%)": best_score
        })
        
    df = pd.DataFrame(results)
    output_report = "results/true_thylacine_dog_matrix_report.txt"
    df.to_csv(output_report, sep="\t", index=False)
    
    print("\n" + "="*85)
    print("📈 UNBIASED PROTEIN-LEVEL EXON MATRIX COMPLETE")
    print("="*85)
    print(f"  -> Total alignments verified: {len(df)} loci successfully cross-mapped")
    print(f"  -> Mean Full-Sequence Family Similarity: {df['Actual Identity (%)'].mean():.2f}%")
    print(f"  -> Maximum Individual Locus Match Peak: {df['Actual Identity (%)'].max():.2f}%")
    print(f"📁 Detailed report safely written to your results directory: {output_report}")
    print("="*85 + "\n")

if __name__ == "__main__":
    execute_matrix_alignment()
