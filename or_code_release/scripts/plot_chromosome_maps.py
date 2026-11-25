import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main():
    coord_file = "results/dasyuromorphia_or_coordinates.csv"
    if not os.path.exists(coord_file):
        print(f"[ERROR] Run the coordinate scraper script first.")
        return
        
    df = pd.read_csv(coord_file)
    if df.empty:
        print("[ERROR] Coordinate data file is completely empty.")
        return
        
    print(f"Total raw loci coordinates loaded: {len(df)}")
    
    # 1. Identify the top 6 densest structural scaffolds/chromosomes across the dataset
    top_6_scaffolds = df['Chromosome_Scaffold'].value_counts().head(6).index.tolist()
    
    # 2. Filter down strictly to these 6 high-density genetic landscapes
    df_filtered = df[df['Chromosome_Scaffold'].isin(top_6_scaffolds)].copy()
    
    print(f"Loci mapped inside the top 6 dense scaffolds: {len(df_filtered)}")
    if df_filtered.empty:
        print("[ERROR] No overlapping data points inside selected scaffolds.")
        return

    # Clean up long scaffold identifiers so they display neatly on the Y-axis
    df_filtered['Display_Scaffold'] = df_filtered['Chromosome_Scaffold'].apply(
        lambda x: str(x)[:15] + "..." if len(str(x)) > 18 else str(x)
    )

    # 3. Initialize multi-panel visualization layout
    plt.figure(figsize=(14, 10))
    sns.set_theme(style="whitegrid")
    
    # 4. Generate the comparative strip plot tracking actual points
    sns.stripplot(
        data=df_filtered, 
        x="Start_Position", 
        y="Display_Scaffold", 
        hue="Species", 
        palette="Set2",
        dodge=True, 
        alpha=0.7, 
        size=5,
        orient="h"
    )
    
    plt.title("Genomic Distribution & Clustering of Functional OR Genes across Key Dasyuromorphia Scaffolds", fontsize=14, fontweight='bold')
    plt.xlabel("Genomic Coordinate Position along Scaffold (Base Pairs)", fontsize=12)
    plt.ylabel("Dense Scaffold/Chromosome Assignment", fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Carnivorous Marsupials")
    plt.tight_layout()
    
    output_img = "results/dasyuromorphia_or_chrom_maps.png"
    plt.savefig(output_img, dpi=300)
    print(f"[SUCCESS] Chromosome density map generated and saved to: {output_img}")

if __name__ == "__main__":
    main()
