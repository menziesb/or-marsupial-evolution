import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def paint_global_chromosomes():
    print("\n📊 Compiling expanded 12-species chromosomal abundance layout...")
    
    # Building the complete, verified locus mapping across all 12 marsupial assemblies
    # Synchronized precisely with your true, uninflated 56-family master matrix data values
    data = {
        "Species": [
            "M. domestica (Opossum)", "P. gunnii (Bandicoot)", "M. lagotis (Bilby)",
            "T. cynocephalus (Thylacine)", "S. harrisii (Devil)", "S. crassicaudata (Dunnart)",
            "A. flavipes (Antechinus)", "My. fasciatus (Numbat)", "P. cinereus (Koala)",
            "V. ursinus (Wombat)", "T. vulpecula (Possum)", "M. eugenii (Wallaby)"
        ],
        "Chr2_Other":,
        "Chr2_OR7":,
        "Chr4_Other":,
        "Chr4_OR7":   [11,  10,  10,  44,  11,  8,   9,   5,   4,   8,   11,  7]
    }
    
    df = pd.DataFrame(data)
    
    # Calculate cumulative totals for clean stacked overlay configurations
    df["Chr2_Total"] = df["Chr2_Other"] + df["Chr2_OR7"]
    df["Chr4_Total"] = df["Chr4_Other"] + df["Chr4_OR7"]
    
    # Establish standard academic plot styling
    sns.set_theme(style="ticks", context="notebook")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7.5), sharey=True)
    plt.subplots_adjust(wspace=0.15)
    
    # --- PANEL 1: CHROMOSOME 2 DISTRIBUTION ---
    # Draw total combined value first in red
    sns.barplot(x="Chr2_Total", y="Species", data=df, color="#de2d26", edgecolor="black", linewidth=0.8, ax=ax1)
    # Draw background other families on top in grey to mask the inner segment
    sns.barplot(x="Chr2_Other", y="Species", data=df, color="#cccccc", edgecolor="black", linewidth=0.8, ax=ax1)
    
    ax1.set_title("Chromosome 2 Locus Architecture", fontsize=12, fontweight="bold", pad=15)
    ax1.set_xlabel("Active Loci Volume Count", fontweight="bold")
    ax1.set_ylabel("")
    ax1.grid(axis="x", linestyle="--", alpha=0.5)
    
    # Specific visual marker to highlight the thylacine breakout row
    ax1.text(365, 3, "★ 49 OR7 Loci", va="center", color="#de2d26", fontweight="bold", fontsize=11)
    
    # --- PANEL 2: CHROMOSOME 4 DISTRIBUTION ---
    sns.barplot(x="Chr4_Total", y="Species", data=df, color="#de2d26", edgecolor="black", linewidth=0.8, ax=ax2)
    sns.barplot(x="Chr4_Other", y="Species", data=df, color="#cccccc", edgecolor="black", linewidth=0.8, ax=ax2)
    
    ax2.set_title("Chromosome 4 Locus Architecture", fontsize=12, fontweight="bold", pad=15)
    ax2.set_xlabel("Active Loci Volume Count", fontweight="bold")
    ax2.set_ylabel("")
    ax2.grid(axis="x", linestyle="--", alpha=0.5)
    
    ax2.text(345, 3, "★ 44 OR7 Loci", va="center", color="#de2d26", fontweight="bold", fontsize=11)
    
    # Position ONE unified, clean global legend safely below the chart axis text frames
    colors = ["#de2d26", "#cccccc"]
    labels = ["Verified OR7 Family Loci", "Other Background OR Families"]
    handles = [plt.Rectangle((0,0),1,1, color=c, edgecolor="black", linewidth=0.8) for c in colors]
    
    fig.legend(handles, labels, loc="upper center", ncol=2, 
               bbox_to_anchor=(0.5, -0.04), fontsize=12, frameon=True, facecolor="white", edgecolor="#d1d5db")
    
    sns.despine(left=False, bottom=False)
    
    os.makedirs("results", exist_ok=True)
    output_plot = "results/or7_chromosomal_abundance.png"
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    print(f"📊 Success! Expanded 12-species chromosome map saved cleanly to: {output_plot}\n")

if __name__ == "__main__":
    paint_global_chromosomes()
