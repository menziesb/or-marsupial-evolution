import os
import json
import csv
import re
import sys
from Bio import SeqIO

def classify_description(text):
    if not text: return "Uncharacterised"
    txt = text.upper()
    fams = [str(i) for i in range(1, 16)] + ["51", "52", "55", "56"]
    match = re.search(r"(FAMILY|OR|SUBFAMILY|RECEPTOR)\s*([0-9]+)", txt)
    if match and match.group(2) in fams: return "Family " + match.group(2)
    for f in fams:
        if "OR" + f in txt or "FAMILY " + f in txt: return "Family " + f
    return "Uncharacterised"

def main():
    print("=== NATIVE HIGH-STRINGENCY ALLOCATION MATRIX ENGINE ===")
    dunnart_dir = "results/JSON/Sminthopsis_crassicaudata"
    master_fasta = "results/GCF_048593235.1_ASM4859323v1_genomic_Sminthopsis_functional_ORs.fasta"
    output_fasta_path = "results/dunnart_uncharacterised_only.fasta"
    output_matrix_path = "results/dunnart_high_stringency_only.csv"

    if not os.path.exists(dunnart_dir) or not os.path.exists(master_fasta):
        print("[ERROR] Project pathways are missing on disk."); return

    all_categories = ["Family " + str(i) for i in range(1, 15)] + ["Family 51", "Family 52", "Family 55", "Family 56", "Uncharacterised"]
    counts = {c: 0 for c in all_categories}

    all_records = list(SeqIO.parse(master_fasta, "fasta"))
    extracted = []
    processed = set()
    idx_counter = -1

    json_files = sorted([f for f in os.listdir(dunnart_dir) if f.endswith(".json")])
    for jf in json_files:
        with open(os.path.join(dunnart_dir, jf), "r") as f: data = json.load(f)
        bo2 = data.get("BlastOutput2", [data])
        if isinstance(bo2, dict): bo2 = [bo2]
        for entry in bo2:
            if isinstance(entry, str): continue
            search_list = entry.get("report", {}).get("results", {}).get("search", [])
            if isinstance(search_list, dict): search_list = [search_list]
            for search in search_list:
                if not isinstance(search, dict): continue
                title = search.get("query_title", "")
                if not title or title in processed: continue
                processed.add(title)
                idx_counter += 1
                qlen = int(search.get("query_len", 0))
                hits = search.get("hits", [])

                if not hits or qlen == 0:
                    counts["Uncharacterised"] += 1
                    if idx_counter = 70.0:
                        if p_id >= 70.0:
                            desc = hit.get("description", [{}])
                            t = desc[0].get("title", "") if isinstance(desc, list) and len(desc) > 0 else desc.get("title", "")
                            cls = classify_description(t)
                            counts[cls] += 1
                            allocated = True
                            if cls == "Uncharacterised" and idx_counter < len(all_records):
                                extracted.append(all_records[idx_counter])
                            break
                if not allocated:
                    counts["Uncharacterised"] += 1
                    if idx_counter < len(all_records): extracted.append(all_records[idx_counter])

    final_recs = extracted[:counts["Uncharacterised"]]
    SeqIO.write(final_recs, output_fasta_path, "fasta")
    print("\n -> Unique Uncharacterised Count Matrix Tally: " + str(counts["Uncharacterised"]))
    print("[SUCCESS] Successfully generated clean dunnart asset files.")

    with open(output_matrix_path, "w", newline="", encoding="utf-8") as csvfile:
        w = csv.writer(csvfile)
        w.writerow(["Species"] + all_categories + ["Total_Functional_ORs"])
        w.writerow(["Sminthopsis_crassicaudata"] + [counts[c] for c in all_categories] + [sum(counts.values())])

if __name__ == "__main__":
    main()