import os
import csv

print('====================================================================')
print('===    LAUNCHING PRECISION ID-MATCHED THYLACINE MAPPER           ===')
print('====================================================================')

true_protein_faa = 'data/genomes/Thylacinus_cynocephalus_translated.faa'
raw_hits_txt = 'results/Thylacinus_cynocephalus_blast_hits.txt'
output_csv = 'results/THYLACINE_TRUE_COORDINATES.csv'

# Fallback path check if the translated file sits in results
if not os.path.exists(true_protein_faa) and os.path.exists('results/Thylacinus_cynocephalus_translated.faa'):
    true_protein_faa = 'results/Thylacinus_cynocephalus_translated.faa'

if not os.path.exists(true_protein_faa):
    print(f"[ERROR] Missing true protein reference file at: {true_protein_faa}"); exit(1)
if not os.path.exists(raw_hits_txt):
    print(f"[ERROR] Missing Thylacine raw hits log at: {raw_hits_txt}"); exit(1)

# 1. Extract the exact tracking keys from your verified protein file headers
# If headers use formatting like >Query_X_CM040574.1_start_end, we extract the numeric footprint
true_coordinate_keys = set()
from Bio import SeqIO

for rec in SeqIO.parse(true_protein_faa, 'fasta'):
    # Capture any numeric chunks in the header text that map to physical positions
    import re
    nums = re.findall(r'\d+', rec.description if rec.description else rec.id)
    for n in nums:
        if len(n) > 5:  # Lock onto large base-pair coordinate integers
            true_coordinate_keys.add(n)

print(f" -> Successfully harvested {len(true_coordinate_keys)} unique structural coordinates from your true protein file.")

csv_headers = ['Query_Seed', 'Scaffold', 'Identity', 'Length', 'E_Value', 'Start_BP', 'End_BP', 'Strand']
mapped_rows = []
seen_starts = set()

# 2. Parse the raw hit rows, matching strictly against your true tree protein keys
with open(raw_hits_txt, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        cols = [c.strip() for c in line.split() if c.strip()]
        if len(cols) < 7: continue
        
        start_str = cols[5]
        end_str = cols[6]
        
        # Cross-reference step: Check if this row's position matches your true protein file digits
        if (start_str in true_coordinate_keys) or (end_str in true_coordinate_keys):
            try:
                start_bp = int(start_str)
                end_bp = int(end_str)
                
                coord_lock = min(start_bp, end_bp)
                if coord_lock in seen_starts: continue
                
                strand = 'forward' if start_bp < end_bp else 'reverse'
                mapped_rows.append([
                    cols[0], cols[1], cols[2], cols[3], cols[4],
                    min(start_bp, end_bp), max(start_bp, end_bp), strand
                ])
                seen_starts.add(coord_lock)
            except ValueError: continue

# Sort chronologically by chromosome position
mapped_rows.sort(key=lambda x: x[5])

with open(output_csv, 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(csv_headers)
    writer.writerows(mapped_rows)

print('\n====================================================================')
print('[SUCCESS] THYLACINE COORDINATES ISOLATED VIA RIGID ID-MATCHING!')
print(f' -> Saved true physical coordinates directly to: {output_csv}')
print(f' -> Total un-scrambled, true sequences mapped: {len(mapped_rows)} genes')
print('====================================================================')
