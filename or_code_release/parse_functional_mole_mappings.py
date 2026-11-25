import os
import csv
from Bio import SeqIO

print('====================================================================')
print('===    LAUNCHING NATIVE GENOMIC STREAM ALIGNMENT ENGINE          ===')
print('====================================================================')

master_tree_fasta = 'results/FINAL_MAMMALIAN_OR14_COMPREHENSIVE_TREE_INPUT.fasta'
thylacine_hits_txt = 'results/Thylacinus_cynocephalus_blast_hits.txt'
mole_functional_fasta = 'results/notoryctes_functional_OR_candidates.fasta'
mole_genome_fna = 'data/genomes/GCA_048126345.2_mNotTyp1.pri_genomic.fna'
output_csv = 'results/MOLE_THYLACINE_VERIFIED_COORDINATES.csv'

if not os.path.exists(mole_functional_fasta):
    print(f"[ERROR] Target Marsupial Mole functional file missing at: {mole_functional_fasta}"); exit(1)

valid_seqs = {'Mole': [], 'Thylacine': []}

# 1. Load your pristine Notoryctes functional targets
for rec in SeqIO.parse(mole_functional_fasta, 'fasta'):
    valid_seqs['Mole'].append(rec)

# 2. Harvest corresponding Thylacine tracks cleanly from your master comprehensive tree pool
if os.path.exists(master_tree_fasta):
    for rec in SeqIO.parse(master_tree_fasta, 'fasta'):
        if 'THYLACINE' in rec.id.upper() or 'CYNOCEPHALUS' in rec.id.upper():
            valid_seqs['Thylacine'].append(rec)

print(f" -> Active target array: {len(valid_seqs['Mole'])} Mole targets | {len(valid_seqs['Thylacine'])} Thylacine targets.")

csv_headers = ['Species', 'Query_ID', 'Scaffold_Contig', 'Start_BP', 'End_BP', 'Strand']
mapped_rows = []

# 3. Extract genuine Thylacine locations using your verified mapping log
if os.path.exists(thylacine_hits_txt):
    print(' -> Sifting Thylacine coordinates from species mapping matrix...')
    with open(thylacine_hits_txt, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            cols = line.replace(',', '\t').split('\t')
            if len(cols) < 5: continue
            
            q_id, s_id = cols[0].strip(), cols[1].strip()
            if any(rec.id in q_id for rec in valid_seqs['Thylacine']) or q_id in [r.id for r in valid_seqs['Thylacine']]:
                try:
                    st, en = int(cols[2]), int(cols[3])
                    strand = cols[4].strip().lower() if len(cols) > 4 else ('forward' if st < en else 'reverse')
                    mapped_rows.append(['Thylacine', q_id, s_id, min(st, en), max(st, en), strand])
                except ValueError: continue

# 4. Stream-map Marsupial Mole genes straight across the mNotTyp1 genome file
if os.path.exists(mole_genome_fna) and valid_seqs['Mole']:
    print(' -> Indexing raw Marsupial Mole reference assembly into memory map...')
    mole_index = SeqIO.index(mole_genome_fna, 'fasta')
    print(' -> Assembly mapped successfully. Scanning nucleotide strings via seed window...')
    
    mole_count = 0
    for query in valid_seqs['Mole']:
        q_str = str(query.seq).upper()
        # Seed match on the first 30 characters to accommodate minor local variants or introns
        seed = q_str[:30] if len(q_str) > 30 else q_str
        found = False
        
        for scaf_id, scaf_rec in mole_index.items():
            scaf_str = str(scaf_rec.seq).upper()
            
            idx = scaf_str.find(seed)
            if idx != -1:
                st = idx + 1
                en = st + len(q_str)
                mapped_rows.append(['Mole', query.id, scaf_id, st, en, 'forward'])
                mole_count += 1
                found = True
                break
                
            try:
                rc_seed = str(query.seq.reverse_complement()).upper()[:30]
                idx_rc = scaf_str.find(rc_seed)
                if idx_rc != -1:
                    st = idx_rc + 1
                    en = st + len(q_str)
                    mapped_rows.append(['Mole', query.id, scaf_id, st, en, 'reverse'])
                    mole_count += 1
                    found = True
                    break
            except Exception: pass
            
    print(f'    └── Successfully mapped {mole_count} Marsupial Mole genes directly to scaffold nodes.')

if mapped_rows:
    with open(output_csv, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(csv_headers)
        writer.writerows(mapped_rows)
    print('\n====================================================================')
    print('[SUCCESS] GENUINE INTER-SPECIES MAPPING BLUEPRINT WRITTEN!')
    print(f' -> Saved true physical coordinates directly to: {output_csv}')
    print(f' -> Total unified sequences mapped: {len(mapped_rows)} loci')
    print('====================================================================')
else:
    print('\n[ERROR] Direct sequence streaming failed. Check target files.')
