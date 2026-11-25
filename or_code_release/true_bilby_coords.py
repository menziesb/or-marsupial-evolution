import os
import re
import csv

print('=== PINPOINTING THE 116 TRUE FUNCTIONAL BILBY GENES ===')

blast_hits_file = 'results/GCF_037893015.1_bilby.v1.9.chrom.fasta_genomic_blast_hits.txt'
bilby_clean_fasta = 'results/Bilby_OR14_True_116.fasta'
output_report = 'results/BILBY_OR14_EXACT_CHROMOSOME_LOCATIONS.txt'
output_csv = 'results/BILBY_OR14_FUNCTIONAL_COORDINATES_ONLY.csv'

if not os.path.exists(blast_hits_file) or not os.path.exists(bilby_clean_fasta):
    print('[ERROR] Missing your raw blast hits file or the clean Bilby 116 FASTA file.')
    exit(1)

# 1. Harvest the exact original row line indices from the headers of your clean 116 file
target_line_numbers = set()
with open(bilby_clean_fasta, 'r') as f:
    for line in f:
        if line.startswith('>'):
            match = re.search(r'line_(\d+)', line)
            if match:
                target_line_numbers.add(int(match.group(1)))

# If headers are sequential, fallback to tracking using your original tree description markers
if not target_line_numbers:
    if os.path.exists('results/PERAMELEMORPH_OR14_FOR_PHYLOGENY.fasta'):
        with open('results/PERAMELEMORPH_OR14_FOR_PHYLOGENY.fasta', 'r') as f:
            for line in f:
                if line.startswith('>') and ('Bilby' in line or 'Macrotis' in line):
                    match = re.search(r'gene_id_(\d+)', line)
                    if match:
                        target_line_numbers.add(int(match.group(1)))

print(f' -> Successfully locked onto your {len(target_line_numbers)} unique functional gene lines.')

functional_rows = []
unique_loci = set()

# 2. Extract coordinates ONLY for those 116 genuine functional lines
with open(blast_hits_file, 'r', encoding='utf-8') as f:
    for idx, line in enumerate(f):
        line_num = idx + 1
        if line_num not in target_line_numbers: 
            continue
        
        parts = line.strip().split('\t')
        if len(parts) < 7: 
            continue
        
        chrom = parts[1].replace('>', '').strip()
        try:
            p1 = int(parts[5])
            p2 = int(parts[6])
            pident = float(parts[2])
            
            start_fixed = min(p1, p2)
            end_fixed = max(p1, p2)
            strand = 'Reverse' if p1 > p2 else 'Forward'
            
            coord_key = (chrom, start_fixed, end_fixed)
            if coord_key not in unique_loci:
                unique_loci.add(coord_key)
                functional_rows.append({
                    'line': line_num, 'chrom': chrom, 'start': start_fixed, 
                    'end': end_fixed, 'strand': strand, 'ident': pident
                })
        except ValueError:
            continue

# Sort rows physically from start to end along the chromosome track
functional_rows.sort(key=lambda x: (x['chrom'], x['start']))

# 3. Write data grids and text directory reports
csv_data = [['Tree_ID', 'Chromosome_Scaffold', 'Start_Position_BP', 'End_Position_BP', 'Strand', 'Bait_Identity']]
with open(output_report, 'w') as out_f:
    out_f.write('====================================================================\n')
    out_f.write('===      GREATER BILBY FAMILY 14 TRUE GENOMIC DIRECTORY          ===\n')
    out_f.write('====================================================================\n\n')
    out_f.write(f'Total Unique High-Stringency Functional Genes Mapped: {len(functional_rows)}\n\n')
    
    current_chr = ''
    chr_genes = []
    
    for idx, g in enumerate(functional_rows):
        gene_label = f'Bilby_OR14_gene_{idx+1}'
        csv_data.append([gene_label, g['chrom'], g['start'], g['end'], g['strand'], g['ident']])
        
        if g['chrom'] != current_chr:
            if current_chr and chr_genes:
                c_min, c_max = min(chr_genes), max(chr_genes)
                out_f.write(f'=== Chromosome Accession: {current_chr} ===\n')
                out_f.write(f' -> Total Verified Functional Genes in Cluster: {len(chr_genes)}\n')
                out_f.write(f' -> Precise Window: {c_min:,} bp to {c_max:,} bp\n')
                out_f.write(f' -> True Functional Locus Span: {(c_max-c_min)/1000000.0:.2f} Mb\n\n')
            current_chr = g['chrom']
            chr_genes = []
        chr_genes.append(g['start'])
        
    if current_chr and chr_genes:
        c_min, c_max = min(chr_genes), max(chr_genes)
        out_f.write(f'=== Chromosome Accession: {current_chr} ===\n')
        out_f.write(f' -> Total Verified Functional Genes in Cluster: {len(chr_genes)}\n')
        out_f.write(f' -> Precise Window: {c_min:,} bp to {c_max:,} bp\n')
        out_f.write(f' -> True Functional Locus Span: {(c_max-c_min)/1000000.0:.2f} Mb\n\n')

with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
    csv.writer(csvfile).writerows(csv_data)

print(f'\n[SUCCESS] Extracted exact chromosomal addresses for all {len(functional_rows)} true genes!')
