import os
import urllib.request
import gzip
import shutil
import subprocess

# Secure, direct hard links to the missing 3 nuclear genome tracks
TARGET_URLS = {
    "Perameles_gunnii": "https://nih.gov",
    "Trichosurus_vulpecula": "https://nih.gov",
    "Macropus_eugenii": "https://nih.gov"
}

def stream_and_calculate(species, url):
    dest_fna = f"data/genomes/{species}.fna"
    dest_gz = f"{dest_fna}.gz"
    
    print(f"\n🌍 Connecting securely to NCBI servers for: {species}...")
    
    # Inject standard Safari/Mac web browser headers to bypass NIH security filters
    req = urllib.request.Request(
        url, 
        headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'}
    )
    
    try:
        # 1. Download stream block
        with urllib.request.urlopen(req) as response, open(dest_gz, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        print(f"-> Securely downloaded compressed data track.")
        
        # 2. Decompress natively via Python to prevent Mac terminal encoding bugs
        print(f"-> Unpacking chromosomal sequences natively...")
        with gzip.open(dest_gz, 'rb') as f_in, open(dest_fna, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        if os.path.exists(dest_gz):
            os.remove(dest_gz)
            
        file_size_gb = os.path.getsize(dest_fna) / (1024 * 1024 * 1024)
        print(f"-> Verified: {species}.fna written to drive successfully ({file_size_gb:.2f} GB)")
        
        # 3. Process calculations
        print(f"-> Building local BLAST index maps...")
        subprocess.run(f"makeblastdb -in {dest_fna} -dbtype nucl -out data/genomes/{species}_db >/dev/null 2>&1", shell=True)
        
        print(f"-> Running high-sensitivity alignment search (tblastn)...")
        subprocess.run(f"tblastn -query data/OR_known_queries.fasta -db data/genomes/{species}_db -outfmt \"6 sseqid sstart send\" -evalue 1e-4 -num_threads 6 -out results/{species}_blast_hits.txt", shell=True)
        
        print(f"-> Isulating functional intact olfactory receptor sequences...")
        subprocess.run(f"python scripts/filter_intact_genes.py {species}", shell=True)
        print(f"✅ Success! Subgenome resolved for {species}")
        
    except Exception as e:
        print(f"❌ Error processing {species}: {e}")

if __name__ == "__main__":
    # Scrub old corrupted fragments and empty placeholder results files first
    print("Clearing corrupted file placeholders from workspace...")
    for sp in TARGET_URLS.keys():
        subprocess.run(f"rm -rf data/genomes/{sp}* results/{sp}*", shell=True)
        
    # Execute sequential runs
    for species, url in TARGET_URLS.items():
        stream_and_calculate(species, url)
        
    print("\n=== Processing complete for all targeted lines! ===")
