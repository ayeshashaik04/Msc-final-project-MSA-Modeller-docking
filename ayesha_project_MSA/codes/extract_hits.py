# extract sequences from BLAST XML and prepare FASTA for MSA

import os
from Bio import SeqIO
from Bio.Blast import NCBIXML

# --- PATH SETUP ---
base_dir = "/media/vanajakshi/PENDRIVE/ayesha_project"
genomes_path = os.path.join(base_dir, "/media/vanajakshi/PENDRIVE/ayesha_project/output/all_genomes.faa")  # master genome sequences
xml_folder = os.path.join(base_dir, "/media/vanajakshi/PENDRIVE/ayesha_project/output/blast_results")      # folder with all BLAST XMLs
template_folder = "/media/vanajakshi/PENDRIVE/ayesha_project/input/rcsb_pdb_5IKQ.fasta"   # folder with template FASTAs
output_folder = os.path.join(base_dir, "/media/vanajakshi/PENDRIVE/ayesha_project/output/msa_inputs_hit_seq")
os.makedirs(output_folder, exist_ok=True)

# --- PARAMETERS ---
identity_cutoff = 90
evalue_cutoff = 1e-5
min_align_len = 70
template_ext = ".fasta"  # your template FASTA extension

# --- Load all genome sequences once ---
print("🔹 Loading all genome sequences...")
genome_records = list(SeqIO.parse(genomes_path, "fasta"))
print(f"✅ Loaded {len(genome_records)} genome protein sequences.")

# --- Function to extract hits from BLAST XML ---
def extract_hits_from_xml(xml_path):
    hit_descriptions = set()
    with open(xml_path) as f:
        blast_records = NCBIXML.parse(f)
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    identity = (hsp.identities / hsp.align_length) * 100
                    if hsp.expect <= evalue_cutoff and identity >= identity_cutoff and hsp.align_length >= min_align_len:
                        hit_descriptions.add(alignment.hit_def.strip())
    return hit_descriptions

# --- MAIN LOOP ---
for template_file in os.listdir(template_folder):
    if not template_file.endswith(template_ext):
        continue

    template_name = os.path.splitext(template_file)[0]
    
    # Find corresponding XML file
    xml_file = None
    for f in os.listdir(xml_folder):
        if f.endswith(".xml") and template_name in f:
            xml_file = os.path.join(xml_folder, f)
            break

    if xml_file is None:
        print(f"⚠️ Skipping {template_name} — XML file not found.")
        continue

    print(f"\n🔹 Processing template: {template_name}")

    # 1️⃣ Extract hit descriptions from XML
    hits = extract_hits_from_xml(xml_file)
    print(f"   Found {len(hits)} BLAST hits passing filters.")

    # 2️⃣ Match full genome sequences using description substring
    matched_records = []
    for record in genome_records:
        for hit_desc in hits:
            if hit_desc in record.description:
                matched_records.append(record)
                break
    print(f"   Retrieved {len(matched_records)} full sequences from genome.")

    # 3️⃣ Read template sequence(s)
    template_path = os.path.join(template_folder, template_file)
    template_records = list(SeqIO.parse(template_path, "fasta"))

    # 4️⃣ Combine template + hits
    combined_records = template_records + matched_records

    # 5️⃣ Save combined FASTA
    output_path = os.path.join(output_folder, f"{template_name}_msa_input.fasta")
    SeqIO.write(combined_records, output_path, "fasta")
    print(f"✅ Saved combined FASTA: {output_path}")

print("\n🎯 All templates processed successfully!")