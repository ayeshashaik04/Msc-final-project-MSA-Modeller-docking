import subprocess

# Input FASTA file (full path)
input_fasta = r"D:\ayesha_project\ayesha_project\input\msa_inputs_hit_seq\5IKQ_msa_input.fasta"

# Output alignment file
output_alignment = r"D:\ayesha_project\ayesha_project\output\aligned_5IKQ.fasta"

# MAFFT command
cmd = f'mafft --auto "{input_fasta}"'

# Run alignment
with open(output_alignment, "w") as outfile:
    subprocess.run(cmd, shell=True, stdout=outfile)

print("Alignment completed successfully.")
print("Output saved at:", output_alignment)
