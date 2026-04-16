import subprocess

input_fasta = r"D:\ayesha_project\ayesha_project\input\msa_inputs_hit_seq - Copy\5IKQ_msa_input2.fasta"
output_alignment = r"D:\ayesha_project\ayesha_project\output\aligned_5IKQ2.fasta"

cmd = f'mafft --auto "{input_fasta}"'

with open(output_alignment, "w") as outfile:
    subprocess.run(cmd, shell=True, stdout=outfile)

print("Alignment completed successfully.")
print("Output saved at:", output_alignment)