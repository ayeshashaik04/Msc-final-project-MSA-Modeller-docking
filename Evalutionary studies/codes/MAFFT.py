import subprocess

input_fasta = r"D:\Evalutionary studies\input\MAGA.fasta"
output_alignment = r"D:\Evalutionary studies\output\aligned.fasta"

cmd = f'mafft --auto "{input_fasta}"'

with open(output_alignment, "w") as outfile:
    subprocess.run(cmd, shell=True, stdout=outfile)

print("Alignment completed successfully.")
print("Output saved at:", output_alignment)