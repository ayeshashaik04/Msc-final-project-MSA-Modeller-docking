import os

input_path = r"D:\swissmodal project\mutation\inputs\COX2_WT.fasta"
output_folder = r"D:\swissmodal project\mutation\outputs"

os.makedirs(output_folder, exist_ok=True)


def read_fasta(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    header = lines[0].strip()
    sequence = "".join(line.strip() for line in lines[1:])
    return header, sequence


def write_fasta(filename, header, sequence):
    with open(filename, "w") as f:
        f.write(header + "\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")


def mutate_sequence(sequence, position, original_aa, new_aa):
    index = position - 1  # convert to 0-based index
    
    if sequence[index] != original_aa:
        raise ValueError(
            f"ERROR: Position {position} contains '{sequence[index]}', "
            f"expected '{original_aa}'. Check your WT sequence numbering."
        )
    
    return sequence[:index] + new_aa + sequence[index+1:]

header, wt_sequence = read_fasta(input_path)

print("WT length:", len(wt_sequence))
print("228:", wt_sequence[227])
print("428:", wt_sequence[427])
print("488:", wt_sequence[487])
print("511:", wt_sequence[510])


mut_228 = mutate_sequence(wt_sequence, 228, "R", "H")
write_fasta(os.path.join(output_folder, "COX2_R228H.fasta"),
            ">COX2_R228H", mut_228)

mut_428 = mutate_sequence(wt_sequence, 428, "P", "A")
write_fasta(os.path.join(output_folder, "COX2_P428A.fasta"),
            ">COX2_P428A", mut_428)

mut_488 = mutate_sequence(wt_sequence, 488, "E", "G")
write_fasta(os.path.join(output_folder, "COX2_E488G.fasta"),
            ">COX2_E488G", mut_488)

mut_511 = mutate_sequence(wt_sequence, 511, "V", "A")
write_fasta(os.path.join(output_folder, "COX2_V511A.fasta"),
            ">COX2_V511A", mut_511)

quad = wt_sequence
quad = mutate_sequence(quad, 228, "R", "H")
quad = mutate_sequence(quad, 428, "P", "A")
quad = mutate_sequence(quad, 488, "E", "G")
quad = mutate_sequence(quad, 511, "V", "A")

write_fasta(os.path.join(output_folder,
            "COX2_Quad_R228H_P428A_E488G_V511A.fasta"),
            ">COX2_Quad_R228H_P428A_E488G_V511A", quad)

print("All 4 single mutants and 1 quadruple mutant created successfully.")