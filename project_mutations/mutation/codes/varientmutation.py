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
            f"ERROR at position {position}: WT has '{sequence[index]}' "
            f"but expected '{original_aa}'"
        )

    return sequence[:index] + new_aa + sequence[index+1:]


header, wt_sequence = read_fasta(input_path)

print("WT Sequence Length:", len(wt_sequence))

mutations = [

(30,'G','D'),
(33,'M','I'),
(38,'D','N'),
(39,'Q','H'),
(51,'G','R'),
(53,'N','I'),
(58,'E','A'),

(72,'N','K'),
(84,'F','I'),
(99,'M','I'),
(102,'V','L'),
(107,'S','L'),
(113,'P','L'),
(131,'L','H'),
(135,'T','I'),

(155,'K','T'),
(177,'P','S'),
(191,'F','V'),
(212,'H','Y'),
(231,'R','H'),
(232,'L','F'),
(249,'P','S'),
(258,'E','K'),

(277,'V','F'),
(283,'G','S'),
(293,'R','L'),
(315,'F','L'),
(325,'E','D'),
(334,'Y','C'),
(337,'H','Y'),
(369,'T','L'),

(393,'F','L'),
(398,'S','F'),
(410,'E','K'),
(423,'G','D'),
(433,'V','E'),
(435,'Q','L'),
(446,'Y','C'),
(460,'P','H'),

(463,'S','L'),
(484,'I','M'),
(498,'P','A'),
(510,'E','K'),
(521,'M','T'),
(526,'C','Y'),
(533,'P','S'),
(534,'S','N')

]


mutated_sequence = wt_sequence

for pos, orig, new in mutations:
    mutated_sequence = mutate_sequence(mutated_sequence, pos, orig, new)


output_file = os.path.join(output_folder, "COX2_All_Mutations.fasta")

write_fasta(
    output_file,
    ">COX2_All_Mutations",
    mutated_sequence
)

print("All mutations successfully applied.")
print("Mutated FASTA saved at:")
print(output_file)