from Bio import AlignIO
import pandas as pd
from collections import Counter

input_path = r"D:\ayesha_project\ayesha_project\output\aligned_5IKQ2.fasta"
output_path = r"D:\ayesha_project\ayesha_project\output\complete_mutation_analysis2.csv"

alignment = AlignIO.read(input_path, "fasta")

num_sequences = len(alignment)
alignment_length = alignment.get_alignment_length()

sequence_ids = [record.id for record in alignment]
reference = alignment[0].seq  # First sequence as reference

results = []

for i in range(alignment_length):

    row = {"Position": i + 1}
    column = alignment[:, i]
    ref_residue = reference[i]

    # Store amino acids for all sequences
    for seq_index in range(num_sequences):
        row[sequence_ids[seq_index]] = alignment[seq_index].seq[i]

    # Count amino acids (excluding gaps for conservation calculation)
    counts = Counter(column)
    count_string = "; ".join([f"{aa}={count}" for aa, count in counts.items()])
    row["AA_Counts"] = count_string

    if "-" in counts:
        non_gap_total = num_sequences - counts["-"]
    else:
        non_gap_total = num_sequences

    if non_gap_total > 0:
        most_common_residue, count = counts.most_common(1)[0]
        conservation_percentage = (count / num_sequences) * 100
    else:
        conservation_percentage = 0

    row["Conservation_Percentage"] = f"{round(conservation_percentage, 2)}%"

    mutation_types = set()

    for seq_index in range(1, num_sequences):
        residue = alignment[seq_index].seq[i]

        if residue != ref_residue:

            if ref_residue == "-" and residue != "-":
                mutation_types.add("Insertion")

            elif ref_residue != "-" and residue == "-":
                mutation_types.add("Deletion")

            elif residue == "*":
                mutation_types.add("Nonsense")

            else:
                mutation_types.add("Missense")
                mutation_types.add("Substitution")
                mutation_types.add("Point mutation")

    if mutation_types:
        row["Mutation_Status"] = "Mutation"
        row["Mutation_Type"] = ", ".join(sorted(mutation_types))
    else:
        row["Mutation_Status"] = "No Mutation"
        row["Mutation_Type"] = "None"

    results.append(row)

df = pd.DataFrame(results)
df.to_csv(output_path, index=False)

print("Complete mutation analysis generated successfully.")