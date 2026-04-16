from Bio import AlignIO
import pandas as pd
from collections import Counter
import os

input_path = r"D:\ayesha_project_MSA\ayesha_project_MSA\output\aligned_5IKQ2.fasta"
output_path = r"D:\ayesha_project_MSA\ayesha_project_MSA\output\mutation_summary_5IKQ.csv"

os.makedirs(os.path.dirname(output_path), exist_ok=True)

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

    for seq_index in range(num_sequences):
        row[sequence_ids[seq_index]] = alignment[seq_index].seq[i]

    counts = Counter(column)
    count_string = "; ".join([f"{aa}={count}" for aa, count in counts.items()])
    row["AA_Counts"] = count_string

    non_gap_counts = {aa: count for aa, count in counts.items() if aa != "-"}
    non_gap_total = sum(non_gap_counts.values())

    if non_gap_total > 0:
        most_common_residue = max(non_gap_counts, key=non_gap_counts.get)
        conservation_percentage = (non_gap_counts[most_common_residue] / num_sequences) * 100
    else:
        conservation_percentage = 0

    row["Conservation_Percentage"] = f"{round(conservation_percentage, 2)}%"

    mutation_details = []

    insertion_count = 0
    deletion_count = 0
    missense_count = 0
    nonsense_count = 0

    for seq_index in range(1, num_sequences):
        residue = alignment[seq_index].seq[i]
        seq_id = sequence_ids[seq_index]

        if residue != ref_residue:

            if ref_residue == "-" and residue != "-":
                mutation_type = "Insertion"
                insertion_count += 1

            elif ref_residue != "-" and residue == "-":
                mutation_type = "Deletion"
                deletion_count += 1

            elif residue == "*":
                mutation_type = "Nonsense"
                nonsense_count += 1

            else:
                mutation_type = "Missense (Substitution/Point)"
                missense_count += 1

            mutation_details.append(f"{seq_id}:{mutation_type}")

    if mutation_details:
        row["Mutation_Status"] = "Mutation"
        row["Mutation_Details"] = " | ".join(mutation_details)
    else:
        row["Mutation_Status"] = "No Mutation"
        row["Mutation_Details"] = "None"

    # Optional summary counts (VERY useful for analysis)
    row["Insertion_Count"] = insertion_count
    row["Deletion_Count"] = deletion_count
    row["Missense_Count"] = missense_count
    row["Nonsense_Count"] = nonsense_count

    results.append(row)

df = pd.DataFrame(results)
df.to_csv(output_path, index=False)

print(" Complete mutation analysis generated successfully.")
print(" Output saved at:", output_path)