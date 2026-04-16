from Bio.PDB import PDBParser, Superimposer
import os

template_file = r"D:\RMSD_MODELLER\PDB structures_swissmodel\inputs\5IKQ.pdb"

model_files = [
    r"D:\RMSD_MODELLER\PDB structures_swissmodel\inputs\COX2_4mutations.pdb",
    r"D:\RMSD_MODELLER\PDB structures_swissmodel\inputs\COX2_All mutations.pdb",
    r"D:\RMSD_MODELLER\PDB structures_swissmodel\inputs\COX2_E488G.pdb",
    r"D:\RMSD_MODELLER\PDB structures_swissmodel\inputs\COX2_P428A.pdb",
    r"D:\RMSD_MODELLER\PDB structures_swissmodel\inputs\COX2_R228H.pdb",
    r"D:\RMSD_MODELLER\PDB structures_swissmodel\inputs\COX2_V511A.pdb"
]

output_file = r"D:\RMSD_MODELLER\PDB structures_swissmodel\outputs\RMSD_chainA.txt"

parser = PDBParser(QUIET=True)
template = parser.get_structure("TEMPLATE", template_file)[0]
template_chain = template["A"]

results = []

for model_file in model_files:
    model_name = os.path.basename(model_file).replace(".pdb", "").replace(" ", "_")
    print(f"\nProcessing: {model_name}")

    try:
        model = parser.get_structure(model_name, model_file)[0]
        model_chain = model["A"]

        template_res = {res.id[1]: res for res in template_chain if res.has_id("CA")}
        model_res = {res.id[1]: res for res in model_chain if res.has_id("CA")}

        common_ids = sorted(set(template_res.keys()) & set(model_res.keys()))

        template_atoms = [template_res[i]["CA"] for i in common_ids]
        model_atoms = [model_res[i]["CA"] for i in common_ids]

        print(f"Matched residues: {len(template_atoms)}")

        if len(template_atoms) < 20:
            print("Skipping due to low matches")
            continue

        sup = Superimposer()
        sup.set_atoms(template_atoms, model_atoms)
        sup.apply(model_chain.get_atoms())

        rmsd = sup.rms
        print(f"{model_name} RMSD = {rmsd:.3f} Å")
        results.append(f"{model_name} : RMSD = {rmsd:.3f} Å")

    except Exception as e:
        print(f"Error processing {model_name}: {e}")

with open(output_file, "w") as f:
    for line in results:
        f.write(line + "\n")

print("\nDone")
print("Saved at:", output_file)