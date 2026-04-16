from Bio.PDB import PDBParser, Superimposer, PDBIO
import os

template_file = r"D:\RMSD_MODELLER\PDB Structures_modeller\Inputs\5IKQ_A.pdb"

model_files = [
    r"D:\RMSD_MODELLER\PDB Structures_modeller\Inputs\COX2_E488G.B99990001.pdb",
    r"D:\RMSD_MODELLER\PDB Structures_modeller\Inputs\COX2_P428A.B99990003.pdb",
    r"D:\RMSD_MODELLER\PDB Structures_modeller\Inputs\COX2_R228H.B99990001.pdb",
    r"D:\RMSD_MODELLER\PDB Structures_modeller\Inputs\COX2_V511A.B99990003.pdb"
]

output_folder = r"D:\RMSD_MODELLER\PDB Structures_modeller\Outputs"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

rmsd_output = os.path.join(output_folder, "RMSD_results.txt")

parser = PDBParser(QUIET=True)
template_structure = parser.get_structure("TEMPLATE", template_file)

results = []

for model_file in model_files:

    model_name = os.path.basename(model_file).split(".")[0]

    model_structure = parser.get_structure(model_name, model_file)

    template_atoms = []
    model_atoms = []

    for temp_res, model_res in zip(template_structure.get_residues(),
                                   model_structure.get_residues()):
        if temp_res.has_id("CA") and model_res.has_id("CA"):
            template_atoms.append(temp_res["CA"])
            model_atoms.append(model_res["CA"])

    super_imposer = Superimposer()
    super_imposer.set_atoms(template_atoms, model_atoms)
    super_imposer.apply(model_structure.get_atoms())

    rmsd_value = super_imposer.rms

    aligned_file = os.path.join(output_folder, model_name + "_aligned.pdb")
    io = PDBIO()
    io.set_structure(model_structure)
    io.save(aligned_file)

    results.append(f"{model_name} : RMSD = {rmsd_value:.3f} Å")

    print(f"{model_name} RMSD =", round(rmsd_value,3), "Å")


with open(rmsd_output, "w") as f:
    for line in results:
        f.write(line + "\n")

print("\nAll RMSD calculations completed")
print("Results saved at:", rmsd_output)
