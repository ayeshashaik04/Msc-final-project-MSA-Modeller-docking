from Bio.PDB import PDBParser, Superimposer, PDBIO
import os

template_file = r"D:\RMSD_MODELLER\PDB Structures_modeller\Inputs\5IKQ_A.pdb"
model_file = r"D:\RMSD_MODELLER\PDB Structures_modeller\Inputs\COX2_4.B99990004.pdb"

output_folder = r"D:\RMSD_MODELLER\PDB Structures_modeller\Outputs"

aligned_output = os.path.join(output_folder, "model_aligned.pdb")
rmsd_output = os.path.join(output_folder, "RMSD_result.txt")

parser = PDBParser(QUIET=True)

template_structure = parser.get_structure("TEMPLATE", template_file)
model_structure = parser.get_structure("MODEL", model_file)

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

io = PDBIO()
io.set_structure(model_structure)
io.save(aligned_output)

with open(rmsd_output, "w") as f:
    f.write(f"RMSD between template and model: {rmsd_value:.3f} Å")

print("RMSD Calculation Completed")
print("RMSD =", round(rmsd_value,3), "Å")
print("Aligned structure saved at:", aligned_output)
print("RMSD result saved at:", rmsd_output)