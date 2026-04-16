from Bio.PDB import PDBParser, Superimposer
import os

template_file = r"D:\RMSD_MODELLER\PDB Structures\Inputs\5IKQ_B.pdb"
model_file = r"D:\RMSD_MODELLER\PDB Structures\Inputs\COX2_All_Mutations.B99990002.pdb"

parser = PDBParser(QUIET=True)

template = parser.get_structure("template", template_file)
model = parser.get_structure("model", model_file)

template_res = {}
model_res = {}

for res in template.get_residues():
    res_id = res.get_id()[1]
    template_res[res_id] = res

for res in model.get_residues():
    res_id = res.get_id()[1]
    model_res[res_id] = res

template_atoms = []
model_atoms = []

for res_id in template_res:
    if res_id in model_res:
        t_res = template_res[res_id]
        m_res = model_res[res_id]

        if t_res.has_id("CA") and m_res.has_id("CA"):
            template_atoms.append(t_res["CA"])
            model_atoms.append(m_res["CA"])


sup = Superimposer()
sup.set_atoms(template_atoms, model_atoms)

rmsd_value = sup.rms

print("RMSD =", round(rmsd_value,3), "Å")
print("Total matched residues:", len(template_atoms))