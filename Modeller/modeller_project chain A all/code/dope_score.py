from modeller import *
from modeller.scripts import complete_pdb

env = Environ()

env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

models = ['COX2_All_Mutations.B99990001.pdb',
          'COX2_All_Mutations.B99990002.pdb',
          'COX2_All_Mutations.B99990003.pdb',
          'COX2_All_Mutations.B99990004.pdb',
          'COX2_All_Mutations.B99990005.pdb']

for m in models:
    mdl = complete_pdb(env, m)
    s = Selection(mdl)
    dope_score = s.assess_dope()
    print("Model:", m, "DOPE score:", dope_score)