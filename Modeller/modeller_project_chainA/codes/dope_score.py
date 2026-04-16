from modeller import *
from modeller.scripts import complete_pdb

env = Environ()

# read topology and parameter libraries
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

models = ['COX2_4.B99990001.pdb',
          'COX2_4.B99990002.pdb',
          'COX2_4.B99990003.pdb',
          'COX2_4.B99990004.pdb',
          'COX2_4.B99990005.pdb']

for m in models:
    mdl = complete_pdb(env, m)
    s = Selection(mdl)
    score = s.assess_dope()
    print(m, score)