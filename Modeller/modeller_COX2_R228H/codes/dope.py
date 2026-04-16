from modeller import *
from modeller.scripts import complete_pdb

env = environ()

# Load topology and parameter libraries
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

models = ['COX2_R228H.B99990001.pdb',
          'COX2_R228H.B99990002.pdb',
          'COX2_R228H.B99990003.pdb',
          'COX2_R228H.B99990004.pdb',
          'COX2_R228H.B99990005.pdb']

for m in models:
    mdl = complete_pdb(env, m)
    s = selection(mdl)
    dope_score = s.assess_dope()
    
    print("Model:", m)
    print("DOPE score:", dope_score)
    print("-"*40)