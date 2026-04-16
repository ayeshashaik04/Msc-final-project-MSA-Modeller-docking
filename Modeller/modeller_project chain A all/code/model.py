from modeller import *
from modeller.automodel import *

env = environ()

a = automodel(
    env,
    alnfile  = 'alignment.ali',            
    knowns   = '5IKQ',                     
    sequence = 'COX2_All_Mutations'        
)

a.starting_model = 1
a.ending_model   = 5  # generate 5 models

# Build the models
a.make()