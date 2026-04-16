from modeller import *
from modeller.automodel import *

env = environ()

a = automodel(env,
              alnfile='alignment.ali',
              knowns='5IKQ',
              sequence='COX2_R228H')   

a.starting_model = 1
a.ending_model   = 5

a.assess_methods = (assess.DOPE,)

a.make()