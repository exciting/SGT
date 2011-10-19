from stool import ElasticDistortion
from ase.calculators import Exciting
from ase import read

exciting=Exciting(bin='excitingser', kpts=(4, 4, 4), xctype='GGArevPBE')

atoms=primcell=read('../StackingFault/00.struct',format='struct')

distortions=ElasticDistortion(structure=atoms, calculator=exciting, order=2, maxstrain=0.03, distortions=11)

distortions.calculate()

distortions.get_elastic_c()

distortions.make_report()