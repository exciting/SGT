from stool import ElasticDistortion,loadObject
from ase.calculators import Exciting
from ase import read

exciting=Exciting(bin='/fshome/chm/git/exciting/bin/excitingser', kpts=(4, 4, 4), xctype='GGArevPBE')

atoms=primcell=read('AlBulk.xml',format='exi')
print atoms.get_cell()
distortions=ElasticDistortion(structure=atoms, 
                              calculator=exciting, 
                              order=3, 
                              maxstrain=0.02, 
                              distortions=5)
distortions.print_distortions()
distortions.calculate()

distortions.save("distwenergy")
