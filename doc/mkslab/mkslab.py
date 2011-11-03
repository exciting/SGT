from stool.mkslab import slab
from ase import read
bulk=read('../Elastic/AlBulk.xml',format='exi')
slab=slab(bulk,layers=3,miller=[2,4,4])
print slab.get_miller()