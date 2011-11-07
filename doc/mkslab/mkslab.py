from stool.mkslab import slab
from ase import read,write
bulk=read('../Elastic/AlBulk.xml',format='exi')
slab=slab(bulk,layers=3,miller=[1,-2,-2])
 
write('slab.pov',slab, display=False, run_povray=True)