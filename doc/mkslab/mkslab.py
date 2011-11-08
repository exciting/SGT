from sgt.mkslab import slab
from ase import read,write
bulk=read('AlBulk.xml',format='exi')

slab=slab(bulk,layers=3,miller=[1,1,0],vacuum=10)
 
write('slab.pov',slab, display=False, run_povray=True)
write('slab.xsf',slab)
write('bulk.xsf',bulk)