from sgt.mkslab import slab
from ase import read,write
bulk=read('AlBulk.xml',format='exi')

slab1=slab(bulk,layers=3,miller=[1,1,0],vacuum=10,method=1)
slab2=slab(bulk,layers=3,miller=[1,1,0],vacuum=10,method=2)
slab3=slab(bulk,layers=3,miller=[1,1,0],vacuum=10,method=3)
write('slab.pov',slab2, display=False, run_povray=True)
write('slab2.xsf',slab2)
write('bulk.xsf',bulk)