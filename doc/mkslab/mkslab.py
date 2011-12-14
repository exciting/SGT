from sgt.mkslab import slab
from ase import read,write
bulk=read('Rutile.struct',format='struct')


slab=slab(bulk,layers=3,miller=[1,0,0],vacuum=10,method=2)

write('slab.pov',slab.get_atoms(), display=False, run_povray=True)
write('slab2.xsf',slab.get_atoms())
write('bulk.xsf',bulk)