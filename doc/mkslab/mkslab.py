from sgt import slab
from ase.io import read,write
bulk=read('Rutile.struct',format='struct')


slab=slab(bulk,layers=3,miller=[1,0,0],vacuum=10,method=2)

write('slab.pov',slab.get_atoms(), display=False, run_povray=True)
write('slab3.xsf',slab.get_atoms())
slab.set_layers(4)
write('slab4.xsf',slab.get_atoms())
write('bulk.xsf',bulk)
