from StackinfFault import  StackingFault
from ase import read, write

primcell=read('00.struct',format='struct')

supercell=StackingFault(primcell,n_slab=3,n_block=1,n_shift=1,shift=1)
 
write('sc.pov',supercell, display=False, run_povray=True)

write("supercell.xml",supercell, show_unit_cell=2,format="exi")