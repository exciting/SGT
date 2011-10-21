#!/usr/bin/env python
# -*- coding: utf8 -*-

from stool import StackingFaultShift
from ase import read, write

primcell=read('00.struct',format='struct')


supercell=StackingFaultShift(primcell,n_slab=3,n_block=2,n_shift=3,shift=1)

write("supercell.xml",supercell, show_unit_cell=2,format="exi") 
write('sc.pov',supercell, display=False, run_povray=True)

write("supercell.xsf",supercell,format="xsf")
write("supercell.struct",supercell,format="struct")
