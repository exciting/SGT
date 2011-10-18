import RandomSQS
from ase.io import read, write

Bulk=read('Tungsten135ideal',format='vasp')

repeat=4
concentration=0.2
suplement="Ir"

WIrSQS=RandomSQS(Bulk,repeat,concentration,suplement)

write('WIr.pov', WIrSQS, show_unit_cell=2, display=False, run_povray=True)