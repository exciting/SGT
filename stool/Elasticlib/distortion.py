from ase import *

class distortion(Atoms):
    def __ini__(self,atoms,cell=None):
        self=atoms
        self.set_cell(cell)