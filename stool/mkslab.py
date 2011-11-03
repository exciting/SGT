from ase import *
from mkslablib import  slab_ducvec
from mkslablib.norm_miller import norm_miller
from mkslablib.fill_with_atoms import fill_with_atoms
from mkslablib.add_vacuum import add_vacuum
from numpy import linalg
import numpy as np
class  slab(Atoms):
    def __init__(self,structure,miller=[1,1,1],method=2,layers=0,vacuum=0.0):
        """ create slab from bulk structure
            arguments:
            structure: Atoms
                Bulk structure to create slab from
            method: integer
                (1) z-vector is perpendicular to surface [DEFAULT]
                (2) z-vector is perpendicular to surface(with depth periodicity)
                (3) z-vector is the nearest to normal in Nth layer
                
        
        """
        super(slab, self).__init__(structure)
        self.miller=norm_miller(miller)
        self.method=method
        self.layers=layers
        self.ducvec=self.get_cell()
        self.rucvec=np.transpose(linalg.inv(self.ducvec))
        h,k,l=self.miller
        self.Ghkl=h*self.rucvec[0][:]+k*self.rucvec[1][:]+l*self.rucvec[2][:]
        if self.layers==0:
            print "zero layers does not define a slab"
        
        else:
            slab_ducvec(self)
            fill_with_atoms(self)
            add_vacuum(self,vacuum)
            
    def get_miller(self):
        return self.miller
    def get_layers(self):
        return self.layers
    def get_ducvec(self):
        return self.ducvec