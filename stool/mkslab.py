from ase import *
from mkslablib import  slab_ducvec
from mkslablib.norm_miller import norm_miller
from mkslablib.fill_with_atoms import fill_with_atoms
from mkslablib.add_vacuum import add_vacuum
from numpy import linalg
import numpy as np 
class  slab(Atoms):
    """ create slab from bulk structure
    
            arguments:
            
            structure: Atoms
                Bulk structure to create slab from
            miller: 3 dimensional integer
                Miller index describing the desired slab
            method: integer
                (1) z-vector is perpendicular to surface [DEFAULT]
                (2) z-vector is perpendicular to surface(with depth periodicity)
                (3) z-vector is the nearest to normal in Nth layer
            layers: integer
                number of layers of constructed slab
            vacuum: float
                vacuum to be added along z direction
                
        
    """
    def __init__(self,structure,miller=[1,1,1],method=2,layers=0,vacuum=0.0):
        
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
            self.slab_vec=slab_ducvec(self)
            fill_with_atoms(self)
            self.repeat(50,50,50)
            self.set_cell(self.slab_vec)
            self.remove_surplus_atoms()
            self.add_vacuum(vacuum)
            
    def get_miller(self):
        """ return Miller indices """
        return self.miller
    def get_layers(self):
        """return number of layers """
        return self.layers
    def get_ducvec(self):
        return self.ducvec
    def remove_surplus_atoms(self):
        """ remove all atoms that are not in the unit cell"""
        positions=self.get_scaled_positions()
        def not_in_box(position):
            position=atom.position
            not_in_box=False
            for i in range (0,2):
                if (position[i]<0 or position[i]>=1):
                    not_in_box=True
            return not_in_box
        del self[[atom.index for position in positions if not_in_box(position)]]
             
    def add_vacuum(self,vacuum):
        """ add vacuum along the 3rd basevector """
        cell=self.get_cell()
        cell[2]+=vacuum
        self.set_cell(cell)