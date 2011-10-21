from ase import *
import copy
class distortion():
    def __init__(self,atoms,calculator=None,cell=None, eta=None,LagrangeS= None):
       # super(distortion, self).__init__()
        self.atoms=atoms.copy()
        self.atoms.set_cell(cell)
        self.atoms.set_calculator(copy.deepcopy(calculator))
        self.eta=eta
        self.LagrangeS=LagrangeS
        

    def print_distortion(self):
        print "LagrangeS:", self.LagrangeS
        print "eta:      ",self.eta
        print self.atoms.get_cell()