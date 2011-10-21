
from ase import *

from Elasticlib  import Elastic_setup,  elastic_select_data
 
import pickle

class ElasticDistortion(Elastic_setup):
    """
    This class can setup the calculations for determining the elastic constants
    
    structure: Atoms object
    
    calculator: calculator object
    
    order: integer
        order of elastic constatns
    maxstrain:number
        strain the structures maximal by this fraction
    distortions:integer
        how many distortion per distortion kid should be generated
    
    """
    def __init__(self,structure=None,calculator=None, order=2, maxstrain=0.03,distortions=11):
       
       
       
        self.setup(structure,calculator, order, maxstrain,distortions)
    
    
    def calculate(self):
        """
        Start The calculations
        """
        for dist in self.distortions:
            if  self.mthd=="Energy":
                dist.atoms.get_potential_energy()
                print "did", dist 
            else:
                dist.atoms.get_stress()
 
    def select_data(self):
        """
        select range an fitorder to eliminate errors
        """
        elastic_select_data(self)
    def get_elastic_c(self):
        """
        calculate the elastic constants
        """
        return 0
    def make_report(self):
        """
        Write Report
        """
        print    'Order of elastic constants      =', self.order         ,\
                    '\nMethod of calculation           =', self.mthd         ,\
                    '\nDFT code name                   =', self.calculator.__class__.__name__ ,\
                    '\nSpace-group number              =', self.SGN          ,\
                    '\nVolume of equilibrium unit cell =', self.structure.get_volume(), '[a.u^3]',\
                    '\nMaximum Lagrangian strain       =', self.maxstrain          ,\
                    '\nNumber of distorted structures  =', self.NPt ,\
                    '\nMethod                          =', self.mthd,\
                    '\nNumber of elastic constatns     =', self.ECs
                   
    def print_distortions(self):
        for dist in self.distortions:
            dist.print_distortion()
            
    def save(self,file):
        """
        save the  ElasticDistortion object with all data
        """
        pickle.dump(self,open(file,"wb",2))


    
        