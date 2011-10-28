
from ase import *

from Elasticlib  import Elastic_setup,  elastic_select_data , ElaSicResult

import numpy as np
 
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
        self.make_report()
    
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
    def get_elastic_c(self,selectdata=None):
        """
        calculate the elastic constants. It takes one argument, the selectdata list.
        It may be a pair of etamax and fit order, which configures all elastic constants,
        or an array of pairs. Each pair configures one of the elastic constants.
        """
        
        self.fitdatasel=[]
        if not(selectdata):
            selectdata=[self.maxstrain, self.order+4]
        else:
            if np.shape(selectdata)==(2,):
                for dist in range(self.ECs):
                    self.fitdatasel.append(selectdata)
            elif np.shape(selectdata)==(self.ECs,2):
                 for dist in range(self.ECs):
                    self.fitdatasel.append(selectdata[dist])
                     
            else:
                print "The selectdata array must have only one eta,order pair for" 
                print "all elastic constants or a list of pairs for each of the", self.ECs," constants"
                return 1
        ElaSicResult(self)
        
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


    
        