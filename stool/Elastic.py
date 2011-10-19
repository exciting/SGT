
from ase import *

from Elasticlib  import Elastic_setup

class ElasticDistortion(Elastic_setup):
    """
    This class can setup the calculations for determining the elastic constants
    structure:Atoms object
    calculator:calculator object
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
        start all the calculations
        """
        return 0
    def select_data(self):
        """
        select range an fitorder to eliminate errors
        """
        return 0
    def get_elastic_c(self):
        """
        calculate the elastic constants
        """
        return 0
    def make_report(self):
        """
        compile a report file
        """
        print "REPORT"