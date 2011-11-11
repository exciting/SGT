from RandomSQS import  SQSscrew
from StackingFaultShift import StackingFaultShift
from ElaStic import   ElaSticDistortion
from mkslab import slab
import pickle

def loadObject(file): 
    return pickle.load(open(file,"rb"))