from RandomSQS import  *
from StackingFaultShift import *
from Elastic import  *
import pickle

def loadObject(file):
    return pickle.load(open(file,"rb"))