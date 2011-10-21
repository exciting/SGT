from RandomSQS import  SQSscrew
from StackingFaultShift import StackingFaultShift
from Elastic import   ElasticDistortion
import pickle

def loadObject(file):
    return pickle.load(open(file,"rb"))