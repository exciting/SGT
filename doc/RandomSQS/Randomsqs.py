from RandomSQS import RandomSQS
from ase.io import read, write

WIrSQS= RandomSQS(LC=3.17,NB=1,E1="W",E2="Ir",x=0.2,trials=100)  

WIrSQS.show(file="SQSscrew.png")