from sgt  import  SQSscrew ,loadObject
  
WIrSQS= SQSscrew(LC=3.17,NB=1,E1="W",E2="Ir",x=0.2,trials=100)  

WIrSQS.show(file="SQSscrew.png")

WIrSQS.save("SQS.pkl")
new=loadObject("SQS.pkl")
new.show(file="SQSscrewnew.png")