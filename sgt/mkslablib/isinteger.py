def isinteger(x):
    eps=1.0e-12
    if(abs(round(x)-x)<eps):
        return True
    else:
        return False
