import fractions
def norm_miller(miller):
    h,k,l=miller
    if (h==0 and k==0 and l==0): 
            print 'ERROR: incorrect (h,k,l)'
            return True
        
         
    x=abs(h)
    y=abs(k)
    z=abs(l)

    #x/=0, y==0, z==0
    if (x!=0 and y==0 and z==0):
        if(x>1): 
            h=h/x
            

    #x==0, y/=0,z==0
    if (x==0 and y!=0 and z==0):
            if(y>1): 
                k=k/y
            

    # x==0, y==0,z/=0
    if (x==0 and y==0 and z!=0):
        if(z>1):
            l=l/z
            
        

    # x==0, y/=0, z/=0
    if (x==0 and y!=0 and z!=0): 
        i=fractions.gcd(y,z)
             
            
        if(i>1):
                
                k=k/i
                l=l/i
       
    # x/=0, y==0,z/=0
    if (x!=0 and y==0 and z!=0):
        i=fractions.gcd(x,z)
        if(i>1):
                
                h=h/i
                l=l/i
                
         
    # x/=0, y/=0,z==0
    if (x!=0 and y!=0 and z==0):
        i= fractions.gcd(x,y)
            
        if(i>1):
                
                h=h/i
                k=k/i
         
    # x/=0, y/=0, z/=0
    if (x!=0 and y!=0 and z!=0):   
        i=fractions.gcd(x,y)
        j=fractions.gcd(i,z)
             
            
        if(j>1):
            
                h=h/j
                k=k/j
                l=l/j
    return [h,k,l]