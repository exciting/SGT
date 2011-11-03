import numpy as np
import math as m
import copy
from numpy import linalg
def slab_ducvec(self):
    self.Dhkl=1.0/m.sqrt(np.dot(self.Ghkl,self.Ghkl))
    N=10
    EPS=1.0e-12 
    points=0
    h,k,l=self.miller
    #    !========  Build lattice grid of (hkl) plane  ==========
    LATGRID1=[]
    for x in range(-N,N):
         for y in range (-N,N):
             for z in range (-N,N):
                 if((x*h+y*k+z*l)==0):   
                     LATGRID1.append(np.array([x,y,z]))
                     points+=1 
    """
    =======         Define surface vectors   ===============
    == Define 1st surface vector (the closest to the origin) ==  
    """
    Svec1=LATGRID1[0]
    for i in range(1,points-1):
         vec1=Svec1[0]*self.ducvec[0,:]+Svec1[1]*self.ducvec[1,:]+Svec1[2]*self.ducvec[2,:]
         vec2=LATGRID1[i][0]*self.ducvec[0,:]+LATGRID1[i][1]*self.ducvec[1,:]+LATGRID1[i][2]*self.ducvec[2,:]
         dist1=np.dot(vec1,vec1)
         dist2=np.dot(vec2,vec2)
         if (dist1>dist2 and dist2!=0):
                Svec1=LATGRID1[i]
    """ !== Define 2nd surface vector   =========================
     ! 2nd surface vector should be:
     ! 1) non-collinear to 1st one
     ! 2) Sistem Svec1,Svec2,Ghkl should be right-hand,
     !    i.e. (Svec1,Svec2,Ghkl)>0
     ! 3) 1st or 2nd nearest to the origin
    """
    flag=False
    for i in range(0,points-1):
        #!---     non-collinearity check     --------------
        vec1=Svec1[0]*self.ducvec[0]+Svec1[1]*self.ducvec[1]+Svec1[2]*self.ducvec[2]
        vec2=LATGRID1[i][0]*self.ducvec[0]+LATGRID1[i][1]*self.ducvec[1]+LATGRID1[i][2]*self.ducvec[2]
        vec3= np.cross(vec1,vec2)
        if (abs(vec3[0])<EPS and abs(vec3[1])<EPS and abs(vec3[2])<EPS):
            pass
        else:
            #!---     right-nahd check     --------------------
            vec3=self.Ghkl
            mat =[vec1,vec2,vec3]
            vol=np.linalg.det(mat)
            if (vol<0):
                pass
            else:
                if (flag==False ): 
                    Svec2=LATGRID1[i]
                    flag=True
                else:
                    vec1=Svec2[0]*self.ducvec[0,:]+Svec2[1]*self.ducvec[1,:]+Svec2[2]*self.ducvec[2,:]
                    vec2=LATGRID1[i][0]*self.ducvec[0,:]+LATGRID1[i][1]*self.ducvec[1,:]+LATGRID1[i][2]*self.ducvec[2,:]
                    dist1=np.dot(vec1,vec1)   
                    dist2=np.dot(vec2,vec2)    
            

                    if(dist1>dist2 and dist2!=0):
                        Svec2=LATGRID1[i]
                
                    
                    
                    
        
 