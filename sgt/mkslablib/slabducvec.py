import numpy as np
import math as m
 
from numpy import linalg
from isinteger import isinteger
def slab_ducvec(self):
    self.Dhkl=1.0/m.sqrt(np.dot(self.Ghkl,self.Ghkl))
    N=10
    EPS=1.0e-12 
    points=0
    N_Layers_Max=1000
    h,k,l=self.miller
    slab_vec=np.zeros((3,3))
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
    for i in range(points):
        vec1=Svec1[0]*self.ducvec[0]+Svec1[1]*self.ducvec[1]+Svec1[2]*self.ducvec[2]
        vec2=LATGRID1[i][0]*self.ducvec[0]+LATGRID1[i][1]*self.ducvec[1]+LATGRID1[i][2]*self.ducvec[2]
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
    for i in range(points):
        #!---     non-collinearity check     --------------
        vec1=Svec1[0]*self.ducvec[0]+Svec1[1]*self.ducvec[1]+Svec1[2]*self.ducvec[2]
        vec2=LATGRID1[i][0]*self.ducvec[0]+LATGRID1[i][1]*self.ducvec[1]+LATGRID1[i][2]*self.ducvec[2]
        vec3= np.cross(vec1,vec2)
        if (abs(vec3[0])<EPS and abs(vec3[1])<EPS and abs(vec3[2])<EPS):
            continue
        
        #!---     right-nahd check     --------------------
        vec3=self.Ghkl
        mat =[vec1,vec2,vec3]
        vol=np.linalg.det(mat)
        if (vol<0):
            continue
         
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
        
                    
    """ !===    Define    3rd  vector  =============================
     ! There are 2 approaches to define 3rd vector:
     !    1) 3rd vector is || Ghkl 
     !      (does'nt take into account number of layers explicitly)
     !
     !            \ 
     !       -----O-----O-----O-----O-------
     !            \
     !       ---O-----O-----O-----O-----O---
     !            \ 
     !       -O-----O-----O-----O-----O-----
     !            \  
     !       -----@-----O-----O-----O-----O-
     !
     !     O - a lattice point , ----  - a plane (hkl)
     !         
     !     @ - the origin, \ - the direction of Ghkl
     !        
     !       
     !    2) 3rd vector is the nearest one to Ghkl in the n-th plane
     !       (takes into account number of layers explicitly)
     !
     !            \ 
     !       -----O-----O-----O-----O-------  (n+1)-th plane
     !            \
     !       ---O-----O-----O-----O-----O---   n-th plane
     !          * \
     !          *
     !       ...............................            
     !
     !           *\
     !       -O-----O-----O-----O-----O-----   1-st plane
     !           *\  
     !       -----@-----O-----O-----O-----O-   0-th plane
     !
     !     O - a lattice point , ----  - a plane (hkl)
     !         
     !     @ - the origin, \ - the direction og Ghkl
     !
     !     * - the direction of the 3rd vector 
     !-----------------------------------------------------------            
     """             
    if (self.method==2):
    #!------   \\ Ghkl   ---------------------------------
        hkl=np.matrix([h,k,l])
        M= np.mat(self.ducvec)*np.transpose(np.mat(self.ducvec))
         
        M2= linalg.inv(M)
      
        V=hkl*M2
            
        V=self.Dhkl**2*V
      
 
        flag=False
            
        for  N_Layers in range(1,N_Layers_Max):
            if (flag==True):
                break
                    
            else:
                    
                temp_vector=np.transpose(N_Layers*V[0])
              
                 
                if all(map(isinteger,temp_vector)):
                    x,y,z=map(int,map(round,temp_vector))
                 
                    
                    print 'Good luck! Perpendicular z-vector is found'
                    print 'The thickness of slab is',N_Layers,'Dhkl'
                    
                    Slab_layers=self.layers
                    Svec3=np.array([x,y,z])
                    Svec3=Slab_layers*Svec3  
                     
                    flag=True
                    
                   
                
        slab_vec[0]=Svec1[0]*self.ducvec[0]+Svec1[1]*self.ducvec[1]+Svec1[2]*self.ducvec[2]
        slab_vec[1]=Svec2[0]*self.ducvec[0]+Svec2[1]*self.ducvec[1]+Svec2[2]*self.ducvec[2]
        slab_vec[2]=Svec3[0]*self.ducvec[0]+Svec3[1]*self.ducvec[1]+Svec3[2]*self.ducvec[2]    
        
    elif self.method==3:
        flag=False
        for x in range(N,-N,-1):
            for y in range (N,-N,-1):
                for z in range(N,-N,-1):
                    vec1=x*self.ducvec[0]+y*self.ducvec[1]+z*self.ducvec[2]
                                                
                    result=np.dot(vec1,self.Ghkl)
                        
                    if (abs(result-self.layers)<EPS):
                            
                        if (flag==False):
                                
                            print 'first',x,y,z
                                
                            Svec3=np.array([x,y,z])
                                
                            flag=True
                                
                        else:
                                
                            vec2=vec1-self.layers*(self.Dhkl**2)*self.Ghkl

                            vec3=Svec3[0]*self.ducvec[0]+Svec3[1]*self.ducvec[1]+Svec3[2]*self.ducvec[2]
                                
                            dist1=np.dot(vec1,vec1)
                                
                            dist2=np.dot(vec3,vec3)
                                
                            if(dist1<dist2):
                                
                                Svec3=np.array([x,y,z])
                                    
        slab_vec[0]=Svec1[0]*self.ducvec[0]+Svec1[1]*self.ducvec[1]+Svec1[2]*self.ducvec[2]
        slab_vec[1]=Svec2[0]*self.ducvec[0]+Svec2[1]*self.ducvec[1]+Svec2[2]*self.ducvec[2]
        slab_vec[2]=Svec3[0]*self.ducvec[0]+Svec3[1]*self.ducvec[1]+Svec3[2]*self.ducvec[2]                        
             
    else:
        slab_vec[0]=Svec1[0]*self.ducvec[0]+Svec1[1]*self.ducvec[1]+Svec1[2]*self.ducvec[2]
        slab_vec[1]=Svec2[0]*self.ducvec[0]+Svec2[1]*self.ducvec[1]+Svec2[2]*self.ducvec[2]
        slab_vec[2]=self.layers*(self.Dhkl**2)*self.Ghkl       
    
    return slab_vec
 