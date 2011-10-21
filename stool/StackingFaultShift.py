#!/usr/bin/env python
# -*- coding: utf8 -*-
from ase import Atoms
import numpy as np
def StackingFaultShift (primunit,n_slab=3,n_block=1,n_shift=3.,shift=0.):
        """
         returns supercell structure containing stacking faults
         
          n_slab: integer   
              number of the slabs (number of the stacking fault)
          n_block: integer    
              number of the blocks in each slab
          n_shift: integer   
              total number of the shifts in special direction
          shift: integer         
              number of the current shift in special direction
         
          
        """
        
        SC = primunit.repeat((1,1,n_slab*n_block))
        positions=SC.get_scaled_positions()
        s1=np.array([1.,-1.,0.])*float(shift)/n_shift
        s2=np.array([-2.,-1.,0.])*float(shift)/n_shift
        s3=np.array([1.,2.,0.])*float(shift)/n_shift
        shiftdir=np.array([s1,s2,s3])
	
        posnew= positions.copy()
        for slab in range(n_slab):  
                
            for i in range(0,len(positions)):
                startslab=1.0/n_slab*slab
                endslab=1.0/n_slab*(slab+1)
                
                if positions[i,2]>=startslab:
                    posnew[i,:]+=shiftdir[slab%3,:]
       
        SC.set_scaled_positions(posnew)
	SC.translate(-SC.get_positions()[0])
        return SC
        
        
