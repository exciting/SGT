#!/usr/bin/env python
# -*- coding: utf8 -*-
from ase import Atoms
def StackingFault (primunit,n_slab=3,n_block=1,n_shift=1,shift=1):
        """
          n_slab = 3      #number of the slabs (number of the stacking fault)
          n_block = 1     #number of the blocks in each slab
          n_shift = 18    #total number of the shifts in special direction
          shift = 1       #number of the current shift in special direction
        """
        
        SC = primunit.repeat((n_slab*n_block,1,1))
        positions=SC.get_scaled_positions()
        shiftdir=[[0,.5,.5],[0,.5,0],[0,0,.5]]
        
        for slab in range(0,n_slab-1):  
                
            for i in range(len(positions)):
                startslab=1.0/n_slab*slab
                endslab=1.0/n_slab*(slab+1)
                
                if positions[i,0]>startslab and positions[i,0]<endslab:
                    positions[i]=positions[i]+shiftdir[slab%3]
                    
        SC.set_scaled_positions(positions)
        return SC
        
        