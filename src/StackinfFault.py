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
        
        SC = primunit.repeat((3,1,1))
        positions=SC.get_positions()
    
        print positions
     
        for i in range(len(positions)):
            if positions[i,0]>1 and positions[i,0]<1.5:
                positions[i]=positions[i]+[0.1,0,0]
        print "after"
        print positions
        
        SC.set_positions(positions)
        return SC
        
        