#!/usr/bin/env python
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%------------------------------------------------- ElaStic_Setup_exciting ------------------------------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#
# AUTHORS:
# Rostam Golesorkhtabar and Pasquale Pavone 
# r.golesorkhtabar@mcl.at
# 
# DATE:
# Sun Jan 01 00:00:00 2012
#
# SYNTAX:
# python ElaStic_Setup_exciting.py
#        ElaStic_Setup_exciting
# 
# EXPLANATION:
# 
#________________________________________________________________________________________________________________________________
from sys   import stdin
from numpy import *
import numpy as np
from distortion import distortion
import shutil
import glob
import math
import sys
try:
    from pyspglib import spglib
except ImportError:
    print "#######################################################"
    print "In order to use the Elastic tools you need  pyspglib"
    print "#######################################################"
    print "http://spglib.sourceforge.net/pyspglibForASE/"
    
class Elastic_setup():
    
    def setup(self,structure,calculator,order, maxstrain,distortions):
        
        self.SGN=spglib.get_symmetry_dataset(structure, symprec=1e-5)["number"]
        self.structure=structure
        self.calculator=calculator
        self.structure.set_calculator(calculator)
        self.order=order
        self.maxstrain=maxstrain
        self.NPt=distortions
        self.distortions=[]
   
        if ( hasattr(self.calculator, 'get_stress') ): self.mthd ='Stress' 
        else:  self.mthd = 'Energy'
        
        #%%%--- DICTIONARIS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        self.Ls_Dic={                       \
        '01':[ 1., 1., 1., 0., 0., 0.],\
        '02':[ 1., 0., 0., 0., 0., 0.],\
        '03':[ 0., 1., 0., 0., 0., 0.],\
        '04':[ 0., 0., 1., 0., 0., 0.],\
        '05':[ 0., 0., 0., 2., 0., 0.],\
        '06':[ 0., 0., 0., 0., 2., 0.],\
        '07':[ 0., 0., 0., 0., 0., 2.],\
        '08':[ 1., 1., 0., 0., 0., 0.],\
        '09':[ 1., 0., 1., 0., 0., 0.],\
        '10':[ 1., 0., 0., 2., 0., 0.],\
        '11':[ 1., 0., 0., 0., 2., 0.],\
        '12':[ 1., 0., 0., 0., 0., 2.],\
        '13':[ 0., 1., 1., 0., 0., 0.],\
        '14':[ 0., 1., 0., 2., 0., 0.],\
        '15':[ 0., 1., 0., 0., 2., 0.],\
        '16':[ 0., 1., 0., 0., 0., 2.],\
        '17':[ 0., 0., 1., 2., 0., 0.],\
        '18':[ 0., 0., 1., 0., 2., 0.],\
        '19':[ 0., 0., 1., 0., 0., 2.],\
        '20':[ 0., 0., 0., 2., 2., 0.],\
        '21':[ 0., 0., 0., 2., 0., 2.],\
        '22':[ 0., 0., 0., 0., 2., 2.],\
        '23':[ 0., 0., 0., 2., 2., 2.],\
        '24':[-1., .5, .5, 0., 0., 0.],\
        '25':[ .5,-1., .5, 0., 0., 0.],\
        '26':[ .5, .5,-1., 0., 0., 0.],\
        '27':[ 1.,-1., 0., 0., 0., 0.],\
        '28':[ 1.,-1., 0., 2., 0., 0.],\
        '29':[ 1.,-1., 0., 0., 0., 2.],\
        '30':[ 1., 0.,-1., 0., 2., 0.],\
        '31':[ 0., 1.,-1., 0., 0., 2.],\
        '32':[ 1., 1.,-1., 2., 2., 2.],\
        '33':[ 1., 0., 0., 2., 2., 0.],\
        '34':[ 0., 1., 0., 2., 2., 0.],\
        '35':[ 1., 1., 0., 2., 2., 0.],\
        '36':[ 1., 1., 0., 2., 0., 0.],\
        '37':[ 1., 1.,-1., 0., 0., 0.],\
        '38':[ 1., 1., 1.,-2.,-2.,-2.],\
        '39':[ 1., 2., 3., 4., 5., 6.],\
        '40':[-2., 1., 4.,-3., 6.,-5.],\
        '41':[ 3.,-5.,-1., 6., 2.,-4.],\
        '42':[-4.,-6., 5., 1.,-3., 2.],\
        '43':[ 5., 4., 6.,-2.,-1.,-3.],\
        '44':[-6., 3.,-2., 5.,-4., 1.]}
        
        self.Ls_str={                                     \
        '01':'(  eta,  eta,  eta,  0.0,  0.0,  0.0)',\
        '02':'(  eta,  0.0,  0.0,  0.0,  0.0,  0.0)',\
        '03':'(  0.0,  eta,  0.0,  0.0,  0.0,  0.0)',\
        '04':'(  0.0,  0.0,  eta,  0.0,  0.0,  0.0)',\
        '05':'(  0.0,  0.0,  0.0, 2eta,  0.0,  0.0)',\
        '06':'(  0.0,  0.0,  0.0,  0.0, 2eta,  0.0)',\
        '07':'(  0.0,  0.0,  0.0,  0.0,  0.0, 2eta)',\
        '08':'(  eta,  eta,  0.0,  0.0,  0.0,  0.0)',\
        '09':'(  eta,  0.0,  eta,  0.0,  0.0,  0.0)',\
        '10':'(  eta,  0.0,  0.0, 2eta,  0.0,  0.0)',\
        '11':'(  eta,  0.0,  0.0,  0.0, 2eta,  0.0)',\
        '12':'(  eta,  0.0,  0.0,  0.0,  0.0, 2eta)',\
        '13':'(  0.0,  eta,  eta,  0.0,  0.0,  0.0)',\
        '14':'(  0.0,  eta,  0.0, 2eta,  0.0,  0.0)',\
        '15':'(  0.0,  eta,  0.0,  0.0, 2eta,  0.0)',\
        '16':'(  0.0,  eta,  0.0,  0.0,  0.0, 2eta)',\
        '17':'(  0.0,  0.0,  eta, 2eta,  0.0,  0.0)',\
        '18':'(  0.0,  0.0,  eta,  0.0, 2eta,  0.0)',\
        '19':'(  0.0,  0.0,  eta,  0.0,  0.0, 2eta)',\
        '20':'(  0.0,  0.0,  0.0, 2eta, 2eta,  0.0)',\
        '21':'(  0.0,  0.0,  0.0, 2eta,  0.0, 2eta)',\
        '22':'(  0.0,  0.0,  0.0,  0.0, 2eta, 2eta)',\
        '23':'(  0.0,  0.0,  0.0, 2eta, 2eta, 2eta)',\
        '24':'( -eta,.5eta,.5eta,  0.0,  0.0,  0.0)',\
        '25':'(.5eta, -eta,.5eta,  0.0,  0.0,  0.0)',\
        '26':'(.5eta,.5eta, -eta,  0.0,  0.0,  0.0)',\
        '27':'(  eta, -eta,  0.0,  0.0,  0.0,  0.0)',\
        '28':'(  eta, -eta,  0.0, 2eta,  0.0,  0.0)',\
        '29':'(  eta, -eta,  0.0,  0.0,  0.0, 2eta)',\
        '30':'(  eta,  0.0, -eta,  0.0, 2eta,  0.0)',\
        '31':'(  0.0,  eta, -eta,  0.0,  0.0, 2eta)',\
        '32':'(  eta,  eta, -eta, 2eta, 2eta, 2eta)',\
        '33':'(  eta,  0.0,  0.0, 2eta, 2eta,  0.0)',\
        '34':'(  0.0,  eta,  0.0, 2eta, 2eta,  0.0)',\
        '35':'(  eta,  eta,  0.0, 2eta, 2eta,  0.0)',\
        '36':'(  eta,  eta,  0.0, 2eta,  0.0,  0.0)',\
        '37':'(  eta,  eta, -eta,  0.0,  0.0,  0.0)',\
        '38':'(  eta,  eta,  eta,-2eta,-2eta,-2eta)',\
        '39':'( 1eta, 2eta, 3eta, 4eta, 5eta, 6eta)',\
        '40':'(-2eta, 1eta, 4eta,-3eta, 6eta,-5eta)',\
        '41':'( 3eta,-5eta,-1eta, 6eta, 2eta,-4eta)',\
        '42':'(-4eta,-6eta, 5eta, 1eta,-3eta, 2eta)',\
        '43':'( 5eta, 4eta, 6eta,-2eta,-1eta,-3eta)',\
        '44':'(-6eta, 3eta,-2eta, 5eta,-4eta, 1eta)'}
        
        self.LC_Dic = {              \
        'CI' :'Cubic I'        ,\
        'CII':'Cubic II'       ,\
        'HI' :'Hexagonal I'    ,\
        'HII':'Hexagonal II'   ,\
        'RI' :'Rhombohedral I' ,\
        'RII':'Rhombohedral II',\
        'TI' :'Tetragonal I'   ,\
        'TII':'Tetragonal II'  ,\
        'O'  :'Orthorhombic'   ,\
        'M'  :'Monoclinic'     ,\
        'N'  :'Triclinic'} 
        #--------------------------------------------------------------------------------------------------------------------------------
        
    
        
        
        #%%%--- Reading the order of the elastic constants ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
     
        if (self.order != 2 and self.order != 3 ):
            sys.exit("\n    ...Oops ERROR: Choose '2' or '3' \n")
        
       
        
        if (1 <= self.SGN and self.SGN <= 2):      # Triclinic
            self.LC = 'N'
            if (self.order == 2): self.ECs = 21
            if (self.order == 3): self.ECs = 56  
        
        elif(3 <= self.SGN and self.SGN <= 15):    # Monoclinic
            self.LC = 'M'
            if (self.order == 2): self.ECs = 13
            if (self.order == 3): self.ECs = 32 
        
        elif(16 <= self.SGN and self.SGN <= 74):   # Orthorhombic
            self.LC = 'O'
            if (self.order == 2): self.ECs =  9
            if (self.order == 3): self.ECs = 20 
        
        elif(75 <= self.SGN and self.SGN <= 88):   # Tetragonal II
            self.LC = 'TII'
            if (self.order == 2): self.ECs =  7
            if (self.order == 3): self.ECs = 16
          
        elif(89 <= self.SGN and self.SGN <= 142):  # Tetragonal I
            self.LC = 'TI'
            if (self.order == 2): self.ECs =  6
            if (self.order == 3): self.ECs = 12  
        
        elif(143 <= self.SGN and self.SGN <= 148): # Rhombohedral II 
            self.LC = 'RII'
            if (self.order == 2): self.ECs =  7
            if (self.order == 3): self.ECs = 20
        
        elif(149 <= self.SGN and self.SGN <= 167): # Rhombohedral I
            self.LC = 'RI'
            if (self.order == 2): self.ECs =  6
            if (self.order == 3): self.ECs = 14
        
        elif(168 <= self.SGN and self.SGN <= 176): # Hexagonal II
            self.LC = 'HII'
            if (self.order == 2): self.ECs =  5
            if (self.order == 3): self.ECs = 12
        
        elif(177 <= self.SGN and self.SGN <= 194): # Hexagonal I
            self.LC = 'HI'
            if (self.order == 2): self.ECs =  5
            if (self.order == 3): self.ECs = 10
        
        elif(195 <= self.SGN and self.SGN <= 206): # Cubic II
            self.LC = 'CII'
            if (self.order == 2): self.ECs =  3
            if (self.order == 3): self.ECs =  8
        
        elif(207 <= self.SGN and self.SGN <= 230): # Cubic I
            self.LC = 'CI'
            if (self.order == 2): self.ECs =  3
            if (self.order == 3): self.ECs =  6
        else: sys.exit('\n     ... Oops ERROR: WRONG Space Group Number !?!?!?    \n')
        
        if (self.order == 2): order = 'second'
        if (self.order == 3): order = 'third'
        
        
        #%%%--- Checking the maximum Lagrangian strain ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        if (1 <  self.maxstrain or  self.maxstrain < 0):
            sys.exit('\n     ... Oops ERROR: The maximum Lagrangian strain is out of range !!!!!!    \n')
        
        self.maxstrain = round( self.maxstrain, 3)
        print '     The maximum Lagrangian strain is '+ str( self.maxstrain) + '\n'
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Checking the number of the distorted structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        
        if (self.NPt < 5):
            sys.exit('\n     ... Oops ERROR: The NUMBER of the distorted structures < 5 !!!!!!    \n')
        if (99 < self.NPt):
            sys.exit('\n     ... Oops ERROR: The NUMBER of the distorted structures > 99 !!!!!!   \n')
        
        if (self.NPt%2 == 0):
            self.NPt   += 1
        print '     The number of the distorted structures is '+ str(self.NPt) + '\n'
        
        ptn = int((self.NPt-1)/2)
        
        if (self.mthd == 'Energy'): interval = 0.0001
        if (self.mthd == 'Stress'): interval = 0.00001
        if ( self.maxstrain/ptn <= interval):
            sys.exit('\n     ... Oops ERROR: The interval of the strain values is less than '+ str(interval) +'\
                      \n                     Choose a larger maximum Lagrangian strain or a less number of distorted structures.\n')
        #--------------------------------------------------------------------------------------------------------------------------------
        
        
        #--------------------------------------------------------------------------------------------------------------------------------
        
     
        
        V0=  self.structure.get_cell()
   
      
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Writing the INFO_ElaStic file ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #-------------------------------------------------------
        
        if (self.mthd == 'Energy'):
            if (self.order == 2):
                if (self.LC == 'CI' or \
                    self.LC == 'CII'):
                    Lag_strain_list = ['01','08','23']
                if (self.LC == 'HI' or \
                    self.LC == 'HII'):
                    Lag_strain_list = ['01','26','04','03','17']
                if (self.LC == 'RI'):
                    Lag_strain_list = ['01','08','04','02','05','10']
                if (self.LC == 'RII'):
                    Lag_strain_list = ['01','08','04','02','05','10','11']
                if (self.LC == 'TI'):
                    Lag_strain_list = ['01','26','27','04','05','07']
                if (self.LC == 'TII'):
                    Lag_strain_list = ['01','26','27','28','04','05','07']
                if (self.LC == 'O'):
                    Lag_strain_list = ['01','26','25','27','03','04','05','06','07']
                if (self.LC == 'M'):
                    Lag_strain_list = ['01','25','24','28','29','27','20','12','03','04','05','06','07']
                if (self.LC == 'N'):
                    Lag_strain_list = ['02','03','04','05','06','07','08','09','10','11',\
                                       '12','13','14','15','16','17','18','19','20','21','22']
        
            if (self.order == 3):
                if (self.LC == 'CI'):
                    Lag_strain_list = ['01','08','23','32','10','11']
                if (self.LC == 'CII'):
                    Lag_strain_list = ['01','08','23','32','10','11','12','09']
                if (self.LC == 'HI'):
                    Lag_strain_list = ['01','26','04','03','17','30','08','02','10','14']
                if (self.LC == 'HII'):
                    Lag_strain_list = ['01','26','04','03','17','30','08','02','10','14','12','31']
                if (self.LC == 'RI'):
                    Lag_strain_list = ['01','08','04','02','05','10','11','26','09','03','17','34','33','35']
                if (self.LC == 'RII'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'TI'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'TII'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'O'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'M'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'N'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
        
        if (self.mthd == 'Stress'):
            if (self.order == 2):
                if (self.LC == 'CI' or \
                    self.LC == 'CII'):
                    Lag_strain_list = ['36']
                if (self.LC == 'HI' or \
                    self.LC == 'HII'):
                    Lag_strain_list = ['36','38']
                if (self.LC == 'RI' or \
                    self.LC == 'RII'):
                    Lag_strain_list = ['36','38']
                if (self.LC == 'TI' or \
                    self.LC == 'TII'):
                    Lag_strain_list = ['36','38']
                if (self.LC == 'O'):
                    Lag_strain_list = ['36','38','40']
                if (self.LC == 'M'):
                    Lag_strain_list = ['36','37','38','39','40']
                if (self.LC == 'N'):
                    Lag_strain_list = ['36','37','38','39','40','41']
        
            if (self.order == 3):
                if (self.LC == 'CI'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'CII'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'HI'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'HII'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'RI'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'RII'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'TI'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'TII'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'O'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'M'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
                if (self.LC == 'N'):
                    sys.exit('\n     ... Oops SORRY: Not implemented yet. \n')
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
        ###---------------------------------------------------- Structures maker -----------------------------------------------------###
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
        fdis = open('Distorted_Parameters','w')
        cont1= 0
        for i in Lag_strain_list:
            Ls_list= self.Ls_Dic[i]
        
            cont1  = cont1 + 1
            if (cont1 < 10):
                Dstn = 'Dst0'+str(cont1)
            else:
                Dstn = 'Dst' +str(cont1)
            
             
             
        
            print>>fdis, Dstn+', Lagrangian strain = ' + self.Ls_str[i]
        
            cont2 = 0
            for s in range(-ptn, ptn+1):
                r =  self.maxstrain*s/ptn
                if (s==0):
                    if (self.mthd == 'Energy'): r = 0.0001
                    if (self.mthd == 'Stress'): r = 0.00001
        
                Ls = zeros(6)
                for xindex in range(6):
                    Ls[xindex] = Ls_list[xindex]
                Lv = r*Ls
        #-------Lagrangian strain to physical strain (eta = eps + 0.5*eps*esp)-----------------------------------------------------------
                eta_matrix      = zeros((3,3))
        
                eta_matrix[0,0] = Lv[0]
                eta_matrix[0,1] = Lv[5]/2.
                eta_matrix[0,2] = Lv[4]/2.
                
                eta_matrix[1,0] = Lv[5]/2.
                eta_matrix[1,1] = Lv[1]
                eta_matrix[1,2] = Lv[3]/2.
        
                eta_matrix[2,0] = Lv[4]/2.
                eta_matrix[2,1] = Lv[3]/2.
                eta_matrix[2,2] = Lv[2]
        
                norm       = 1.0
        
                eps_matrix = eta_matrix
                if (linalg.norm(eta_matrix) > 0.7):
                    sys.exit('\n     ... Oops ERROR: Too large deformation!\n') 
        
                while( norm > 1.e-10 ):
                    x          = eta_matrix - dot(eps_matrix, eps_matrix)/2.
                    norm       = linalg.norm(x - eps_matrix)      
                    eps_matrix = x
        
        #--------------------------------------------------------------------------------------------------------------------------------
                i_matrix   = array([[1., 0., 0.],
                                    [0., 1., 0.], 
                                    [0., 0., 1.]])
                def_matrix = i_matrix + eps_matrix
                M_new      = dot(self.structure.get_cell(), def_matrix)
        
                cont2 = cont2 + 1
                if (cont2 < 10):
                    Dstn_cont2 = Dstn+ '_0'+str(cont2)
                else:
                    Dstn_cont2 = Dstn+ '_' +str(cont2)
        
                print>>fdis, Dstn_cont2 + ',  eta = ' + str(r)
        
                 
                 
                self.distortions.append(
                    distortion(self.structure, calculator=calculator,cell=M_new,eta=[cont2,r], LagrangeS= [cont1, self.Ls_str[i]])
                    )
        
       
        #--------------------------------------------------------------------------------------------------------------------------------
                 
        print    "total distortions", len(self.distortions)
        fdis.close()
         #--------------------------------------------------------------------------------------------------------------------------------
    
    