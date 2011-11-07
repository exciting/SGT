#!/usr/bin/env python 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%% ---------------------------------------------- ElaStic_Result_Energy_2nd ---------------------------------------------- %%%#
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
# python ElaStic_Result_Energy_2nd.py
#        ElaStic_Result_Energy_2nd
#
# EXPLANATION:
# 
#________________________________________________________________________________________________________________________________

from sys   import stdin
from numpy import *
import numpy as np
import subprocess
import os.path
import shutil
import math
import time
import sys
import os
def ElaSicResult(self):
    #%%%%%%%%--- CONSTANTS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    _e     = 1.602176565e-19              # elementary charge
    Bohr   = 5.291772086e-11              # a.u. to meter
    Ryd2eV = 13.605698066                 # Ryd to eV
    cnvrtr = (_e*Ryd2eV)/(1e9*Bohr**3)    # Ryd/[a.u.^3] to GPa
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Dictionaries ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    #--------------------------------------------------------------------------------------------------------------------------------
    
   
    lineseparator=' '
    for i in range(0,79):
        lineseparator=lineseparator+'%'
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    ### ------------------------------ Calculating the second derivative and Cross-Validation error ----------------------------- ###
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
   
    
    # Range of Deformation
   
    
   
    
    A2 = []
    for i in range(1, self.ECs+1):
        if (i<10):
            Dstn = 'Dst0'+str(i)
        else:
            Dstn = 'Dst' +str(i)
    
       
        mdri   = abs( self.fitdatasel[i-1][0])
        ordri  = int(abs(self.fitdatasel[i-1][1]))
    
        
    
        strain = []
        energy = []
        for k in range(0, self.NPt-1, 2):
            if (-mdri <= self.distortions[(i-1)*self.NPt+k].eta[1] 
                and self.distortions[(i-1)*self.NPt+k].eta[1] <= mdri):
                 
                strain.append(self.distortions[(i-1)*self.NPt+k].eta[1])
                energy.append(self.distortions[(i-1)*self.NPt+k].atoms.get_potential_energy()) 
          
        if (len(strain) < ordri+1):  
            sys.exit('\n     ... Oops ERROR: NOT enough energy points in '+Dstn+'_Energy.dat for '+str(ordri)+' order polynomial fit.\n')
    
        coeffs = np.polyfit(strain, energy, ordri)
        A2.append(coeffs[ordri-2])
    
    A2 = np.array(A2)
    if (len(A2) != self.ECs):
        sys.exit("\n     ... Oops ERROR: The number of data in the 'ElaStic_2nd.in' file is NOT equal to "+str(self.ECs)+"\n")
    
    self.C = zeros((6,6))
    #%%%--- Cubic structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (self.LC == 'CI' or \
        self.LC == 'CII'):
        self.C[0,0] =-2.*(A2[0]-3.*A2[1])/3.
        self.C[1,1] = self.C[0,0]
        self.C[2,2] = self.C[0,0]
        self.C[3,3] = A2[2]/6.
        self.C[4,4] = self.C[3,3]
        self.C[5,5] = self.C[3,3]
        self.C[0,1] = (2.*A2[0]-3.*A2[1])/3.
        self.C[0,2] = self.C[0,1]
        self.C[1,2] = self.C[0,1]
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Hexagonal structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (self.LC == 'HI' or \
        self.LC == 'HII'):
        self.C[0,0] = 2.*A2[3]
        self.C[0,1] = 2./3.*A2[0] + 4./3.*A2[1] - 2.*A2[2] - 2.*A2[3]
        self.C[0,2] = 1./6.*A2[0] - 2./3.*A2[1] + 0.5*A2[2]
        self.C[1,1] = self.C[0,0]
        self.C[1,2] = self.C[0,2]
        self.C[2,2] = 2.*A2[2]
        self.C[3,3] =-0.5*A2[2] + 0.5*A2[4]
        self.C[4,4] = self.C[3,3]
        self.C[5,5] = .5*(self.C[0,0] - self.C[0,1])
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Rhombohedral I structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (self.LC == 'RI'):
        self.C[0,0] = 2.*A2[3]
        self.C[0,1] = A2[1]- 2.*A2[3]
        self.C[0,2] = .5*( A2[0] - A2[1] - A2[2])
        self.C[0,3] = .5*(-A2[3] - A2[4] + A2[5])
        self.C[1,1] = self.C[0,0]
        self.C[1,2] = self.C[0,2]
        self.C[1,3] =-self.C[0,3]
        self.C[2,2] = 2.*A2[2]
        self.C[3,3] = .5*A2[4]
        self.C[4,4] = self.C[3,3]
        self.C[4,5] = self.C[0,3]
        self.C[5,5] = .5*(self.C[0,0] - self.C[0,1])
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Rhombohedral II structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (self.LC == 'RII'):
        self.C[0,0] = 2.*A2[3]
        self.C[0,1] = A2[1]- 2.*A2[3]
        self.C[0,2] = .5*( A2[0] - A2[1] - A2[2])
        self.C[0,3] = .5*(-A2[3] - A2[4] + A2[5])
        self.C[0,4] = .5*(-A2[3] - A2[4] + A2[6])
        self.C[1,1] = self.C[0,0]
        self.C[1,2] = self.C[0,2]
        self.C[1,3] =-self.C[0,3]
        self.C[1,4] =-self.C[0,4]    
        self.C[2,2] = 2.*A2[2]
        self.C[3,3] = .5*A2[4]
        self.C[3,5] =-self.C[0,4]
        self.C[4,4] = self.C[3,3]
        self.C[4,5] = self.C[0,3]
        self.C[5,5] = .5*(self.C[0,0] - self.C[0,1])
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Tetragonal I structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (self.LC == 'TI'):
        self.C[0,0] = (A2[0]+2.*A2[1])/3.+.5*A2[2]-A2[3]
        self.C[0,1] = (A2[0]+2.*A2[1])/3.-.5*A2[2]-A2[3]
        self.C[0,2] = A2[0]/6.-2.*A2[1]/3.+.5*A2[3]
        self.C[1,1] = self.C(0,0)
        self.C[1,2] = self.C(0,2)
        self.C[2,2] = 2.*A2[3]
        self.C[3,3] = .5*A2[4]
        self.C[4,4] = self.C[3,3]
        self.C[5,5] = .5*A2[5]
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Tetragonal II structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (self.LC == 'TII'):
        self.C[0,0] = (A2[0]+2.*A2[1])/3.+.5*A2[2]-A2[4]
        self.C[1,1] = self.C[0,0]
        self.C[0,1] = (A2[0]+2.*A2[1])/3.-.5*A2[2]-A2[4]
        self.C[0,2] = A2[0]/6.-(2./3.)*A2[1]+.5*A2[4]
        self.C[0,5] = (-A2[2]+A2[3]-A2[6])/4.
        self.C[1,2] = self.C[0,2]
        self.C[1,5] =-self.C[0,5]
        self.C[2,2] = 2.*A2[4]
        self.C[3,3] = .5*A2[5]
        self.C[4,4] = self.C[3,3]
        self.C[5,5] = .5*A2[6]
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Orthorhombic structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (self.LC == 'O'):
        self.C[0,0] = 2.*A2[0]/3.+4.*A2[1]/3.+A2[3]-2.*A2[4]-2.*A2[5]
        self.C[0,1] = 1.*A2[0]/3.+2.*A2[1]/3.-.5*A2[3]-A2[5]
        self.C[0,2] = 1.*A2[0]/3.-2.*A2[1]/3.+4.*A2[2]/3.-.5*A2[3]-A2[4]
        self.C[1,1] = 2.*A2[4]
        self.C[1,2] =-2.*A2[1]/3.-4.*A2[2]/3.+.5*A2[3]+A2[4]+A2[5]
        self.C[2,2] = 2.*A2[5]
        self.C[3,3] = .5*A2[6]
        self.C[4,4] = .5*A2[7]
        self.C[5,5] = .5*A2[8]
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Monoclinic structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (self.LC == 'M'):
        self.C[0,0] = 2.*A2[0]/3.+8.*(A2[1]+A2[2])/3.-2.*(A2[5]+A2[8]+A2[9])
        self.C[0,1] = A2[0]/3.+4.*(A2[1]+A2[2])/3.-2.*A2[5]-A2[9]
        self.C[0,2] =(A2[0]-4.*A2[2])/3.+A2[5]-A2[8]
        self.C[0,5] =-1.*A2[0]/6.-2.*(A2[1]+A2[2])/3.+.5*(A2[5]+A2[7]+A2[8]+A2[9]-A2[12])
        self.C[1,1] = 2.*A2[8]
        self.C[1,2] =-4.*(2.*A2[1]+A2[2])/3.+2.*A2[5]+A2[8]+A2[9]+A2[12]
        self.C[1,5] =-1.*A2[0]/6.-2.*(A2[1]+A2[2])/3.-.5*A2[3]+A2[5]+.5*(A2[7]+A2[8]+A2[9])
        self.C[2,2] = 2.*A2[9]
        self.C[2,5] =-1.*A2[0]/6.+2.*A2[1]/3.-.5*(A2[3]+A2[4]-A2[7]-A2[8]-A2[9]-A2[12])
        self.C[3,3] = .5*A2[10]
        self.C[3,4] = .25*(A2[6]-A2[10]-A2[11])
        self.C[4,4] = .5*A2[11]
        self.C[5,5] = .5*A2[12]
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Triclinic structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (self.LC == 'N'):
        self.C[0,0] = 2.*A2[0]
        self.C[0,1] = 1.*(-A2[0]-A2[1]+A2[6])
        self.C[0,2] = 1.*(-A2[0]-A2[2]+A2[7])
        self.C[0,3] = .5*(-A2[0]-A2[3]+A2[8]) 
        self.C[0,4] = .5*(-A2[0]+A2[9]-A2[4])
        self.C[0,5] = .5*(-A2[0]+A2[10]-A2[5])
        self.C[1,1] = 2.*A2[1]
        self.C[1,2] = 1.*(A2[11]-A2[1]-A2[2])
        self.C[1,3] = .5*(A2[12]-A2[1]-A2[3])
        self.C[1,4] = .5*(A2[13]-A2[1]-A2[4])
        self.C[1,5] = .5*(A2[14]-A2[1]-A2[5])
        self.C[2,2] = 2.*A2[2] 
        self.C[2,3] = .5*(A2[15]-A2[2]-A2[3])
        self.C[2,4] = .5*(A2[16]-A2[2]-A2[4])
        self.C[2,5] = .5*(A2[17]-A2[2]-A2[5])
        self.C[3,3] = .5*A2[3]
        self.C[3,4] = .25*(A2[18]-A2[3]-A2[4])
        self.C[3,5] = .25*(A2[19]-A2[3]-A2[5])
        self.C[4,4] = .5*A2[4]
        self.C[4,5] = .25*(A2[20]-A2[4]-A2[5])
        self.C[5,5] = .5*A2[5]
    #--------------------------------------------------------------------------------------------------------------------------------
    
    #%%%--- Calculating the elastic moduli ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    CONV = cnvrtr * 1.
    
    for i in range(5):
        for j in range(i+1,6):
            self.C[j,i] = self.C[i,j] 
    
    self.C = self.C * CONV/ self.structure.get_volume()
    
    BV = (self.C[0,0]+self.C[1,1]+self.C[2,2]+2*(self.C[0,1]+self.C[0,2]+self.C[1,2]))/9
    GV = ((self.C[0,0]+self.C[1,1]+self.C[2,2])-(self.C[0,1]+self.C[0,2]+self.C[1,2])+3*(self.C[3,3]+self.C[4,4]+self.C[5,5]))/15
    EV = (9*BV*GV)/(3*BV+GV)
    nuV= (1.5*BV-GV)/(3*BV+GV)
    S  = linalg.inv(self.C)
    BR = 1/(S[0,0]+S[1,1]+S[2,2]+2*(S[0,1]+S[0,2]+S[1,2]))
    GR =15/(4*(S[0,0]+S[1,1]+S[2,2])-4*(S[0,1]+S[0,2]+S[1,2])+3*(S[3,3]+S[4,4]+S[5,5]))
    ER = (9*BR*GR)/(3*BR+GR)
    nuR= (1.5*BR-GR)/(3*BR+GR)
    BH = 0.50*(BV+BR)
    GH = 0.50*(GV+GR)
    EH = (9.*BH*GH)/(3.*BH+GH)
    nuH= (1.5*BH-GH)/(3.*BH+GH)
    AVR= 100.*(GV-GR)/(GV+GR)
    
    self.BV=BV
    self.GV= GV
    self.EV= EV
    self.nuV= nuV
    self.S= S
    self.BR= BR
    self.GR= GR
    self.ER= ER
    self.nuR= nuR
    self.BH=BH
    self.GH=GH
    self.EH=EH
    self.nuH=nuH
    self.AVR=AVR
    fo = open('ElaStic_2nd.out','w')
    self.havecontants=True
    print_elastic_moduli(self,fo)
    return self.C
    #--------------------------------------------------------------------------------------------------------------------------------
def print_elastic_moduli(self,fo):    
    #%%%--- Writing the output file ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lineseparator=' '
    for i in range(0,79):
        lineseparator=lineseparator+'%' 
    print >>fo,  '    The output of ElaStic code                                              \n'\
                 '    Today is '+ time.asctime() +                                           '\n'\
                                                                                             '\n'\
                 '    Symmetry of the second-order elastic constant matrix in Voigt notation. \n'\
                    + head[self.LC] +                                                             '\n'\
                 '    Elastic constant (stiffness) matrix in GPa :                            \n'
    
    for i in range(0,6):
        print >>fo, '',
        for j in range(0,6):
            print >>fo, '%11.1f'%(self.C[i,j]),
        print >>fo
    
    print >>fo,'\n\n    Elastic compliance matrix in 1/GPa : \n'
    
    for i in range(0,6):
        print >>fo, '',
        for j in range(0,6):
            print >>fo, '%11.5f'%(self.S[i,j]),
        print >>fo
    
    print >>fo, '\n'+ lineseparator +'\n'
    
    print >>fo, '    Voigt bulk  modulus, K_V = {0}  GPa'.format('%8.2f'%(self.BV))
    print >>fo, '    Voigt shear modulus, G_V = {0}  GPa'.format('%8.2f'%(self.GV)) + '\n'
    
    print >>fo, '    Reuss bulk  modulus, K_R = {0}  GPa'.format('%8.2f'%(self.BR))
    print >>fo, '    Reuss shear modulus, G_R = {0}  GPa'.format('%8.2f'%(self.GR)) + '\n'
    
    print >>fo, '    Hill bulk  modulus, K_H  = {0}  GPa'.format('%8.2f'%(self.BH))
    print >>fo, '    Hill shear modulus, G_H  = {0}  GPa'.format('%8.2f'%(self.GH))
    
    print >>fo, '\n'+ lineseparator +'\n'
    
    print >>fo, '    Voigt Young modulus,  E_V = {0}  GPa'.format('%8.2f'%(self.EV))
    print >>fo, '    Voigt Poisson ratio, nu_V = {0}'     .format('%8.2f'%(self.nuV)) + '\n'
    
    print >>fo, '    Reuss Young modulus,  E_R = {0}  GPa'.format('%8.2f'%(self.ER))
    print >>fo, '    Reuss Poisson ratio, nu_R = {0}'     .format('%8.2f'%(self.nuR)) + '\n'
    
    print >>fo, '    Hill Young modulus,  E_H  = {0}  GPa'.format('%8.2f'%(self.EH))
    print >>fo, '    Hill Poisson ratio, nu_H  = {0}'     .format('%8.2f'%(self.nuH))
    
    print >>fo, '\n'+ lineseparator +'\n'
    
    print >>fo, '    Elastic Anisotropy in polycrystalline, AVR = {0} %'.format('%8.3f'%(self.AVR))
    
    print >>fo, '\n'+ lineseparator +'\n'
    
    print >>fo, '    Eigenvalues of elastic constant (stiffness) matrix:   \n'
    
    eigval=linalg.eig(self.C)
    for i in range(6):
        print >>fo,'%16.1f' % float(eigval[0][i])
    
    print >>fo,'\n    ... Have a G00D Day, Week, Month, Year, and Century (if you are lucky) ...    '\
               '\n               Bye-Bye! Tschuess! Ciao! Poka! Zia Jian! KhodaHafez!             \n'
    if (fo):fo.close()
    #--------------------------------------------------------------------------------------------------------------------------------
    
head = {                                                                     \
    'CI':'\
        for, space group-number between 207 and 230, cubic structure.        \n\n\
                   C11     C12     C12      0       0       0                  \n\
                   C12     C11     C12      0       0       0                  \n\
                   C12     C12     C11      0       0       0                  \n\
                    0       0       0      C44      0       0                  \n\
                    0       0       0       0      C44      0                  \n\
                    0       0       0       0       0      C44                 \n',\
    'CII':'\
        for, space group-number between 195 and 206, cubic structure.        \n\n\
                   C11     C12     C12      0       0       0                  \n\
                   C12     C11     C12      0       0       0                  \n\
                   C12     C12     C11      0       0       0                  \n\
                    0       0       0      C44      0       0                  \n\
                    0       0       0       0      C44      0                  \n\
                    0       0       0       0       0      C44                 \n',\
    'HI':'\
        for, space group-number between 177 and 194, hexagonal structure.    \n\n\
                   C11     C12     C13      0       0       0                  \n\
                   C12     C11     C13      0       0       0                  \n\
                   C13     C13     C33      0       0       0                  \n\
                    0       0       0      C44      0       0                  \n\
                    0       0       0       0      C44      0                  \n\
                    0       0       0       0       0   (C11-C12)/2            \n',\
    'HII':'\
        for, space group-number between 168 and 176, hexagonal structure.    \n\n\
                   C11     C12     C13      0       0       0                  \n\
                   C12     C11     C13      0       0       0                  \n\
                   C13     C13     C33      0       0       0                  \n\
                    0       0       0      C44      0       0                  \n\
                    0       0       0       0      C44      0                  \n\
                    0       0       0       0       0   (C11-C12)/2            \n',\
    'RI':'\
        for, space group-number between 149 and 167, rhombohedral structure. \n\n\
                   C11     C12     C13     C14      0       0                  \n\
                   C12     C11     C13    -C14      0       0                  \n\
                   C13     C13     C33      0       0       0                  \n\
                   C14    -C14      0      C44      0       0                  \n\
                    0       0       0       0      C44     C14                 \n\
                    0       0       0       0      C14  (C11-C12)/2            \n',\
    'RII':'\
        for, space group-number between 143 and 148, rhombohedral structure. \n\n\
                   C11     C12     C13     C14     C15      0                  \n\
                   C12     C11     C13    -C14    -C15      0                  \n\
                   C13     C13     C33      0       0       0                  \n\
                   C14    -C14      0      C44      0     -C15                 \n\
                   C15    -C15      0       0      C44     C14                 \n\
                    0       0       0     -C15     C14  (C11-C12)/2            \n',\
    'TI':'\
        for, space group-number between 89 and 142, tetragonal structure.    \n\n\
                   C11     C12     C13      0       0       0                  \n\
                   C12     C11     C13      0       0       0                  \n\
                   C13     C13     C33      0       0       0                  \n\
                    0       0       0      C44      0       0                  \n\
                    0       0       0       0      C44      0                  \n\
                    0       0       0       0       0      C66                 \n',\
    'TII':'\
        for, space group-number between 75 and 88, tetragonal structure.     \n\n\
                   C11     C12     C13      0       0      C16                 \n\
                   C12     C11     C13      0       0     -C16                 \n\
                   C13     C13     C33      0       0       0                  \n\
                    0       0       0      C44      0       0                  \n\
                    0       0       0       0      C44      0                  \n\
                   C16    -C16      0       0       0      C66                 \n',\
    'O':'\
        for, space group-number between 16 and 74, orthorhombic structure.   \n\n\
                   C11     C12     C13      0       0       0                  \n\
                   C12     C22     C23      0       0       0                  \n\
                   C13     C23     C33      0       0       0                  \n\
                    0       0       0      C44      0       0                  \n\
                    0       0       0       0      C55      0                  \n\
                    0       0       0       0       0      C66                 \n',\
    'M':'\
        for, space group-number between 3 and 15, monoclinic structure.      \n\n\
                   C11     C12     C13      0       0      C16                 \n\
                   C12     C22     C23      0       0      C26                 \n\
                   C13     C23     C33      0       0      C36                 \n\
                    0       0       0      C44     C45      0                  \n\
                    0       0       0      C45     C55      0                  \n\
                   C16     C26     C36      0       0      C66                 \n',\
    'N':'\
        for, space group-number between 1 and 2, triclinic structure.        \n\n\
                   C11     C12     C13     C14      C15    C16                 \n\
                   C12     C22     C23     C24      C25    C26                 \n\
                   C13     C23     C33     C34      C35    C36                 \n\
                   C14     C24     C34     C44      C45    C46                 \n\
                   C15     C25     C35     C45      C55    C56                 \n\
                   C16     C26     C36     C46      C56    C66                 \n'}
    #--------------------------------------------------------------------------------------------------------------------------------

