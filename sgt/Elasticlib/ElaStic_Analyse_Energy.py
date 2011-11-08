#!/usr/bin/env python
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%% ------------------------------------------------ ElaStic_Analyse_Energy ----------------------------------------------- %%%#
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
# python ElaStic_Analyse_Energy.py
#        ElaStic_Analyse_Energy
# 
# EXPLANATION:
# 
#________________________________________________________________________________________________________________________________

from sys   import stdin
from numpy import *
from math import *
import numpy as np
import subprocess
import warnings
import os.path
import shutil
import math
import sys
import os
from gracepar import gracepar
import copy

#%%%%%%%%--- CONSTANTS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
_e     = 1.602176565e-19              # elementary charge
Bohr   = 5.291772086e-11              # a.u. to meter
Ryd2eV = 13.605698066                 # Ryd to eV
cnvrtr = (_e*Ryd2eV)/(1e9*Bohr**3)    # Ryd/[a.u.^3] to GPa
#--------------------------------------------------------------------------------------------------------------------------------

#%%%--- Subroutins ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def sortstrain(s,e):
    ss=[]
    ee=[]
    ww=[]
    for i in range(len(s)): ww.append(s[i])
    ww.sort()
    for i in range(len(s)):
        ss.append(s[s.index(ww[i])])
        ee.append(e[s.index(ww[i])])
    return ss, ee
#--------------------------------------------------------------------------------------------------------------------------------

def elastic_select_data(self):
    #%%%--- Reading the INFO_ElaStic file ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    ### ----------------------------- Calculating the second derivative and Cross-Validation error ------------------------------ ###
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    
    CONV = cnvrtr * factorial(self.order)*1.
    
    for i in range(1, self.ECs+1):
        if (i<10):
            Dstn = 'Dst0'+str(i)
        else:
            Dstn = 'Dst' +str(i)
    
        if (self.order == 2):
            fD = open(Dstn+'_d2E.dat', 'w')
        if (self.order == 3):
            fD = open(Dstn+'_d3E.dat', 'w')
    
        fE = open(Dstn+'_CVe.dat', 'w')
        print >> fD, '# Max. eta    SUM(Cij) \n#'
        print >> fE, '# Max. eta    Cross-Validation error   \n#'
    
        for j in range(self.order+4, self.order-1, -2):
            if  (j == 2): nth = '2nd'
            elif(j == 3): nth = '3rd'
            else:
                nth = str(j) + 'th'
    
            print >> fD, '\n# '+ nth +' order fit.'
            print >> fE, '\n# '+ nth +' order fit.'
    
    #------ Get Energies-------------------------------------------------------------------------------------------------
             
    
            nl     = 0
            strain = []
            energy = []
            while (nl < self.NPt):
                    
                strain.append(self.distortions[(i-1)*self.NPt+nl].eta[1])
                energy.append(self.distortions[(i-1)*self.NPt+nl].atoms.get_potential_energy()) 
                nl+=1
            strain, energy = sortstrain(strain,energy)
            straintmp, energytmp = copy.copy(strain),copy.copy(energy)
    #--------------------------------------------------------------------------------------------------------------------------------
            while (len(straintmp) > j): 
                emax  = max(straintmp)
                emin  = min(straintmp)
                emax  = max(abs(emin),abs(emax))
                coeffs= polyfit(straintmp, energytmp, j)
                if (self.order == 2):
                    Cij  = coeffs[j-2]*CONV/self.structure.get_volume()         # in GPa unit 
                if (self.order == 3):
                    Cij  = coeffs[j-3]*CONV/self.structure.get_volume() * 0.001 # in TPa unit
    
                print >>fD, '%13.10f'%emax, '%18.6f'%Cij
    
                if (abs(straintmp[0]+emax) < 1.e-7):
                    straintmp.pop(0); energytmp.pop(0)
                if (abs(straintmp[len(straintmp)-1]-emax) < 1.e-7):
                    straintmp.pop()
                    energytmp.pop()
          
    #------ Cross-Validation error calculations -------------------------------------------------------------------------------------
            straintmp, energytmp = copy.copy(strain),copy.copy(energy)
            while (len(straintmp) > j+1): 
                emax = max(straintmp)
                emin = min(straintmp)
                emax = max(abs(emin),abs(emax))
    
                S = 0
                for k in range(len(straintmp)):
                    Y      = energytmp[k]
                    etatmp = []
                    enetmp = []
    
                    for l in range(len(straintmp)):
                        if (l==k): pass
                        else:            
                            etatmp.append(straintmp[l])
                            enetmp.append(energytmp[l])
    
                    Yfit = polyval(polyfit(etatmp,enetmp, j), strain[k])
                    S    = S + (Yfit-Y)**2
    
                CV = sqrt(S/len(straintmp))
                print >>fE, '%13.10f'%emax, CV
    
                if (abs(straintmp[0]+emax) < 1.e-7):
                    straintmp.pop(0)
                    energytmp.pop(0)
                if (abs(straintmp[len(straintmp)-1]-emax) < 1.e-7):
                    straintmp.pop()
                    energytmp.pop()
               
        fD.close()
        fE.close()
        
    #-- Plotting --------------------------------------------------------------------------------------------------------------------
        
        Glines=gracepar.splitlines(True)
      
        TMP = []
        if (self.order == 2):
            for k in range(1, 45):
                TMP.append(Glines[k])
    
        if (self.order == 3):
            for k in range(48, 92):
                TMP.append(Glines[k])
    
        for k in range(164, 219):
            TMP.append(Glines[k])
    
        TMP.insert(99,'    s2 legend  " n = '+str(self.order+0)+'"\n')
        TMP.insert(91,'    s1 legend  " n = '+str(self.order+2)+'"\n')
        TMP.insert(83,'    s0 legend  " n = '+str(self.order+4)+'"\n')
        TMP.insert(46,'    subtitle "Plot for '+ Dstn +' deformation, n = Order of polynomial fit"\n')
    
        GdE = open(Dstn+'_d'+str(self.order)+'E.par', 'w')
        for l in range(len(TMP)):
            print >>GdE, TMP[l],
        GdE.close()
    
        os.system('xmgrace '+ Dstn +'_d'+str(self.order)+'E.dat -param '+Dstn+'_d'+str(self.order)+'E.par -saveall '+Dstn+'_d'+str(self.order)+'E.agr &')
    
        TMP = []
        for k in range(154, 162):
            TMP.append(Glines[k])
    
        for k in range(164, 219):
            TMP.append(Glines[k])
    
        TMP.insert(63,'    s2 legend  " n = '+str(self.order+0)+'"\n')
        TMP.insert(55,'    s1 legend  " n = '+str(self.order+2)+'"\n')
        TMP.insert(47,'    s0 legend  " n = '+str(self.order+4)+'"\n')
        TMP.insert(10,'    subtitle "Plot for '+ Dstn +' deformation, n = Order of polynomial fit"\n')
    
        CVe = open(Dstn+'_CVe.par', 'w')
        for k in range(len(TMP)):
            print >>CVe, TMP[k],
        CVe.close()
    
        os.system('xmgrace '+ Dstn +'_CVe.dat -param '+ Dstn +'_CVe.par -saveall '+ Dstn +'_CVe.agr &')
    
    #os.chdir('../')
    
    #--- Writing the ElaStic_???.in file --------------------------------------------------------------------------------------------
    if (self.order == 2): orth  = '2nd'
    if (self.order == 3): orth  = '3rd'
    
    fri = open('ElaStic_'+ orth +'.in', 'w')
    for i in range(1, self.ECs+1):
        if (i<10):
            Dstn = 'Dst0'+str(i)
        else:
            Dstn = 'Dst' +str(i)
        print >>fri, Dstn+'    eta_max    Fit_order'
    fri.close()
    #--------------------------------------------------------------------------------------------------------------------------------