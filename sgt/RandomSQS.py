#!/usr/bin/env python
# -*- coding: utf8 -*-
from ase.lattice.cubic import *
 
from numpy import *
from math import *
from ase import *
from pylab import *
import random
import pickle
import inspect
from ase.io.trajectory import  PickleTrajectory
#print '**** Create the 5x9 supercell structure ****'
#LC      = input('1. What is the lattice constant in Angstron? \n')
#NB    = input('2. What is the length of dislocation line in Bergers vector? \n')
#E1      = raw_input('3. What is the first element symbol? \n')
#E2      = raw_input('4. What is the second element symbol? \n')
#x    = input('5. What is the concentration of second alloying element? i.e. 0.25 \n')

class SQSscrew(object):
    r"""
         Create the object that contains the 5x9 supercell structure with screw dislocations
         
         LC: number
                Lattice constant in Ångström.
         NB: number   
                Length of dislocation line in Bergers vector.
         E1: string
                First element symbol.
         E2: string   
                Second element symbol.
         x: number   
                Concentration of second alloying element. i.e. 0.25.
    """
    
    def   __init__(self,LC=None,NB=None,E1=None,E2=None,x=.2,trials=10000 ):
        
        
        self.NB=NB
        self.concentration=x
        self.E1=E1
        self.E2=E2
        self.lat=Atoms()
        ### single slice of 135-atom unitcell
        SS      = BodyCenteredCubic(directions=[[1,1,-2],[-2,7,-5],[1,1,1]],size=(5,1,1),symbol=E1,latticeconstant=LC)
        cell    = SS.get_cell()
        ### define the two dislocation centers
        c    = [[cell[0][0]/2,cell[1][1]*(1/3.+2/27.),0],[cell[0][0],cell[1][1]*(1/3.+2/27.-1/27.),0]]
        ### NB-slice supercell
        self.SC    = SS.repeat((1,1,NB))
        print 'Total number of atoms is:', len(self.SC)
        ### ideal pair-wise correlation coefficient: c0
        c0    = (2*x-1)**2
        print 'The ideal pair-wise correlation coefficient is:', c0
        
        
        print '**** Get the nearest neighbor list: NN1 to NN6 ****'
        dist = [round(sqrt(3)/2*LC,6), round(LC,6), round(sqrt(2)*LC,6), round(sqrt(11)/2*LC,6), round(sqrt(3)*LC,6), round(2*LC,6)]
        NN1 = [[] for i in range(len(self.SC))]
        NN2 = [[] for i in range(len(self.SC))]
        NN3 = [[] for i in range(len(self.SC))]
        NN4 = [[] for i in range(len(self.SC))]
        NN5 = [[] for i in range(len(self.SC))]
        NN6 = [[] for i in range(len(self.SC))]
        ### repeat the supercell structure to include the neighbors of the bondary atoms
        NSC = self.SC.repeat((3,3,3))
        for i in range(13*len(self.SC),14*len(self.SC)):
            pos1 = NSC[i].get_position()
            for j in range(len(NSC)):
                pos2 = NSC[j].get_position()
                diff  = pos1 - pos2
                bond  = sqrt(dot(diff,diff))
                if   round(bond,6)==dist[0]:
                    NN1[i-13*len(self.SC)].append(j)
                elif round(bond,6)==dist[1]:
                    NN2[i-13*len(self.SC)].append(j)
                elif round(bond,6)==dist[2]:
                    NN3[i-13*len(self.SC)].append(j)
                elif round(bond,6)==dist[3]:
                    NN4[i-13*len(self.SC)].append(j)
                elif round(bond,6)==dist[4]:
                    NN5[i-13*len(self.SC)].append(j)
                elif round(bond,6)==dist[5]:
                    NN6[i-13*len(self.SC)].append(j)
        
        
        print '**** Insert two dislocations in single slice ****'
        Nump    = 20
        dens    = 0.4
        s0,r0    = 0,0
        u    = [[[(r-r0)/dens,(s-s0)/dens] for s in range(1,Nump+1)] for r in range(1,Nump+1)]
        Num    = 10
        rot    = atan((c[1][1]-c[0][1])/(c[1][0]-c[0][0]))
        ucen    = [[[[c[g][0]+cell[0][0]*j+cell[1][0]*i,c[g][1]+cell[0][1]*j+cell[1][1]*i,0] for g in range(len(c))] for i in range(-Num,Num+1)] for j in range(-Num,Num+1)]
        L    = SS.get_positions()
        phi    = [[sum([sum([sum([(-1)**g*atan2(-sin(rot)*(-u[r][s][0]+ucen[j][i][g][0])+cos(rot)*(-u[r][s][1]+ucen[j][i][g][1]),cos(rot)*(-u[r][s][0]+ucen[j][i][g][0])+sin(rot)*(-u[r][s][1]+ucen[j][i][g][1])) for g in range(len(c))]) for i in range(2*Num+1)]) for j in range(2*Num+1)]) for s in range(Nump)] for r in range(Nump)]
        Lphi    = [sum([sum([sum([(-1)**g*atan2(-sin(rot)*(-L[r][0]+ucen[j][i][g][0])+cos(rot)*(-L[r][1]+ucen[j][i][g][1]),cos(rot)*(-L[r][0]+ucen[j][i][g][0])+sin(rot)*(-L[r][1]+ucen[j][i][g][1])) for g in range(len(c))]) for i in range(2*Num+1)]) for j in range(2*Num+1)]) for r in range(len(L))]
        Lzero    = sum([sum([sum([(-1)**g*atan2(-sin(rot)*ucen[j][i][g][0]+cos(rot)*ucen[j][i][g][1],cos(rot)*ucen[j][i][g][0]+sin(rot)*ucen[j][i][g][1]) for g in range(len(c))]) for i in range(2*Num+1)]) for j in range(2*Num+1)])
        Lsp1    = sum([sum([sum([(-1)**g*atan2(-sin(rot)*(-cell[0][0]+ucen[j][i][g][0])+cos(rot)*(-cell[0][1]+ucen[j][i][g][1]),cos(rot)*(-cell[0][0]+ucen[j][i][g][0])+sin(rot)*(-cell[0][1]+ucen[j][i][g][1])) for g in range(len(c))]) for i in range(2*Num+1)]) for j in range(2*Num+1)])
        Lsp2    = sum([sum([sum([(-1)**g*atan2(-sin(rot)*(-cell[1][0]+ucen[j][i][g][0])+cos(rot)*(-cell[1][1]+ucen[j][i][g][1]),cos(rot)*(-cell[1][0]+ucen[j][i][g][0])+sin(rot)*(-cell[1][1]+ucen[j][i][g][1])) for g in range(len(c))]) for i in range(2*Num+1)]) for j in range(2*Num+1)])
        T    = (mat(cell).T).I
        LLLLint    =[dot(T,L[i]) for i in range(len(L))]
        Lcorrect=[LLLLint[i][0,0]*(Lsp1-Lzero)+LLLLint[i][0,1]*(Lsp2-Lzero) for i in range(len(LLLLint))]
        Lscrew    =[L[i]+[0,0,(Lphi[i]-Lcorrect[i])*cell[2][2]/2/pi] for i in range(len(L))]
        self.screw    = SS.copy()
        self.screw.set_positions(Lscrew)
        
        
        ### get the nearest and next nearest neighbor atom indexs for both dislocations 
        cnn1, cnnn1, cnn2, cnnn2 = [], [], [], []
        for i, p in enumerate(L):
            diff1 = p[:2]-c[0][:2]
            diff2 = p[:2]-c[1][:2]
            r1 = sqrt(dot(diff1, diff1))
            r2 = sqrt(dot(diff2, diff2))
            if r1 < dist[0]:
                cnn1.append(i)
            elif dist[0] < r1 < dist[1]: 
                cnnn1.append(i)
            if r2 < dist[0]:
                cnn2.append(i)
            elif dist[0] < r2 < dist[1]:
                cnnn2.append(i)
        
        ### apply the polarity to obtain the asymmetric core
        cpos = self.screw.get_positions()
        for pl in [1.0]:
            cposn = cpos.copy()
            for i in cnn1:
                cposn[i] = cposn[i]+pl*cell[2]*2/18
            for i in cnnn1:
                cposn[i] = cposn[i]-pl*cell[2]*1/18
            for i in cnn2:
                cposn[i] = cposn[i]-pl*cell[2]*2/18
            for i in cnnn2:
                cposn[i] = cposn[i]+pl*cell[2]*1/18
            self.screwn = self.screw.copy()
            self.screwn.set_positions(cposn)
        
        
        print '**** Generate random alloys and calculate the correlation coefficients ****'
        cc,ll,rmse,l=[],[],1,[]
        
        traj=PickleTrajectory('Random.traj','w')
        for no in range(trials):
            index= []
            ln = range(len(self.SC))
            ### replace E1 atoms randomly
            for i in range(int(len(self.SC)*x)):
                j = random.randint(i,len(self.SC)-1)
                (ln[i],ln[j])=(ln[j],ln[i])
                index.append(ln[i])
            AB = Atoms()
            AA = self.SC.copy()
            AA.set_tags([i+1 for i in range(len(AA))])
            for i in index:
                AA[i].set_symbol(E2)
                AB.append(AA[i])
            ### specify the spin value: 1 for E1 and -1 for E2 
            NA = AA.repeat((3,3,3))
            s = [0]*len(NA)
            for i in range(len(NA)):
                symbol=NA[i].get_symbol()
                if symbol==E1:
                    s[i]=1.
                else:
                    s[i]=-1.
            ### sum up the products of pair-wise
            c1,c2,c3,c4,c5,c6=0,0,0,0,0,0
            for i in range(13*len(self.SC),14*len(self.SC)):
                c1=c1+sum([s[i]*s[j] for j in NN1[i-13*len(self.SC)]])/len(NN1[i-13*len(self.SC)])/len(self.SC)
                c2=c2+sum([s[i]*s[j] for j in NN2[i-13*len(self.SC)]])/len(NN2[i-13*len(self.SC)])/len(self.SC)
                c3=c3+sum([s[i]*s[j] for j in NN3[i-13*len(self.SC)]])/len(NN3[i-13*len(self.SC)])/len(self.SC)
                c4=c4+sum([s[i]*s[j] for j in NN4[i-13*len(self.SC)]])/len(NN4[i-13*len(self.SC)])/len(self.SC)
                c5=c5+sum([s[i]*s[j] for j in NN5[i-13*len(self.SC)]])/len(NN5[i-13*len(self.SC)])/len(self.SC)
                c6=c6+sum([s[i]*s[j] for j in NN6[i-13*len(self.SC)]])/len(NN6[i-13*len(self.SC)])/len(self.SC)
            c1,c2,c3,c4,c5,c6=round(c1,4),round(c2,4),round(c3,4),round(c4,4),round(c5,4),round(c6,4)
            cn = [c1,c2,c3,c4,c5,c6]
            cc.append(cn)
            ### standart deviation 
            rmsen = sqrt(sum([(cn[i]-c0)**2 for i in range(len(cn))])/len(cn))
            print no, c1, c2, c3, c4, c5, c6, rmsen
            ### put E2 atoms to the end of the list for VASP
            del(AA[index])
            for atom in AB:
                AA.append(atom)
            traj.write(AA)
            ### get the atomic order for new random alloy
            order = AA.get_tags()
            dictvar  = [[i+1,order[i]] for i in range(len(order))]
            for i in range(len(dictvar)):
                (dictvar[i][0],dictvar[i][1])=(dictvar[i][1],dictvar[i][0])
            dictvar.sort()
            norder = [dictvar[i][1] for i in range(len(dictvar))]
            ll.append(norder)
            ### compare the new structure with the old one
            if rmsen<rmse:
                self.rmse = rmsen
                self.cp   = cn
                self.l    = norder
                self.lat  = AA
            
 
        
        out = file('Random.dat','a')
        print self.cp[0],self.cp[1],self.cp[2],self.cp[3],self.cp[4],self.cp[5], self.rmse, '==> final result!'
        print >>out, self.cp[0],self.cp[1],self.cp[2],self.cp[3],self.cp[4],self.cp[5], self.rmse
        for i in self.l:
            print >>out, i
        dictvar = [[i,self.l[i]] for i in range(len(self.l))]
        for i in range(len(dictvar)):
            (dictvar[i][0],dictvar[i][1])=(dictvar[i][1],dictvar[i][0])
        dictvar.sort()
        norder = [dictvar[i][1] for i in range(len(dictvar))]
        SC0 = self.screw.repeat((1,1,self.NB))
        SC1 = self.screwn.repeat((1,1,self.NB))
        self.lat0 = Atoms()
        self.lat1 = Atoms()
        self.lat0.set_cell(self.SC.get_cell())
        self.lat1.set_cell(self.SC.get_cell())
        for i in norder:
            self.lat0.append(SC0[i])
            self.lat1.append(SC1[i])
        for i in range(int(len(self.SC)*(1-self.concentration)), len(self.SC)):
            self.lat0[i].set_symbol(self.E2)
            self.lat1[i].set_symbol(self.E2)
        tcell = self.SC.get_cell()
        tcell[1][2] = -tcell[2][2]/self.NB/2
        self.lat0.set_cell(tcell,scale_atoms=True)
        self.lat1.set_cell(tcell,scale_atoms=True)
        #write('BCC.poscar',self.lat,format='vasp',direct=True)
        #write('SYM.poscar',self.lat0,format='vasp',direct=True)
        #write('ASY.poscar',self.lat1,format='vasp',direct=True)
        print '############# Done! ##############' 
                
    def show(self, file=None):
        """
        The show method displays the screw dislocations.
         
        file: string
            If given, the picture is written to a file. The format is determined by the suffix.
            Else the picture is displayed to the user.
        
        """
    
        cell  = self.lat.get_cell()
        pos_i = self.lat.get_positions()
        pos_f = self.lat0.get_positions()
        pos_d = pos_f - pos_i
        
        nn = []
        for i1 in range(len(pos_i)):
            for i2 in range(i1,len(pos_i)):
                diff = pos_i[i1][:2] - pos_i[i2][:2]
                r = sqrt(np.dot(diff, diff))
                if 0.1 < r and r < 2.8: nn.append([i1,i2])
        
        X0= [p[0] for p in pos_i]
        Y0= [p[1] for p in pos_i]
        X = [(pos_i[i[0]][0]+pos_i[i[1]][0])/2 for i in nn]
        Y = [(pos_i[i[0]][1]+pos_i[i[1]][1])/2 for i in nn]
        
        d = [pos_d[nn[i][0]][2]-pos_d[nn[i][1]][2] for i in range(len(nn))]
        for i in range(len(d)):
            while True:
                try:
                    if d[i]/cell[2][2] > 1./2:
                        d[i] = d[i]-cell[2][2]
                    elif d[i]/cell[2][2] < -1./2:
                        d[i] = d[i]+cell[2][2]
                    else:
                        break
                except:
                    print "Retrying!"
        
        U0= [pos_i[nn[i][0]][0]-pos_i[nn[i][1]][0] for i in range(len(nn))]
        V0= [pos_i[nn[i][0]][1]-pos_i[nn[i][1]][1] for i in range(len(nn))]
        U = [(pos_i[nn[i][0]][0]-pos_i[nn[i][1]][0])*d[i]/cell[2][2]*3 for i in range(len(nn))]
        V = [(pos_i[nn[i][0]][1]-pos_i[nn[i][1]][1])*d[i]/cell[2][2]*3 for i in range(len(nn))]
        
        dc1, dc2 = [cell[0][0]/2,cell[1][1]*(1/3.+2/27.)], [cell[0][0],cell[1][1]*(1/3.+1/27.)]
        symbols = self.lat.get_chemical_symbols()
        colors  = ['w']*len(self.lat)
        for i, symbol in enumerate(symbols):
            if symbol != 'W': colors[i]='g'
        
        scatter(X0,Y0,s=80,marker='o',c=colors,lw=2,label='W atom')
        scatter(dc1[0],dc1[1],s=30,marker='^',c='m',label='disl cent 1')
        scatter(dc2[0],dc2[1],s=30,marker='^',c='m',label='disl cent 2')
        plot((0,cell[0][0],cell[0][0]+cell[1][0],cell[1][0],0),(0,0,cell[1][1],cell[1][1],0),'g-',lw=2,label='unitcell')
        ax = gca()
        
        com_tuple = []
        for i in range(len(nn)):
            if 2*cell[1][1]/27 < sqrt((X[i]-dc1[0])**2+(Y[i]-dc1[1])**2) < cell[0][0]/15 or 2*cell[1][1]/27 < sqrt((X[i]-dc2[0])**2+(Y[i]-dc2[1])**2) < cell[0][0]/15:
                com_tuple.append((X[i],abs(d[i])))
                color = 'r'
            else:
                color = 'k'
            arr = Arrow(X[i]-U[i]/2,Y[i]-V[i]/2,U[i],V[i],width=d[i],color=color)
            ax.add_patch(arr)
        
        com_sort = sorted(com_tuple)
        polarity = (com_sort[1][1]-com_sort[4][1]+com_sort[5][1]-com_sort[3][1]+com_sort[2][1]-com_sort[0][1])/cell[2][2]
        print 'Core polarity p = ', polarity
        #arr.set_label('diff disp')
        legend()
        xlabel('[11-2] --->')
        ylabel('[110] --->')
        title('1/2<111> screw dislocation')
        axis([-2.5, cell[0][0]+cell[1][0]+2.5, -(((cell[0][0]+cell[1][0])+5)*3/4-cell[1][1])/2, cell[1][1]+(((cell[0][0]+cell[1][0])+5)*3/4-cell[1][1])/2])
        if file :
            savefig(file)
        else:
            show()
    
    def get_bcc(self):
        """
        Returns the undisturbed lattice.
        """
        return self.lat
    def get_sym(self):
        """
        Returns the lattice with screw dislocation.
        """
        return self.lat0
    def get_asy(self):
        """
        Returns the lattice with asymmetric screw dislocation.
        """
        return self.lat1
    
    def save(self,filename):
        out=open(filename,"wb")
        pickle.dump(self,out,2)
        out.close()
      
