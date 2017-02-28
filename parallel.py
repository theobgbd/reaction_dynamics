# -*- coding: utf-8 -*-
"""
Théo Beigbeder
Dynamique de la réaction AB+C 
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import arange, cos, exp
from scipy.integrate import odeint
from math import *
import random
import os,sys
import time
from mpl_toolkits.mplot3d import *
from matplotlib import cm
from decimal import *
from multiprocessing import Pool
from multiprocessing import Process

import routineABC as routine
import potABC as pot

global isurf
global ma
global mb
global mc
global Evib
global Ecoll
global RABin
global ncoup
global dt
global tsup
global eps
global isurf
global DU
global itraj
global iprint

#Mute output
save_stdout = sys.stdout
sys.stdout = open('trash', 'w')

#Fichiers de sortie
fileEvib = open("EvibProd.dat", "w+")
fileRAB = open("RAB.dat", "w+")
fileEvibreac=open("EvibReac.dat","w+")
filestat=open("statsparall.dat","w")
#filestat.close()
#fileinitial=open("initial.dat","w+")
filetraj=open("traj.dat","w+")
filevibphase=open("vibphase.dat","w+")
#fileweighttraj=open("weighttraj.dat","w+")
#filewheightpes=open("wheightpes.dat","w+")
#fileresults=open("Results.dat","w+")
#filePES=open("PES.dat","w+")


#=======Declaration des tableaux===================
R=[0,0,0]
y=[0,0,0,0]
RC=[0,0,0]
RT=[0,0,0]
RP=[0,0,0]
traj=[[],[],[]]

#===========Lecture du fichier d'entrée==========================
infile=open("./DataReacABC.dat",'r')
ma=float(infile.readline().split()[0])
mb=float(infile.readline().split()[0])
mc=float(infile.readline().split()[0])
Evib=float(infile.readline().split()[0])
[isample,rminus,rmaxus]=(infile.readline().split())
Ecollmax=float(infile.readline().split()[0])
RABin=float(infile.readline().split()[0])
ncoup=int(infile.readline().split()[0])
dt=float(infile.readline().split()[0])
tsup=float(infile.readline().split()[0])
eps=float(infile.readline().split()[0])
isurf=int(infile.readline().split()[0])
DU=float(infile.readline().split()[0])
itraj=float(infile.readline().split()[0])
iprint=float(infile.readline().split()[0])
plotPES=float(infile.readline().split()[0])
v=float(infile.readline().split()[0])

infile.close()
#Import des potentiels
[SAB,DAB,RAB,aAB,SAC,DAC,RAC,aAC,SBC,DBC,RBC,aBC]=pot.LectPotABC(isurf)

#Facteurs
uma_ua=1836
Ang_ua=1.89
ua_Ang=0.529
eV_ua=0.0367512
ua_eV=27.21
timesec_ua=4.13411E4

#Switch
ma=ma*uma_ua
mb=mb*uma_ua
mc=mc*uma_ua
mrbc=(mb*mc)/(mb+mc)
mrab=(mb*ma)/(mb+ma)
mrabc=(ma*(mb+mc))/(ma+mb+mc)
xk=2*DBC*(aBC**2)
Evib=Evib*eV_ua
Ecollmax=Ecollmax*eV_ua
eps=eps*eV_ua
DU=DU*eV_ua
tsup=tsup*timesec_ua
dt=dt*timesec_ua
isample=int(isample)
Etot=0

debut=0
fin=0
getcontext().prec = 5
dEcoll=(0.001)
fin=0
fin_i=0
#RP=pot.plotPES(isurf) 

#print(Evib,Ecollmax)

#Compteurs
nrxn=0
nrflx=0
nexpl=0
Evibmean=0
lEcoll=[]
lpreac=[]
#Evib=6.64E-34*(v+1/2)

t0 = time.time()
n=0
Ecoll=0.0009

#=============Boucle principale (trajectoires)==================== 
#Ecoll=Ecollmax-dEcoll*eV_ua
def simul(Ecoll):
    R=[0,0,0]
    y=[0,0,0,0]
    RC=[0,0,0]
    RT=[0,0,0]
    RP=[0,0,0]
    traj=[[],[],[]]
    nrxn=0
    nrflx=0
    nexpl=0
    Evibmean=0
    lEcoll=[]
    lpreac=[]
    nrxn=0
    nrflx=0
    nexpl=0
    Evibmean=0
    print('#################################################################################################')
    print('Ecoll=',Ecoll*ua_eV,"eV")
    for ntraj in range(0,ncoup):
        print('===================================')
        print("Trajectoire",ntraj)
        debut = time.time()
        xk=2*DBC*(aBC**2)
        dr=sqrt(2*mrbc*Evib)
        
        print('Ecoll i=',Ecoll*ua_eV,"eV")
        print("Evibi:", Evib*ua_eV)
        
        prmax=sqrt(2*mrbc*Evib)
        theta=random.random()*pi
        
        if isample==0:
            Rmax=RBC-(1/aBC)*log(1-sqrt(Evib/DBC))
            Rmin=RBC-(1/aBC)*log(1+sqrt(Evib/DBC))
        else:
            Rmax=float(rmaxus)
            Rmin=float(rminus)
        R[0]=RABin*Ang_ua
        
        if Evib==0:
            y[0]=RBC
            y[1]=0
            y[2]=R[0]+((mc*RBC)/(mb+mc))
            R[1]=y[0]
            R[2]=R[0]+R[1]
        else:
            y[0]=Rmin+float(ntraj)/float(ncoup+1)*(Rmax-Rmin)
            y[2]=R[0]+((mc*RBC)/(mb+mc))
            RC=routine.coord(y[2],y[0])
            Vleps=pot.PotABC(RC)
            b=random.random()
            print(Evib,Vleps,DBC)
            if (b<0.5):
                y[1]=sqrt((2*mrbc*(Evib-Vleps-DBC)))
            else:
                y[1]=-sqrt((2*mrbc*(Evib-Vleps-DBC)))
            
            R[1]=y[0]
            R[2]=R[0]+R[1]
            
        y[3]=-sqrt(2*Ecoll*mrabc)
        yin1=y[0]
        yin2=y[1]    
        #========= Integration par RK4 =========
        t=0
        #========= Boucle des temps ============
        print(tsup)
        while (t<tsup):
            RC=routine.coord(y[2],y[0])    
            Vleps=pot.PotABC(RC)
            Etotp=(y[1]**2)/(2*mrbc)+(y[3]**2)/(2*mrabc)+Vleps+DBC
            #---------------------------
            #print(y)
            [y,t]=routine.RK4(y,4,t,dt)
                #print(t,y)
    
            #---------------------------
            RC=routine.coord(y[2],y[0])    
            Vleps=pot.PotABC(RC)
            Etot=(y[1]**2)/(2*mrbc)+(y[3]**2)/(2*mrabc)+Vleps+DBC
            DeltaE=abs(Etotp-Etot)
                #print("\n",Etotp,Etot)
            if (DeltaE>eps):
                print("WARNING : Conservation of energy error !!!!","DE=",DeltaE*ua_eV," eV")
                break
            
            if ntraj==itraj:
                filevibphase.write(str(y[0])+str(" ")+str(y[1])+str(" ")+str(Vleps*ua_eV)+str("\n"))
                #print(str(y[0])+str(" ")+str(y[1])+str(" ")+str(Vleps*ua_eV)+str("\n"))
                filetraj.write(str(RC[0]*ua_Ang)+" "+str(RC[1]*ua_Ang)+" "+str(Vleps*ua_eV)+str("\n"))
                traj[0]=traj[0]+[RC[0]*ua_Ang]
                traj[1]=traj[1]+[RC[1]*ua_Ang]
                traj[2]=traj[2]+[Vleps*ua_eV]
                
                #fileweighttraj.write(str((y[0]*ua_Ang/mrbc))+str(" ")+str(y[2]*ua_Ang/mrabc)+str(" ")+str(Vleps*ua_eV)+str("\n"))
        
            if RC[1]>(4*RBC) :
                if RC[0]<4*RAB :
                    nrxn=nrxn+1
                    #print(ntraj,yin1,yin2)
                    if iprint==1:
                        print('*** REAC ***')
                RT=routine.coord(y[2],y[0])
                Vtest=pot.PotABC(RT)
                xmuab=(ma*mb)/(ma+mb)
                pab=xmuab/mrabc*y[3]-xmuab*(mc/(mc+mb))/mrbc*y[1]
                Evibtest=(pab**2/(2*xmuab))+Vtest+DAB
                Evibmean=Evibmean+Evibtest
                print('Evibmean: ',Evibmean/Etot*100)
                #fileinitial.write(yin1*ua_Ang,Evibtest*ua_eV)
                break
                
            if (RC[0])>(2*RABin*Ang_ua):
                print('*** REFl ***')
                nrflx=nrflx+1
                print(ntraj,yin1,yin2)
                RT=routine.coord(y[2],y[0])
                Vtest=pot.PotABC(RT)
                Evibtot=y[1]**2/(2*mrbc)+Vtest+DBC
                #print(yin1*ua_Ang,Evibtt*ua_eV)
                break    

        #fin=time.time()
        #dtime=fin-debut
        
       # print("Temps moyen de calcul/trajectoire:",round(dtime/60),'m',round(dtime%60),'s',"\nTemps restant estimé:",round((dtime*(ncoup-ntraj))/60),'minutes')
        #print (y)

    #==========Stats================#
    nINDET=ncoup-nrxn-nrflx
    prxn=nrxn/ncoup
    prflx=nrflx/ncoup
    pnind=nINDET/ncoup
    PEvibmean=(Evibmean/float(nrxn+1))
    Epc=((PEvibmean/(Etot))*100)
    dP=sqrt(prxn*(1-prxn)/ncoup)    
    lEcoll=lEcoll+[Ecoll]
    lpreac=lpreac+[prxn]
    filestat=open("statsparall.dat","a")
    filestat.write(str(Ecoll*ua_eV)+str(" ")+str(prxn)+str(" ")+str(dP)+str("\n"))
    fileRAB=open("RABparall.dat","a")
    fileRAB.write(str(yin2*ua_Ang)+str(" ")+str(prxn)+str("\n"))
    print("===============================================================================")        
    print("Evib=",Evib*ua_eV,'PREAC=',prxn,'PREF=',prflx,'PATHOS=',pnind)
    print("===============================================================================")        
    print('PREAC=',prxn*100,"%")
    print('PREF=',prflx*100,"%")
    print('PATHOS=',pnind*100,"%")
    print('Evib-final=',PEvibmean*ua_eV,"eV")
    print('Etot:',Etot*ua_eV)
    print('Bilan energetique:',Epc,"%")
    print("===============================================================================")        

    #frequences=frequences+[Ecoll*ua_eV,prxn]
    
    #print(filestat.write("hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh"))
    print(str(Ecoll*ua_eV)+str(" ")+str(prxn)+str(" ")+str(dP)+str("\n"))
    return [prxn,yin2,traj]
    
#print (frequences)
#plt.plot(lEcoll,lpreac)
stats=[]
with Pool(6) as p:
        print (p)
        stats=(p.map(simul, np.linspace(0.001*eV_ua,2*eV_ua,100)))
#print(stats)

sys.stdout = save_stdout  
os.system("gnuplot -p plotparall.plt")
if plotPES==1:
    #RP=pot.plotPES(isurf)
    RP=pot.plotPEScontour(isurf,stats[2])

    