# -*- coding: utf-8 -*-
"""
PES -- 3 Atoms Reaction
"""
from math import *
from __main__ import *
import routineABC as routine

def LectPotABC(surf):
    from __main__ import DU
    if surf==1:
        SAB=0.05
        DAB=0.17445
        RAB=1.4
        aAB=1.02730
        SAC=0.05
        DAC=0.12660
        RAC=1.4
        aAC=1.02730
        SBC=0.3
        DBC=0.17445
        RBC=1.4
        aBC=1.02730
        
    elif surf==2:
        SAB=0.30
        DAB=0.17445
        RAB=1.40
        aAB=1.02730
        SAC=0.050
        DAC=0.1266
        RAC=1.40
        aAC=1.02730
        SBC=0.050
        DBC=0.17445
        RBC=1.40
        aBC=1.02730
    
    elif surf==3:
        SAB=0
        DAB=0.17445
        RAB=1.4
        aAB=1.3174
        SAC=0
        DAC=0
        RAC=0
        aAC=0
        SBC=0
        DBC=DU+DAB
        RBC=1.4
        aBC=1.3174
        
    return [SAB,DAB,RAB,aAB,SAC,DAC,RAC,aAC,SBC,DBC,RBC,aBC] 
    
#==========================Potentiel============================
def PotABC(R):
    from __main__ import isurf
    if (isurf==1) or (isurf==2):
        Vpot=LEPS(R)
        return Vpot
    elif (isurf==3):
        Vpot=DMorse(R)
        return Vpot
        
        
def DpotABC(Rg,Rp):
    from __main__ import isurf
    if(isurf==1) or (isurf==2):        
        return DPotLEPS(Rg,Rp) #[dVdRg,dvdRp] 
    elif (isurf==3):
        return DMorsederiv(Rg,Rp)

#===========================LEPS=================================

def LEPS(R): 
    from __main__ import isurf
    [SAB,DAB,RAB,aAB,SAC,DAC,RAC,aAC,SBC,DBC,RBC,aBC]=LectPotABC(isurf)

    #Delta R
    DRAB=R[0]-RAB
    DRBC=R[1]-RBC
    DRAC=R[2]-RAC
    
    #Alpha*Delta R
    aDRAB=-aAB*DRAB
    aDRAC=-aAC*DRAC
    aDRBC=-aBC*DRBC
    
    #Coeff
    CAB=DAB/(4*(1+SAB))
    CAC=DAC/(4*(1+SAC))
    CBC=DBC/(4*(1+SBC))

    #Integrales Q
    QAB=CAB*((3+SAB)*exp(2*aDRAB)-(2+6*SAB)*exp(aDRAB))    
    QAC=CAC*((3+SAC)*exp(2*aDRAC)-(2+6*SAC)*exp(aDRAC))
    QBC=CBC*((3+SBC)*exp(2*aDRBC)-(2+6*SBC)*exp(aDRBC))
    
    #Integrales J
    JAB=CAB*((1+3*SAB)*exp(2*aDRAB)-(6+2*SAB)*exp(aDRAB))
    JAC=CAC*((1+3*SAC)*exp(2*aDRAC)-(6+2*SAC)*exp(aDRAC))
    JBC=CBC*((1+3*SBC)*exp(2*aDRBC)-(6+2*SBC)*exp(aDRBC))
    
    #Expression du potentiel
    Vleps=QAB+QAC+QBC-sqrt(0.5*(((JAB-JBC)**2)+((JBC-JAC)**2)+(JAC-JAB)**2))
    return Vleps
    
def DPotLEPS(RG,RP):
    h=0.5E-6
#Dérivée de RG
    #Rg+h
    RGplus=RG+h
    R=routine.coord(RGplus,RP)   
    V1=LEPS(R)     
    #Rg-h
    RGmoins=RG-h
    R=routine.coord(RGmoins,RP)   
    V2=LEPS(R)
    #dV/dRG       
    dVdRG=(V1-V2)/(2*h)

#Dérivée de RP
    #Rg+h
    RPplus=RP+h
    R2=routine.coord(RG,RPplus)
    V3=LEPS(R2)     
    #Rg-h
    RPmoins=RP-h
    R2=routine.coord(RG,RPmoins)
    V4=LEPS(R2) 
    #dV/dRG       
    dVdRP=(V3-V4)/(2*h)
    
    return [dVdRG,dVdRP]
    
#==========================Morse=================================
def DMorse(R):
    from __main__ import isurf
    [SAB,DAB,RAB,aAB,SAC,DAC,RAC,aAC,SBC,DBC,RBC,aBC]=LectPotABC(isurf)
    VR1R2=DAB*(1-exp(-aAB*(R[0]-RAB)))**(2)-DAB+DBC*(1-exp(-aBC*(R[1]-RBC)))**2-DBC
    return VR1R2

def DMorsederiv(Rg,Rp):
    h=0.5E-6
    
    #====Derivee par Rg
    #Rg+h
    Rgplus=Rg+h
    R=routine.coord(Rgplus,Rp)
    Vleps=DMorse(R)
    V1=Vleps
    #Rg-h
    Rgmoins=Rg-h
    R=routine.coord(Rgmoins,Rp)
    Vleps=DMorse(R)
    V2=Vleps
    dVdRG=(V1-V2)/(2*h)
    
    #====Derivee par Rp
    #Rg+h
    Rpplus=Rp+h
    R=routine.coord(Rg,Rpplus)
    Vleps=DMorse(R)
    V3=Vleps
    #Rg-h
    Rpmoins=Rp-h
    R=routine.coord(Rg,Rpmoins)
    Vleps=DMorse(R)
    V4=Vleps
    dVdRP=(V3-V4)/(2*h)
   
    return [dVdRG,dVdRP]

#======================Plot de la PES ============================
def plotPES(isurf,**plotsettings):
    LectPotABC(isurf)
    
    if plotsettings.get("range") == None:
        xi=0
        xf=100
        yi=-5
        yf=100
    if plotsettings.get("depth") == None:
        depth=75
        
    X=np.linspace(xi,xf,depth)
    Y=np.linspace(yi,yf,depth)
    Px=[]
    Py=[]
    Pz=[]
    if(isurf==1) or (isurf==2):        
        for i in X:
            for j in Y:
                RP[0]=0.5+(i-1)/100*RABin*Ang_ua
                RP[1]=0.5+(j-1)/100*RBC*Ang_ua
                RP[2]=RP[0]+RP[1]
                Vpes=pot.PotABC(RP)
                Px=Px+[RP[0]*ua_Ang]
                Py=Py+[RP[1]*ua_Ang]
                Pz=Pz+[Vpes*ua_eV]
                filePES.write(str(RP[0]*ua_Ang)+str(" ")+str(RP[1]*ua_Ang)+str(" ")+str(Vpes*ua_eV)+"\n")
        

    elif (isurf==3):
        for i in X:
            for j in Y:
                RP[0]=0.5+(i-1)/100*RABin*Ang_ua
                RP[1]=0.5+(j-1)/100*RBC*Ang_ua
                RP[2]=RP[0]+RP[1]
                Vpes=pot.PotABC(RP)
                Px=Px+[RP[0]*ua_Ang]
                Py=Py+[RP[1]*ua_Ang]
                Pz=Pz+[Vpes*ua_eV]
                filewheightpes.write(str(RP[1]*ua_Ang/mrbc)+str(" ") + str ((RP[0]+((mc*RP[1])/(mb+mc)))*ua_Ang/mrabc)+str(" ") +str (Vpes*ua_eV)+"\n")
                filePES.write(str(RP[0]*ua_Ang)+str(" ")+str(RP[1]*ua_Ang)+str(" ")+str(Vpes*ua_eV)+"\n")
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #print(Px,Py,Pz)
    ax.plot_trisurf(Px, Py, Pz, cmap=cm.jet, linewidth=0.1)
    #ax.plot(Px,Py)
    plt.show()
    
    return R 


def plotPEScontour(isurf,traj,**plotsettings):
    uma_ua=1836
    Ang_ua=1.89
    ua_Ang=0.529
    eV_ua=0.0367512
    ua_eV=27.21
    timesec_ua=4.13411E4
    LectPotABC(isurf)
    filePES=open("PES.dat","w+")
    filewheightpes=open("wheightpes.dat","w+")
    from __main__ import Ecoll
    from __main__ import Evib 
    #from __main__ import Epc
    #from parallel import Epc
    from __main__ import RABin
    mpl.rcParams['text.usetex']=True
    mpl.rcParams['text.latex.unicode']=True

    if plotsettings.get("range") == None:
        xi=0
        xf=3.5
        yi=0
        yf=3.5
    if plotsettings.get("depth") == None:
        depth=60
        
    x=np.linspace(xi,xf,depth)
    y=np.linspace(yi,yf,depth)
    Px=[]
    Py=[]
    Pz=[]
    
    Pxt=[]
    Pyt=[]
    Pzt=[]
    X, Y = np.meshgrid(x, y)
    
    if(isurf==1) or (isurf==2):
        for i in x:
            Pz=[]
            for j in y:
                RP[0]=j*Ang_ua
                RP[1]=i*Ang_ua
                RP[2]=RP[0]+RP[1]
                Vpes=pot.PotABC(RP)
                Px=Px+[RP[0]]
                Py=Py+[RP[1]]
                Pz=Pz+[Vpes*ua_eV]
                filewheightpes.write(str(Px)+str(" ") + str (Py)+str(" ") +str (Pz)+"\n")
                filePES.write(str(RP[1]*ua_Ang)+str(" ")+str(RP[0]*ua_Ang)+str(" ")+str(Vpes*ua_eV)+"\n")
            Pxt=Pxt+[Px]
            Pyt=Pyt+[Py]
            Pzt=Pzt+[Pz]
        
        Z=np.array(Pzt)
        CS=plt.contour(X, Y, Z, np.arange(-4.5,-2,0.5),zdir='z',inline=0.5 , linewidth=0.1)
        plt.axis([0,3.5, 0, 3.5])
    
    elif (isurf==3):
        for i in x:
            Pz=[]
            for j in y:
                RP[0]=(i)*Ang_ua
                RP[1]=(j)*Ang_ua
                RP[2]=RP[0]+RP[1]
                Vpes=pot.PotABC(RP)
                Px=Px+[RP[0]]
                Py=Py+[RP[1]]
                Pz=Pz+[Vpes*ua_eV]
                filewheightpes.write(str(RP[1]*ua_Ang/mrbc)+str(" ") + str ((RP[0]+((mc*RP[1])/(mb+mc)))*ua_Ang/mrabc)+str(" ") +str (Vpes*ua_eV)+"\n")
                filePES.write(str(RP[1]*ua_Ang)+str(" ")+str(RP[0]*ua_Ang)+str(" ")+str(Vpes*ua_eV)+"\n")
            Pxt=Pxt+[Px]
            Pyt=Pyt+[Py]
            Pzt=Pzt+[Pz]
        Z=np.array(Pzt)   
        CS=plt.contour(X, Y, Z, np.arange(-5,0,0.5),zdir='z',inline=0.5 , linewidth=0.1)
        plt.axis([0,3.5, 0, 3.5])
        
    
    
       
    #Plot de la surface
    
    
    plt.plot(traj[0],traj[1], '-b',linewidth=2) 
    plt.xlabel(r'$r_1$ $(A)$')
    plt.ylabel(r'$r_2$ $(A)$')
    
    plt.clabel(CS, inline=1, fontsize=10)
    
    plt.text(2.5, 3, r'$E_{coll}^i=%s$ eV'%(np.round(Ecoll*ua_eV,3)), fontsize=18)
    plt.text(2.5, 2.75, r'$E_{vib}^i=%s$ eV'%(np.round(Evib*ua_eV,3)), fontsize=18)
    #plt.text(2.5, 2.5, r'$E_{vib}^f=%s$ '%(np.round(Epc,3)), fontsize=18)
    
    #Plot de la trajectoire
    
  
    plt.show()
    return R 
      
      

    #cset = ax.contourf(Px, Py, Pz,zdir='z',offset=-10,cmap=cm.winter)    # winter colour map
