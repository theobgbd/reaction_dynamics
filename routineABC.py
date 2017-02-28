# -*- coding: utf-8 -*-
"""
Fonctions routines pour l'etude de la dynamique de la réaction à trois atomes
"""
from __main__ import *
from potABC import *

def RK4(Y,N,T,H):
    nmax=4
    #DYDT=[0,0,0,0]
    DYDT = [None]*nmax
    YT=[None]*nmax
    DYT=[None]*nmax
    DYM=[None]*nmax
    
    if N>nmax:
        return
    
    DYDT=sysdiff(Y,T)

    HH=H*0.5
    H6=H/6
    TH=T+HH
    
    for i in range(0,N):
        YT[i]=Y[i]+HH*DYDT[i]
    
    DYT=sysdiff(YT,TH)
    
    for j in range(0,N):
        YT[j]=Y[j]+HH*DYT[j]
    
    DYM=sysdiff(YT,TH)
    
    for k in range(0,N):
        YT[k]=Y[k]+H*DYM[k]
        DYM[k]=DYT[k]+DYM[k]
    
    DYT=sysdiff(YT,T+H)
    
    for l in range(0,N):
        Y[l]=Y[l]+H6*(DYDT[l]+DYT[l]+2*DYM[l])
    T=T+H
    return [Y,T]

def sysdiff(Y,t):
    from __main__ import mrbc,mrabc
   # [ma,mb,mc,mrab,mrbc,mrabc]=mass
    dydt=[0,0,0,0]
    #Collision
    dydt[0]=Y[1]/mrbc
    dydt[2]=Y[3]/mrabc

    [dVdy3,dVdy1]=DpotABC(Y[2],Y[0])
    #Vibration
    dydt[1]=-dVdy1
    dydt[3]=-dVdy3
    #print(dydt)
    return dydt
        
    
def coord(RG,RP):
    from __main__ import mc,mb,ma
    R=[0,0,0]
    R[0]=RG-((mc*RP)/(mb+mc))
    R[1]=RP
    R[2]=R[0]+R[1]    
    return R

        
    
    
    
    
    
