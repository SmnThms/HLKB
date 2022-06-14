# -*- coding: utf-8 -*-
"""
Created on Mon Oct 2 9:00:00 2017, @author: Simon
"""

from ..config import LEVELS
from ..constant import S,I
from .Level import Level
from .wigner_functions import clebsch
import numpy as np

# included_levels = [(n,L)]

def LmSmLmI(included_levels=LEVELS):
    levels_list = []
    for (n,L) in included_levels:
        for mL in np.arange(-L,L+1):
            for mS in [-S,S]:
                for mI in np.arange(-I,I+1):
                    levels_list.append(Level(n=n,L=L,mS=mS,mL=mL,mI=mI,I=I,S=S))
    return levels_list
            
def LJmJmI(included_levels=LEVELS):
    levels_list = []
    for (n,L) in included_levels:
        for J in np.unique(np.abs(L+np.arange(-S,S+1))):
            for mJ in np.arange(-J,J+1):
                for mI in np.arange(-I,I+1):
                    levels_list.append(Level(n=n,L=L,J=J,mJ=mJ,mI=mI,I=I,S=S))
    return levels_list

def LJFmF(included_levels=LEVELS):
    levels_list = []
    for (n,L) in included_levels:
        for J in np.unique(np.abs(L+np.arange(-S,S+1))):
            for F in np.unique(np.abs(J+np.arange(-I,I+1))):
                for mF in np.arange(-F,F+1):
                    levels_list.append(Level(n=n,L=L,J=J,F=F,mF=mF,I=I,S=S))
    return levels_list

N = len(LJmJmI()) # dimension of our state space
# N1S = len([a for a in LJmJmI() if a.n==1]) # dimension of the 1S state space

def balmer_basis(): # to calculate the Balmer-series fluorescence
    levels_list = []
    for (n,L) in [(2,0),(2,1)]:
        for J in np.unique(np.abs([L+S,L-S])):
            for F in np.unique(np.abs(J+np.arange(-I,I+1))):
                for mF in np.arange(-F,F+1):
                    levels_list.append(Level(n=n,L=L,J=J,F=F,mF=mF,I=I,S=S))
    return levels_list # in the {|L,J,F,mF>} basis


# =============================================================================
# Change of basis
# =============================================================================

def convert(H,P): return np.dot(P,np.dot(H,P.T))

def LSI_to_LJI(): 
    P = np.zeros((N,N))
    for j, a in enumerate(LmSmLmI()): # initial level
        for i, b in enumerate(LJmJmI()): # final level
            if a.L!=b.L or a.mI!=b.mI or a.n!=b.n:
                P[i,j] = 0 # same mI and L required for both levels
            else:
                P[i,j] = clebsch(j1=a.L,m1=a.mL,j2=a.S,m2=a.mS,J=b.J,M=b.mJ)
    return P

def LJI_to_LSI():
    return LSI_to_LJI().transpose()
            
def LJI_to_LJF(): 
    P = np.zeros((N,N))
    for j, a in enumerate(LJmJmI()): # initial level
        for i, b in enumerate(LJFmF()): # final level
            if a.L!=b.L or a.J!=b.J or a.n!=b.n:
                P[i,j] = 0 # same J and L required for both levels
            else:
                P[i,j] = clebsch(j1=a.J,m1=a.mJ,j2=a.I,m2=a.mI,J=b.F,M=b.mF)
    return P

def LJF_to_LJI():
    return LJI_to_LJF().transpose()