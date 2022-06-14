# -*- coding: utf-8 -*-
"""
Created on Mon Oct 2 9:00:00 2017, @author: Simon
"""

from ..config import ISOTOPE 
from ..constant import h,c,alpha,rydberg,mu_en,an,I,S
from ..basis import LJFmF,wigner3j,wigner6j,wigner9j
from .hyperfine_structure_corrections import hfs_corrections
from mpmath import mpf, sqrt
import numpy as np

E0 = 2*h*c*rydberg/alpha**2         # Énergie propre de l'électron, J
En = E0/mu_en                       # Énergie propre du noyau, J
Er = E0/(1+mu_en)                   # Énergie propre réduite de l'électron, J

An = (1+an)*alpha**4*Er**3/E0/En/h # en Hz
N = len(LJFmF())


def H_hyperfin(): # Base {|LJFmF>}
    return An*2/3*M_delta_r_sur_r2(IS()) + An*M_sur_r3(IL()+IX())

def H_hyperfin_diagonal(): # Base {|LJFmF>}
    return np.diag([E_Fermi(a.n,a.L,a.J,a.F) for a in LJFmF()])

def H_hyperfin_non_diagonal(): # Base {|LJFmF>}
    M = M_sur_r3(IL())
    return An*3/2*(M - np.diag(np.diag(M)))

def E_Fermi(n,L,J,F):
    theta_nr = 1/(n**3*J*(J+1)*(2*L+1))
    K = F*(F+1)-I*(I+1)-J*(J+1)
    return An*theta_nr*K # en Hz

def E_hfs(n,L,J,F):
    return E_Fermi(n,L,J,F) + hfs_corrections(n,L,J,F)

def theoretical_hyperfine_splitting(n=1):
    x = E_hfs(n,0,1/2,I+1/2) - E_hfs(n,0,1/2,I-1/2)
    print(f'Écart hyperfin théorique pour n={n} :',f'{float(x*1e-6):.6f} MHz')
    
def theoretical_transition_frequency(n1,l1,j1,f1,n2,l2,j2,f2):
    pass


# =============================================================================
# Matrix elements
# =============================================================================

def IX():
    H = np.zeros((N,N),dtype=mpf)
    for i, a in enumerate(LJFmF()):
        for j, b in enumerate(LJFmF()):
            if a.F==b.F and a.mF==b.mF and a.n==b.n:
                _X_ = -sqrt(10*(2*1+1)*(2*b.J+1)*(2*a.J+1)*S*(S+1)*(2*S+1)) \
                    *sqrt((2*b.L+1)*(2*a.L+1))*wigner3j(b.L,2,a.L,0,0,0) \
                    *wigner9j(S,1,S,b.L,2,a.L,b.J,1,a.J)*(-1)**b.L   
                H[i,j] = (-1)**(a.I+b.J+a.F)*wigner6j(a.F,b.J,a.I,1,a.I,a.J) \
                          *sqrt(I*(I+1)*(2*I+1))*_X_
    return H
                         
def IL():                
    H = np.zeros((N,N),dtype=mpf)
    for i, a in enumerate(LJFmF()):
        for j, b in enumerate(LJFmF()):
            if a.F==b.F and a.mF==b.mF and a.L==b.L and a.n==b.n:
                _L_ = (-1)**(S+a.L+b.J+1)*sqrt((2*b.J+1)*(2*a.J+1)) \
                    *wigner6j(b.L,b.J,S,a.J,a.L,1)*sqrt(a.L*(a.L+1)*(2*a.L+1))
                H[i,j] = (-1)**(I+b.J+a.F)*sqrt(I*(I+1)*(2*I+1)) \
                         *wigner6j(a.F,b.J,I,1,I,a.J)*_L_
    return H
             
def IS():                
    H = np.zeros((N,N),dtype=mpf)
    for i, a in enumerate(LJFmF()):
        for j, b in enumerate(LJFmF()):
            if a.F==b.F and a.mF==b.mF and a.L==b.L and a.n==b.n:
                _S_ = (-1)**(S+a.L+b.J+1)*sqrt((2*b.J+1)*(2*a.J+1)) \
                    *wigner6j(S,b.J,a.L,a.J,S,1)*sqrt(S*(S+1)*(2*S+1))
                H[i,j] = (-1)**(I+b.J+a.F)*sqrt(I*(I+1)*(2*I+1)) \
                         *wigner6j(a.F,b.J,I,1,I,a.J)*_S_
    return H
                        
def M_sur_r3(M):
    H = np.zeros((N,N),dtype=mpf)
    for i, a in enumerate(LJFmF()):
        for j, b in enumerate(LJFmF()):
            if a.n==b.n and a.L==b.L and a.L!=0:
                H[i,j] = M[i,j]/(a.n**3*a.L*(a.L+1)*(a.L+1/2))
    return H

def M_delta_r_sur_r2(M):
    H = np.zeros((N,N),dtype=mpf)
    for i, a in enumerate(LJFmF()):
        for j, b in enumerate(LJFmF()):
            if a.n==b.n and a.L==b.L and a.L==0:
                H[i,j] = M[i,j]*4/a.n**3
    return H