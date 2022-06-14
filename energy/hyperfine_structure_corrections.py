# -*- coding: utf-8 -*-
"""
Created on Mon Oct 2 9:00:00 2017, @author: Simon
"""

from ..config import ISOTOPE 
from ..constant import h,c,alpha,rydberg,mu_en,ae,an,I,S
from ..basis import LJFmF,wigner3j,wigner6j,wigner9j
from .hyperfine_structure import E_Fermi
from mpmath import mpf, sqrt, pi, log, psi

N = len(LJFmF())
E0 = 2*h*c*rydberg/alpha**2         # Énergie propre de l'électron, J
En = E0/mu_en                       # Énergie propre du noyau, J
Er = E0/(1+mu_en)                   # Énergie propre réduite de l'électron, J

An = (1+an)*alpha**4*Er**3/E0/En/h # en Hz
hyperfine_splitting = {'PROTIUM':mpf('1420405751.768'), 
                       'DEUTERIUM':mpf('327384335.2522')}[ISOTOPE]

def correction_relativiste(n,L,J):
    k = (-1)**(J-L+1/2)*(J+1/2)
    gam = sqrt(k**2-alpha**2)
    N = sqrt((n-abs(k))**2 + 2*(n-abs(k))*gam + k**2)
    theta_rel = abs(k)/(J*(J+1))*abs(N-2*k*(gam+n-abs(k)))/(N**4*gam*(4*gam**2-1))
    theta_nr = 1/(n**3*J*(J+1)*(2*L+1))
    return theta_rel/theta_nr - 1
    # return alpha**2*((12*k**2-1)/(2*k**2*(2*k-1)*(2*k+1)) + 3/(2*n*abs(k)) \
    #                   + (3-8*k)/(2*n**2*(2*k-1)))

def correction_ae(n,L,J):
    if L==0: return ae
    else: return ae/(2*(-1)**(J-L+1/2)*(J+1/2))
    
def correction_recul(n,L,J):
    if L==0: return 0
    k = (-1)**(J-L+1/2)*(J+1/2)
    return (1-1/2/k)*(an+1/2)/(an+1)*mu_en

def correction_QE(n,L,J,F):
    if J<1 or I<1: return 0
    Q = 2.86e-31 # Moment quadrupolaire électrique du deutéron (Reid75)
    K = F*(F+1) - I*(I+1) - J*(J+1)
    rmoins3 = (c/alpha/E0)**(-3) / (n**3*L*(L+1)*(L+1/2))
    return c*alpha*Q*rmoins3*(2*J-1)/(2*J+2) \
           *(3*K*(K+1) - 4*I*(I+1)*J*(J+1))/(8*I*(2*I-1)*J*(2*J-1)) 

def correction_radiative_variable_S(n,L):
    if L!=0: return 0
    N = {1:17.85567203, 2:12.03214158, 3:10.449809, 4:9.722413,
         5:9.304114,    6:9.031832,    7:8.840123,  8:8.697639}
    return alpha**3/pi*(N[n] + 15/8/n - 9/8/n**2 \
                           + 8/3*(log(alpha**-2/2)+439/480) \
                                *(psi(0,n)-log(n)-1/n+1/4/n**2))
        
def correction_radiative_variable_ordre_2_S(n,L):
    if L!=0: return 0
    return alpha**2*E0/En*(3/2/n**2 \
            - 7*(1+an)/4*(-9/14/n + 1/28/n**2 + psi(0,n) - log(n)) \
            + (1/4/(1+an)-1/2)*(1/2/n + 5/12/n**2 + psi(0,n) - log(n)))
        
def correction_radiative_constante_S(L):
    if L!=0: return 0
    theta_1S = 1/(1**3*1/2*(1/2+1)*(2*0+1))
    return hyperfine_splitting / (An*theta_1S*(2*I+1)) - 1 \
           - correction_relativiste(1,0,1/2) \
           - correction_ae(1,0,1/2) \
           - correction_recul(1,0,1/2) \
           - correction_radiative_variable_S(1,0) \
           - correction_radiative_variable_ordre_2_S(1,0)
           
def correction_radiative_hors_S(n,L,J):
    if L==0: return 0
    k = (-1)**(J-L+1/2)*(J+1/2)
    if L==1:
        a21 = {1/2:-2*(1-1/n**2), 3/2:0}
        gJ  = {1/2:3.437610,      3/2:0.12609}
        return alpha**3/pi*(a21[J]*log(1/alpha**2)+gJ[J])
    if L>1:
        low = {(3,3/2):-2.06839e-2, (3,5/2):-3.45522e-2}
        high = (4*k+1)*(6*k+1)*(6*k**2+3*k-1)/(8*k**3*(2*k+1)**2*(2*k-1)*(k+1)) \
               + 3*k*(6*k+1)/(n*8*abs(k)*k**2*(2*k+1)) + (4*k-1)/(n**2*2*k*(1-2*k))
        return alpha**3/pi*(low[(n,J)] + high)
    
def hfs_corrections(n,L,J,F):
    return E_Fermi(n,L,J,F) \
           *(correction_relativiste(n,L,J) \
           + correction_ae(n,L,J) \
           + correction_recul(n,L,J) \
           + correction_radiative_variable_S(n,L) \
           + correction_radiative_variable_ordre_2_S(n,L) \
           + correction_radiative_constante_S(L) \
           + correction_radiative_hors_S(n,L,J)) \
           + correction_QE(n,L,J,F)           