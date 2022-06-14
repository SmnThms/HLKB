# -*- coding: utf-8 -*-
"""
Created on Mon Oct 2 9:00:00 2017, @author: Simon

Le zéro d'énergie est celle du système {électron + proton} sans interaction.
Les formules donnent un résultat en fréquence, si multipliées par c*rydberg.
(Dans le CODATA, c'est un résultat en énergie, si multipliées par me*c**2)
On rappelle que me*c**2/h = 2/alpha**2 * c*rydberg.
On rappelle que la masse réduite vérifie me/mr = 1 + me/mn.
"""

from ..config import ISOTOPE
from .codata_values import c,alpha,rydberg,mu_en
from .fine_structure_corrections import Lamb_shift

from mpmath import mpf, sqrt

d = {True:mpf(1), False:mpf(0)} # Kroenecker delta

def E_Ry(n): # Rydberg energy levels (Hz)
    return c*rydberg/n**2/(1+mu_en)

def f(n,j): # Dirac's f function
    return 1/sqrt(1 + alpha**2/( n -j -1/2 +sqrt((j+1/2)**2-alpha**2) )**2)

def Dirac(n,j): 
    return (f(n,j)-1)/(1+mu_en) * 2/alpha**2

def Dirac_2(n,j): 
    return -(f(n,j)-1)**2/(1+mu_en)**2/(1+1/mu_en) / alpha**2

def Dirac_2bis(n,l,j):
    if l==0 and ISOTOPE==1: return 0
    kappa = (-1)**(j-l+1/2)*(j+1/2)
    return alpha**2/kappa/(2*l+1)/n**3/(1+mu_en)/(1+1/mu_en)**2

def E_Dirac(n,l,j): # Dirac energy levels (Hz)
    n,l,j = mpf(n),mpf(l),mpf(j)
    return c*rydberg*(Dirac(n,j) + Dirac_2(n,j) + Dirac_2bis(n,l,j))

def E_fs(n,l,j): # Total fine structure (Hz)
    return E_Dirac(n,l,j) + Lamb_shift(n,l,j)