# -*- coding: utf-8 -*-
"""
Created on Mon Oct 2 9:00:00 2017, @author: Simon

Le zéro d'énergie est celle du système {électron + proton} sans interaction.
Les formules donnent un résultat en fréquence, si multipliées par c*rydberg.
(Dans le CODATA, c'est un résultat en énergie, si multipliées par me*c**2)
On rappelle que me*c**2/h = 2/alpha**2 * c*rydberg.
On rappelle que la masse réduite vérifie me/mr = 1 + me/mn.

Précautions avec mpmath : il ne faut pas qu'une expression calculée par python 
(par exemple, l'intérieur d'une parenthèse) commence par une fraction illimitée
qui ne soit pas de type mpf.
Ainsi, 2/3*alpha est incorrect, il faut écrire mpf(2)/3*alpha ou alpha*2/3.
"""

from .codata_values import c,rydberg
from .fine_structure import Dirac,Dirac_2,Dirac_2bis,E_fs
from .fine_structure_corrections import self_energy,vacuum_polarization, \
                  nuclear_size,two_photons,three_photons,nuclear_self_energy, \
                  nuclear_size_corrections,relativistic_recoil,radiative_recoil
from .hyperfine_structure import E_Fermi
from .hyperfine_structure_corrections import                 
from mpmath import mpf

def table_of_fine_structure_corrections(n,l,j):
    n,l,j = mpf(n),mpf(l),mpf(j)
    def recoil(n,l,j):
        return Dirac_2(n,j) + Dirac_2bis(n,l,j) + relativistic_recoil(n,l,j) \
               + radiative_recoil(n,l)
    def total_corrections(n,l,j):
        return E_fs(n,l,j)/c/rydberg - Dirac(n,j)
    K = c*rydberg*1e-3
    def to_str(x): return f'{float(K*x):.3f} kHz'
    print(f'Corrections pour n={n}, L={l}, J={j} :')
    print('Auto-énergie :          ',to_str(self_energy(n,l,j)))
    print('Polarisation du vide :  ',to_str(vacuum_polarization(n,l,j)))
    print('Recul :                 ',to_str(recoil(n,l,j)))
    print('Taille du noyau :       ',to_str(nuclear_size(n,l,j)))
    print('Deux et trois photons : ',to_str(two_photons(n,l,j)+three_photons(n,l,j)))
    print('Auto-énergie nucléaire :',to_str(nuclear_self_energy(n,l)))
    print('Corr. taille du noyau : ',to_str(nuclear_size_corrections(n,l,j)))
    print('Total :                 ',to_str(total_corrections(n,l,j)),'\n')
    
def table_of_hyperfine_structure_corrections(n,L,J,F):
    K = E_Fermi(n,L,J,F)
    def radiatives(n,L,J): 
        return correction_radiative_variable_S(n,L) \
           + correction_radiative_variable_ordre_2_S(n,L) \
           + correction_radiative_constante_S(L) \
           + correction_radiative_hors_S(n,L,J)
    def to_str(x): return f'{float(x*1e-3):.3f} kHz'
    print(f'Corrections pour n={n}, L={L}, J={J}, F={F} :')
    print('Relativistes :         ',to_str(correction_relativiste(n,L,J)*K))
    print('Moment anomal e- :     ',to_str(correction_ae(n,L,J)*K))
    print('Recul :                ',to_str(correction_recul(n,L,J)*K))
    print('Radiatives :           ',to_str(radiatives(n,L,J)*K))
    print('QE :                   ',to_str(correction_QE(n,L,J,F)))
    print('Total :                ',to_str(correction_structure_hyperfine(n,L,J,F)),'\n')