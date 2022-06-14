# -*- coding: utf-8 -*-
"""
Created on Mon Oct 2 9:00:00 2017, @author: Simon
"""

from ..config import CODATA_YEAR,ISOTOPE

from mpmath import mp, mpf, pi
DPS = 25 # Decimal precision 
mp.dps = DPS

h = mpf('6.62607015e-34')          # Planck constant, J.s
c = mpf('299792458')               # Speed of light in vacuum, m/s
e = mpf('1.602176634e-19')         # Elementary charge, C
k = mpf('1.380649e-23')            # Boltzmann constant, J/K

CODATA = {2018:{'alpha':      mpf('7.2973525693e-3'),
                'mu_ep':      mpf('5.44617021487e-4'),
                'mu_ed':      mpf('2.724437107462e-4'),
                'mu_emu':     mpf('4.83633169e-3'),
                'Rinfty':     mpf('10973731.568160'),
                'rp':         mpf('8.414e-16'),
                'rd':         mpf('2.12799e-15'),
                'ae':         mpf('1.15965218128e-3'),
                'gp':         mpf('5.5856946893'),
                'gd':         mpf('0.8574382338')},
          2014:{'alpha':      mpf('7.2973525664e-3'),
                'mu_ep':      mpf('5.44617021352e-4'),
                'mu_ed':      mpf('2.724437107484e-4'),
                'mu_emu':     mpf('4.83633170e-3'),
                'Rinfty':     mpf('10973731.568508'),
                'rp':         mpf('8.751e-16'),
                'rd':         mpf('2.1413e-15')}}

alpha = CODATA[CODATA_YEAR]['alpha']
rydberg = CODATA[CODATA_YEAR]['Rinfty']
rp = CODATA[CODATA_YEAR]['rp']
rd = CODATA[CODATA_YEAR]['rd']
mu_ep = CODATA[CODATA_YEAR]['mu_ep']
mu_ed = CODATA[CODATA_YEAR]['mu_ed']
mu_emu = CODATA[CODATA_YEAR]['mu_emu']
ae = CODATA[CODATA_YEAR]['ae']
ap = CODATA[CODATA_YEAR]['gp']/2 - 1
ad = CODATA[CODATA_YEAR]['gd']*mu_ep/mu_ed/2 - 1
hbar = h/(2*pi)

rn = {'PROTIUM':rp,'DEUTERIUM':rd}[ISOTOPE]
an = {'PROTIUM':ap,'DEUTERIUM':ad}[ISOTOPE]
mu_en = {'PROTIUM':mu_ep,'DEUTERIUM':mu_ed}[ISOTOPE]
I = {'PROTIUM':1/2,'DEUTERIUM':1}[ISOTOPE]
S = 1/2