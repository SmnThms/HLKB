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

from ..config import ISOTOPE
from .codata_values import c,alpha,rydberg,rn,mu_en,mu_emu

from mpmath import mpf, pi, log, psi, euler, nsum, nprod, zeta

d = {True:mpf(1), False:mpf(0)} # Kroenecker delta

# =============================================================================
# Numerical values taken from CODATA 2018
# =============================================================================

# Bethe logarithm (Drake, PRA 41 1243, 1990) (n,L)
Bethe_log = {(1,0):mpf('2.984128555765498'),   (2,0):mpf('2.811769893120563'),
             (3,0):mpf('2.767663612491822'),   (4,0):mpf('2.749811840454057'),
             (6,0):mpf('2.735664206935105'),   (8,0):mpf('2.730267260690589'),
             (2,1):mpf('-0.030016708630213'),  (3,1):mpf('-0.038190229385312'),
             (3,2):mpf('-0.005232148140883'),  (4,1):mpf('-0.041954894598086'),
             (4,2):mpf('-0.006740938876975'),  (6,2):mpf('-0.008147203962354'),
             (8,2):mpf('-0.008785042984125'), (12,2):mpf('-0.009342953986099')}

# Relativistic recoil (n,l,j)
G_rec = {(1,0,1/2):mpf('9.720')/pi,       (2,0,1/2):mpf('14.899')/pi, 
         (2,1,1/2):mpf('1.5097')/pi,      (2,1,3/2):mpf('-2.1333')/pi,
         (3,0,1/2):mpf('15.242')/pi,      (4,0,1/2):mpf('15.115')/pi,
         (5,0,1/2):mpf('14.941')/pi,      (6,0,1/2):mpf('14.8')/pi,
         (8,0,1/2):mpf('14.7')/pi}

# Self-energy (n,l,j)
GSE = {(1,0,1/2): mpf('-30.290240'),      (2,0,1/2): mpf('-31.185150'),
       (3,0,1/2): mpf('-31.04770'),       (4,0,1/2): mpf('-30.9120'),
       (6,0,1/2): mpf('-30.711'),         (8,0,1/2): mpf('-30.606'),
       (2,1,1/2): mpf('-0.97350'),        (4,1,1/2): mpf('-1.1640'),
       (2,1,3/2): mpf('-0.48650'),        (4,1,3/2): mpf('-0.6090'),
       (8,2,3/2): mpf('0.007940'),        (12,2,3/2):mpf('0.009130'),
       (4,2,5/2): mpf('0.03163'),         (6,2,5/2): mpf('0.03417'),
       (8,2,5/2): mpf('0.03484'),         (12,2,5/2):mpf('0.03512')}

# Vacuum polarization
GVP = {(1,0,1/2): mpf('-0.618724'),       (2,0,1/2): mpf('-0.808872'),
       (3,0,1/2): mpf('-0.814530'),       (4,0,1/2): mpf('-0.806579'),
       (6,0,1/2): mpf('-0.791450'),       (8,0,1/2): mpf('-0.781197'),
       (2,1,1/2): mpf('-0.064006'),       (4,1,1/2): mpf('-0.080007'),
       (2,1,3/2): mpf('-0.014132'),       (4,1,3/2): mpf('-0.017666')}
ratio_VP_hadrons = 0.671

# Two-photon corrections
B60 = {(1,0,1/2): mpf('-78.7'),           (2,0,1/2): mpf('-63.6'),
       (3,0,1/2): mpf('-60.5'),           (4,0,1/2): mpf('-58.9'),
       (6,0,1/2): mpf('-56.9'),           (8,0,1/2): mpf('-55.9'),
       (2,1,1/2): mpf('-1.8'),            (4,1,1/2): mpf('-2.5'),
       (2,1,3/2): mpf('-1.8'),            (4,1,3/2): mpf('-2.5'),
       (8,2,3/2): mpf('0.245'),           (12,2,3/2):mpf('0.259'),
       (4,2,5/2): mpf('-0.178'),          (6,2,5/2): mpf('-0.207'),
       (8,2,5/2): mpf('-0.221'),          (12,2,5/2):mpf('-0.235')} # (n,l,j)
N   = {(1,0): mpf('17.85567203'),         (2,0): mpf('12.03214158'),
       (3,0): mpf('10.449809'),           (4,0): mpf('9.722413'),
       (6,0): mpf('9.031832'),            (8,0): mpf('8.697639'),
       (2,1): mpf('0.003300635'),         (4,1): mpf('-0.000394332')} # (n,l)

# Nuclear polarizability (Hz) (isotope)
pol_nucl = {'PROTIUM': mpf('0.393e3'),
            'DEUTERIUM': mpf('-21.78e3')+mpf('-0.541e3')}
# Effective Friar radius (m) (isotope)
r_Friar  = {'PROTIUM': mpf('1.947e-15'),
            'DEUTERIUM': mpf('1.43e-15')}


# =============================================================================
# Lamb Shift
# =============================================================================
    
def relativistic_recoil(n,l,j):
    def a_n(n,l):
        if l!=0: return 1/(l*(l+1)*(2*l+1))
        return -2*log(2/n) - 2 + 1/n - 2*nsum(lambda k:1/k,[1,n])
    ES = 2/(1+mu_en)**2*alpha**3/(1+1/mu_en)/pi/n**3 \
       *(d[l==0]*log(1/alpha**2)/3 - Bethe_log.get((n,l),0)*8/3 \
         - d[l==0]/9 - a_n(n,l)*7/3 - 2*d[l==0]*(log(1+mu_en)/(1-mu_en**2) \
                                              + log(1+1/mu_en)/(1-1/mu_en**2)))
    if l==0: D60 = 4*log(2)-7/2
    else:    D60 = 2*(3-(l*(l+1))/n**2)/((2*l-1)*(2*l+1)*(2*l+3))
    return ES + 2*alpha**4*mu_en/n**3*(D60 + alpha*G_rec.get((n,l,j),0))

def radiative_recoil(n,l):
    E_RR = 6*zeta(3) - 2*pi**2*log(2) + pi**2*35/36 - 448/27 \
         + 2/3*pi*alpha*log(1/alpha**2)**2
    return d[l==0]*E_RR/(1+1/mu_en)/(1+mu_en)**2*2*alpha**4/pi**2/n**3

def self_energy(n,l,j):
    A40 = -Bethe_log.get((n,l),0)*4/3 + 10/9*d[l==0] \
        - (1-d[l==0])/2/(2*l+1)/(-1)**(j-l+1/2)/(j+1/2)*(1+mu_en)
    A50 = (mpf(139)/32-2*log(2))*pi*d[l==0]
    A61 = d[l==0]*(4*nsum(lambda k:1/k,[1,n])+28/3*log(2)-4*log(n) \
                 -601/180-77/45/n**2) \
        + d[l==1]*(1-1/n**2)*(mpf(2)/15+1/3*d[j==1/2])
    if l!=0: A61 += (96*n**2-32*l*(l+1))/3/n**2/nprod(lambda k:(2*l+k),[-1,3])
    F = d[l==0]*4/3*log((1+mu_en)/alpha**2) + A40 + alpha*A50 \
      - alpha**2*log((1+mu_en)/alpha**2)**2*d[l==0] \
      + alpha**2*log((1+mu_en)/alpha**2)*A61 + GSE.get((n,l,j),0)*alpha**2
    return F/(1+mu_en)**3*2*alpha**3/pi/n**3

def vacuum_polarization(n,l,j):
    H = d[l==0]*(-mpf(4)/15+5/48*pi*alpha-2/15*alpha**2*log((1+mu_en)/alpha**2)) \
      + GVP.get((n,l,j),0)*alpha**2 \
      + d[l==0]*(mpf(19)/45-pi**2/27+(mpf(1)/16-31*pi**2/2880)*pi*alpha)*alpha**2
    H_muons_hadrons = -d[l==0]*4/15*mu_emu**2 * (1 + ratio_VP_hadrons)
    return (H + H_muons_hadrons)/(1+mu_en)**3*2*alpha**3/pi/n**3

def two_photons(n,l,j):
    B40 = d[l==0]*(pi**2*log(2)*3/2 - 10/27*pi**2 - 2179/648 - 9/4*zeta(3)) \
        + (1-d[l==0])*(pi**2*log(2)/2 - pi**2/12 - 197/144 - 3/4*zeta(3)) \
           /(-1)**(j-l+1/2)/(j+1/2)/(2*l+1)*(1+mu_en)       
    B50 = d[l==0]*mpf('-21.554467196279')
    B63 = d[l==0]*(-mpf(8)/27)
    B62 = d[l==0]*16/9*(mpf(71)/60-log(2)+euler+psi(0,n)-log(n)-1/n+1/4/n**2) \
        + d[l==1]*4/27*(1-1/n**2)
    B61 = d[l==0]*(mpf(413581)/64800+4/3*N.get((n,l),0)+2027/864*pi**2-616/135*log(2) \
             -2/3*log(2)*pi**2+40/9*log(2)**2+zeta(3)+(304/135-32/9*log(2)) \
             *(mpf(3)/4+euler+psi(0,n)-log(n)-1/n+1/4/n**2)-43/36+133/864*pi**2) \
        + d[l==1]*(N.get((n,l),0)*4/3+(1-1/n**2)*(31/405+d[j==1/2]/3-8/27*log(2)))
    B72 = d[l==0]*pi*(-mpf(139)/48+4/3*log(2)-5/72)
    B71 = d[l==0]*(-116 + pi*(mpf(427)/36-16/3*log(2)) \
                            *(3/4-1/n+1/4/n**2+euler+psi(0,n)-log(n))) \
        + d[l==1]*pi*(mpf(139)/144-4/9*log(2)+5/216)*(1-1/n**2)
    F = B40 + B50*alpha + B60.get((n,l,j),0)*alpha**2 \
        + B61*alpha**2*log((1+mu_en)/alpha**2)  \
        + B62*alpha**2*log((1+mu_en)/alpha**2)**2 \
        + B63*alpha**2*log((1+mu_en)/alpha**2)**3 \
        + B71*alpha**3*log((1+mu_en)/alpha**2) \
        + B72*alpha**3*log((1+mu_en)/alpha**2)**2
    return F/(1+mu_en)**3*2*alpha**4/pi**2/n**3

def three_photons(n,l,j):
    # a4 = nsum(lambda n:1/(2**n*n**4),[1,inf])
    a4 = mpf('0.51747906167389938')
    if l==0: C40 = -568*a4/9 + 85*zeta(5)/24 - 121*pi**2*zeta(3)/72 \
          - 84071*zeta(3)/2304 - 71*log(2)**4/27 - 239*pi**2*log(2)**2/135 \
          + 4787*pi**2*log(2)/108 + 1591*pi**4/3240 - 252251*pi**2/9720 \
          + 679441/93312
    else: C40 = (-100*a4/3 + 215*zeta(5)/24 - 83*pi**2*zeta(3)/72 \
          - 139*zeta(3)/18 - 25*log(2)**4/18 + 25*pi**2*log(2)**2/18 \
          + 298*pi**2*log(2)/9 + 239*pi**4/2160 -17101*pi**2/810 -28259/5184) \
          /(-1)**(j-l+1/2)/(j+1/2)/(2*l+1)*(1+mu_en)
    X = (-mpf(1523)/648 - 10/27*pi**2 + 3/2*pi**2*log(2) - 9/4*zeta(3) - 82/81)
    C61 = d[l==1]*2/9*(1-1/n**2)*X*(1+mu_en)**3
    C62 = d[l==0]*-2/3*X
    F = C40 + C61*alpha**2*log((1+mu_en)/alpha**2) \
        + C62*alpha**2*log((1+mu_en)/alpha**2)**2
    return F/(1+mu_en)**3*2*alpha**5/pi**3/n**3

def nuclear_size(n,l,j):  
    if l!=0: return 0
    lC = alpha**2/(4*pi*rydberg)
    return mpf(2)/3*(rn/lC)**2*2*alpha**2/(1+mu_en)**3/n**3

def nuclear_size_corrections(n,l,j):
    lC = alpha**2/(4*pi*rydberg)
    kappa = (-1)**(j-l+1/2)*(j+1/2)
    E5 = d[l==0]*(-1/3)*alpha**3*(r_Friar[1]/lC)**3
    if ISOTOPE=='DEUTERIUM': 
        E5 += d[l==0]*((-1/3)*alpha**3*(r_Friar[2]/lC)**3)
    E6 = alpha**4*(rn/lC)**2 \
       * ( d[l==0]*-2/3*(mpf(9)/4/n**2 - 3 - 1/n + 2*euler - log(n/2) + psi(0,n) \
                           + log(1.068497*alpha*rn/lC/(1+mu_en)) \
                           - (4*log(2) - 5)) \
         + d[kappa==1]*(1-1/n**2)/6 )
    if n!=1 and l==1: delta_nP = 4*N.get((n,l),0)/(1-1/n**2)
    else: delta_nP = 0
    E7 = alpha**5*(rn/lC)**2/pi \
       * ( d[l==0]*2/3*(-log(alpha**-2)**2*2/3 + log(rn/lC/(1+mu_en)**2)) \
         + d[l==1]*(1-1/n**2)/6*(log(2*alpha**-2)*8/9 + 11/27 + d[kappa==1] \
                                 + delta_nP))
    elastic = (E5+E6+E7)*2/(1+mu_en)**3/n**3
    inelastic = d[l==0]*pol_nucl.get(ISOTOPE)/n**3 / (c*rydberg)
    return elastic + inelastic

def nuclear_self_energy(n,l):
    E_SEN = d[l==0]*log((1+1/mu_en)/alpha**2) - Bethe_log.get((n,l),0)
    return E_SEN*4/3/(1+1/mu_en)**2/(1+mu_en)*2*alpha**3/pi/n**3

def Lamb_shift(n,l,j):
    n,l,j = mpf(n),mpf(l),mpf(j)
    return c*rydberg*(relativistic_recoil(n,l,j) + self_energy(n,l,j) \
           + vacuum_polarization(n,l,j) + two_photons(n,l,j) \
           + three_photons(n,l,j) + nuclear_size(n,l,j) \
           + nuclear_size_corrections(n,l,j) \
           + radiative_recoil(n,l) + nuclear_self_energy(n,l))