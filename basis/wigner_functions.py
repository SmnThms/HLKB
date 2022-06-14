# -*- coding: utf-8 -*-
"""
Created on Mon Oct 2 9:00:00 2017, @author: Simon
"""

from mpmath import sqrt, factorial as fact
import numpy as np

def clebsch(j1,m1,j2,m2,J,M):
    return (-1)**(j1-j2+M)*np.sqrt(2*J+1)*wigner3j(j1,j2,J,m1,m2,-M)
    
def wigner3j(j1,j2,j3,m1,m2,m3):
    if m1+m2+m3!=0: return 0
    if j1-m1!=np.floor(j1-m1) or j2-m2!=np.floor(j2-m2) \
       or j3-m3!=np.floor(j3-m3): return 0
    if j3>j1+j2 or j3<abs(j1-j2): return 0
    t1, t2, t3, t4, t5 = j2-m1-j3, j1+m2-j3, j1+j2-j3, j1-m1, j2+m2
    tmin, tmax = max(0, max(t1,t2)), min(t3, min(t4,t5))
    wigner = 0
    for t in np.arange(tmin,tmax+1,1):
        wigner = wigner + (-1)**t/(fact(t)*fact(t-t1)*fact(t-t2) \
                                  *fact(t3-t)*fact(t4-t)*fact(t5-t))
    return wigner * (-1)**(j1-j2-m3) * sqrt(fact(j1+j2-j3)*fact(j1-j2+j3) \
              *fact(-j1+j2+j3)/fact(j1+j2+j3+1)*fact(j1+m1)*fact(j1-m1) \
              *fact(j2+m2)*fact(j2-m2)*fact(j3+m3)*fact(j3-m3))
           
def wigner6j(j1,j2,j3,J1,J2,J3):
    if abs(j1-j2)>j3 or j1+j2<j3 or abs(j1-J2)>J3 or j1+J2<J3 \
    or abs(J1-j2)>J3 or J1+j2<J3 or abs(J1-J2)>j3 or J1+J2<j3: return 0
    if 2*(j1+j2+j3)!=round(2*(j1+j2+j3)) or 2*(j1+J2+J3)!=round(2*(j1+J2+J3)) \
    or 2*(J1+j2+J3)!=round(2*(J1+j2+J3)) or 2*(J1+J2+j3)!=round(2*(J1+J2+j3)):
        return 0
    t1, t2, t3, t4 = j1+j2+j3, j1+J2+J3, J1+j2+J3, J1+J2+j3
    t5, t6, t7 = j1+j2+J1+J2, j2+j3+J2+J3, j1+j3+J1+J3
    tmin, tmax = max(0,max(t1,max(t2,max(t3,t4)))), min(t5,min(t6,t7))
    wigner = 0
    for t in np.arange(tmin,tmax+1,1):
        wigner = wigner + (-1)**t*fact(t+1)/(fact(t-t1)*fact(t-t2)*fact(t-t3) \
                                 *fact(t-t4)*fact(t5-t)*fact(t6-t)*fact(t7-t))
    return wigner*sqrt(TriCoef(j1,j2,j3)*TriCoef(j1,J2,J3) \
                 *TriCoef(J1,j2,J3)*TriCoef(J1,J2,j3))
        
def wigner9j(j11,j12,j13,j21,j22,j23,j31,j32,j33):
    return np.sum([(-1)**(2*k)*(2*k+1)*wigner6j(j11,j21,j31,j32,j33,k) \
                   *wigner6j(j12,j22,j32,j21,k,j23)\
                   *wigner6j(j13,j23,j33,k,j11,j12) for k in range(10)])

def TriCoef(a,b,c): return fact(a+b-c)*fact(a-b+c)*fact(-a+b+c)/(fact(a+b+c+1))

def WignEck(ja,ma,jb,mb,k,q): return (-1)**(jb-mb)*wigner3j(jb,k,ja,-mb,q,ma)