# -*- coding: utf-8 -*-
"""
Created on Mon Oct 2 9:00:00 2017, @author: Simon
"""

from mpmath import mpf

class Level:
    def __init__(self,n=None, L=None, J=None, F=None, I=None, S=None,
                      mS=None, mL=None, mI=None, mJ=None, mF=None):
        for k in ['n','L','J','F','mS','mL','mI','mJ','mF','I','S']: 
            if locals()[k] is None: setattr(self,k,None)
            else: setattr(self,k,mpf(str(locals()[k])))
    def __repr__(self):
        def formattage(x):
            if x%1==0.5: s = str(int(2*x))+'/2  '
            else : s = str(int(x))+'  '
            if x<0: return s
            else: return s+' '
        return '| '+''.join(['='.join([k,formattage(v)]) for k,v \
                               in self.__dict__.items() if v is not None])+'>'