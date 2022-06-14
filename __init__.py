# -*- coding: utf-8 -*-
"""
Created on Mon Oct 2 9:00:00 2017, @author: Simon
"""

import os,sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from config import ISOTOPE,CODATA_YEAR

if __name__ == '__main__':
    print(f'{ISOTOPE = }\t{CODATA_YEAR = }')