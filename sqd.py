# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 18:09:26 2016

@author: kezhili
"""

import sys

def squared(x):
    y = x * x
    return y

if __name__ == '__main__':
    x = float(sys.argv[1])
    sys.stdout.write(str(squared(x)))