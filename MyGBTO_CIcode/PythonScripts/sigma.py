#!/usr/bin/env python3.10
"""
documentation to be added

"""

import sys
import os

import numpy as np

from xml.dom import minidom

if __name__ == '__main__':
    esta = np.loadtxt('eigval')
    nsta = len(esta)
    nbb = int(sys.argv[1])
    prob = np.zeros((nbb,nsta+1))
    bprob = np.zeros((nbb,nsta+1))
    b = np.zeros(nbb)
    for ib in range(nbb):
       dat = np.loadtxt('prob'+str(ib))
       b[ib] = dat[0]
       prob[ib] = dat[1:nsta+2]
       bprob[ib] = dat[1:nsta+2]*dat[0]

    sig = []
    for ista in range(nsta):
       sig.append(np.trapz(bprob[:,ista], x=b))
    print(sig[0]*2.0*np.pi*0.28)
    print(np.sum(sig[1:])*2.0*np.pi*0.28)
    print(np.sum(sig)*2.0*np.pi*0.28)
       
