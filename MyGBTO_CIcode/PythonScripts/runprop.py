#!/usr/bin/env python3.10
"""
documentation to be added

run Prop

"""

run_prop = '/home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTO_CIcode/lib_prop/Prop'

import sys
import os

import numpy as np

from xml.dom import minidom

if __name__ == '__main__':
    nbb = int(sys.argv[1])
    for ib in range(nbb):
       os.system('cp matcoll'+str(ib)+' matcoll') 
       os.system(run_prop) 
       os.system('cp fort.100 '+' prob'+str(ib)) 
