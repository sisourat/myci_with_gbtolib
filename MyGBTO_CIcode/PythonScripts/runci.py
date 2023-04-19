#!/usr/bin/env python3.10

# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import sys

import itertools
import numpy as np

from xml.dom import minidom
from collections import Counter

from spineigenfunctions import *

###################################################################

def ormas_like():
# finp is the xml file
    finp = sys.argv[1]
#    finp = 'input.xml'
    doc = minidom.parse(finp)

# general info
    general = doc.getElementsByTagName("general")
    outfile = general[0].getElementsByTagName('outfile')[0].childNodes[0].data
    elec = general[0].getElementsByTagName('electron')
    spin = int(general[0].getElementsByTagName('spin')[0].childNodes[0].data)
    na = elec[0].getAttribute("na")
    nb = elec[0].getAttribute("nb")

# Get all the blocks
    blocks = doc.getElementsByTagName("block")

    nCSFs = 0

    for block in blocks:
        spaces = block.getElementsByTagName('space')
        seqs = []
        occspace = []
        for space in spaces:
            ne = int(space.getAttribute("ne"))
            ip = int(space.getAttribute("imo"))
            fp = int(space.getAttribute("fmo"))

            seq = np.arange(ip, fp + 1)
            seqtot = np.concatenate((seq, seq), axis=0)
            seqs = np.concatenate((seqs, seq), axis=0)

            s = []
            ind = []
            for l in list(itertools.combinations(seqtot, ne)):
                count = []
                for i in range(ip, fp + 1):
                    count.append(l.count(i))
                s.append(count)

        # remove duplicates
            occ = []
            for i in s:
                if i not in occ:
                    occ.append(i)

            occspace.append(occ)

        if (len(occspace) == 1):
            seqocc = occspace[0]

        for i in range(len(occspace) - 1):
            if (i == 0):
                seqocc = occspace[i]
            s = []
            for j in seqocc:
                for k in occspace[i + 1]:
                    s.append(j + k)
                seqocc = s

        for occ in seqocc:
            p = []
            unp = []
            for i in range(len(occ)):
                if (occ[i] == 2):
                    p.append(int(seqs[i]))
                if (occ[i] == 1):
                    unp.append(int(seqs[i]))

            nCSFs += printCSF(spin, p, unp)

# print info in CI code format
    fhead = open('header.txt', 'a')
    print(nCSFs, '0.000000000000001', file=fhead)
    print(outfile, file=fhead)
    print(na, nb, file=fhead)
    print(nCSFs, file=fhead)
    fhead.close()

    return 0

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    import os

    if os.path.exists("header.txt"):
        os.remove("header.txt")
    if os.path.exists("list.txt"):
        os.remove("list.txt")
    if(ormas_like()==0):
        os.system("cat header.txt list.txt > input_cicode.txt")
        os.remove("header.txt")
        os.remove("list.txt")
        os.system("/home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTO_CIcode/myCI.exe input_cicode.txt")
        print('Job Done')
