import sys
from spineigenfunctions import *

finp = open(sys.argv[1],'r')

for l in finp:
    w = l.split()
    ip = w.index("p")
    iup = w.index("up")

    spin = float(w[1])

    paired = []
    unpaired = []
    for i in range(ip+1,iup):
        paired.append(int(w[i]))
    for i in range(iup+1,len(w)):
        unpaired.append(int(w[i]))

    iprint = printCSF(spin,paired,unpaired)
    
