import sys 
import numpy as np
from numpy import linalg as LA

fmo = open("moints_vp","r")
nmo = int(fmo.readline())+1

moints_k = np.zeros((nmo,nmo))
moints_v = np.zeros((nmo,nmo))
moints_vasymp = np.zeros((nmo,nmo))
for l in fmo:
 dat = l.split()
 i = int(dat[0])
 j = int(dat[1])
 moints_k[i,j] = float(dat[3])
 moints_v[i,j] = float(dat[4])

fmo.close()

fmo = open("moints_vp_asymp","r")
nmo_asymp = int(fmo.readline())+1
for l in fmo:
 dat = l.split()
 i = int(dat[0])
 j = int(dat[1])
 moints_vasymp[i,j] = float(dat[4])

fmo.close()

#Hproj = moints_k + moints_vasymp # - moints_vasymp
Hproj = moints_k + moints_v - moints_vasymp

w, v = LA.eig(Hproj)
print(np.sort(w))
