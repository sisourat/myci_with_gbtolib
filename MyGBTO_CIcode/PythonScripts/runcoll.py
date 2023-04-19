#!/usr/bin/env python3.10
"""
documentation to be added

run scatci_integrals_vp

implement the collision code (from myFano.f90 see CIcode)

"""

run_scatci_integrals_vp = '/home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/bin/scatci_integrals_vp'
run_coll_coupling_cistates = '/home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTO_CIcode/run_coll_coupling_cistates.exe < input_runcoll.txt'

# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import sys
import os

import numpy as np

from xml.dom import minidom

###################################################################

def read_coll_input():
# finp is the xml file
    finp = sys.argv[1]
    doc = minidom.parse(finp)

    states = doc.getElementsByTagName("States")
    ci_fname = states[0].getAttribute("filename")
    stamin = int(states[0].getAttribute("stamin"))
    stamax = int(states[0].getAttribute("stamax"))
    print(stamin,stamax)

    bimp = doc.getElementsByTagName("ImpactParam")
    bimp_type = str(bimp[0].getAttribute("type"))
    bimp_min = float(bimp[0].getAttribute("bmin"))
    bimp_max = float(bimp[0].getAttribute("bmax"))
    nb = int(bimp[0].getAttribute("nb"))
    print(bimp_type,bimp_min,bimp_max,nb)
    if(bimp_type=='linear' or bimp_type=='Linear' or bimp_type=='LINEAR'):
        bgrid = [bimp_min + i*(bimp_max-bimp_min)/nb for i in range(nb)]
    else:
        raise Exception("Sorry, the impact param. grid must be linear") 
    print(bgrid)  

    vimp = doc.getElementsByTagName("ImpactVel")
    vx = vimp[0].getAttribute("vx")
    vy = vimp[0].getAttribute("vy")
    vz = vimp[0].getAttribute("vz")
    print(vx,vy,vz)

    zgrid = doc.getElementsByTagName("Zgrid")
    zgrid_type = str(zgrid[0].getAttribute("type"))
    zgrid_min = float(zgrid[0].getAttribute("zmin"))
    zgrid_max = float(zgrid[0].getAttribute("zmax"))
    nbz = int(zgrid[0].getAttribute("nzgrid"))
    print(zgrid_type,zgrid_min,zgrid_max,nbz)

    if(zgrid_type=='linear' or zgrid_type=='Linear' or zgrid_type=='LINEAR'):
        zgrid = [zgrid_min + i*(zgrid_max-zgrid_min)/nbz for i in range(nbz)]
    elif(zgrid_type=='exp' or zgrid_type=='Exp' or zgrid_type=='EXP'):
        zgrid = np.zeros(nbz+1)
        spac=0.0
        for i in range(1,nbz//2+1):
            spac+=1.1**i
        spac=(zgrid_max-zgrid_min)/(2.0*spac)
        j=0
        for i in range(nbz//2+1,nbz):
            j+=1
            zgrid[i] = zgrid[i-1]+1.1**j*spac
            zgrid[nbz-i] = -zgrid[i]
        zgrid[nbz]=zgrid_max
        zgrid[0]=-zgrid_max
    else:
        raise Exception("Sorry, the Z grid must be linear or exp") 
    print(zgrid)
   
    istate = doc.getElementsByTagName("InitState")
    ista = istate[0].getAttribute("state")
    print(ista)

    pot = doc.getElementsByTagName("Potential")
    charge = pot[0].getAttribute("charge")
    print(charge)

    molden = doc.getElementsByTagName("Molden")
    molden_fname = molden[0].getAttribute("filename")
    print(molden_fname)

    inp_param = {'ci_fname': ci_fname, 'stamin': stamin, 'stamax': stamax, 'vx' : vx, 'vy': vy, 'vz': vz}
    inp_param['zgrid'] = zgrid
    inp_param['bgrid'] = bgrid
    inp_param['ista'] = int(ista)
    inp_param['charge'] = float(charge)
    inp_param['molden'] = str(molden_fname)

    return inp_param


def runcoll(inp_param):
#    print(inp_param)

    molden_fname = inp_param['molden']

    os.system('rm -f tmp moints_vp moints fort.10 log_file.0')
    l = 'X 100 ' + str(inp_param['charge']) + ' 0.0  0.0 ' + '100000.0'
    lsed = "sed 's/X/"+l+"/g'  "+molden_fname+' > tmp'
    os.system(lsed)
    os.system(run_scatci_integrals_vp)
    os.system('mv moints_vp moints_vp_asymp') 

    ib=0
    for b in inp_param['bgrid']:
       os.system('rm -f matcoll') 
 #      os.system('echo '+str(b) +'> log_file.tot')
       os.system('echo '+ str(b) + ' ' + str(inp_param['vz']) +  '> matcoll')
#       os.system('echo '+  str(inp_param['stamax']-inp_param['stamin']) + ' ' + str(inp_param['stamax']-inp_param['stamin'])  +' >> matcoll')
       os.system('echo '+  str(inp_param['stamax']-inp_param['stamin']) + ' ' + str(len(inp_param['zgrid']))  +' >> matcoll')

       for z in inp_param['zgrid']:
          os.system('echo '+str(z) +'>> matcoll')
          os.system('rm -f tmp moints_vp moints fort.10')
          l = 'X 100 ' + str(inp_param['charge']) + ' ' + str(b) + ' 0.0 ' + str(z) 
          lsed = "sed 's/X/"+l+"/g'  "+molden_fname+' > tmp'
          os.system(lsed)
          os.system(run_scatci_integrals_vp)

          lecho = 'echo ' + str(inp_param['stamax']-inp_param['stamin']) + ' ' + str(inp_param['stamax']-inp_param['stamin']) + ' '  + inp_param['ci_fname'] + ' ' + inp_param['ci_fname'] + ' matcoll_ > input_runcoll.txt' 
          os.system(lecho) #compute the coupling matrix between CI states
          os.system(run_coll_coupling_cistates) #compute the coupling matrix between CI states
          os.system('cat fort.10 >> matcoll') 
       os.system('cp matcoll '+'matcoll'+str(ib)) 
       ib+=1
##          os.system('cat log_file.0 >> log_file.tot') 
##          os.system('head tmp') 

#       os.system(run_coll_prop)  # propagate the TDSE for b

    return 0

if __name__ == '__main__':
    inp_param = read_coll_input()
    print('Read Coll Input Done')
    if(runcoll(inp_param) == 0):
       print('Run Coll Done')
