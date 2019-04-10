"""
Script to generate the initial configuration of polymers tethered to a sphere
for use in BDpack

Output files: q.st.dat, CoM.st.dat, rfin.st.dat

Written by: Tiras Y. Lin
Date: 7/13/18
"""

import numpy as np
import math
import sys

#reading in user inputs
mode = input("Random(1) or ideal(2): ")

if mode == 1:
    nchain_pp = input("Number of chains per particle (nchain_pp): ")
elif mode == 2:
    shape = input("Platonic solid (1-5): ")
    if shape == 1:
        nchain_pp = 4
    elif shape == 2:
        nchain_pp = 6
    elif shape == 3:
        nchain_pp = 8
    elif shape == 4:
        nchain_pp = 12
    elif shape == 5:
        nchain_pp = 20
    else:
        print('There are only 5 platonic solids! :(')
        sys.exit()
else:
    print 'not an option :('
    sys.exit()

nchain = input("Number of particles (nchain): ")
#nchain_pp = input("Number of chains per particle (nchain_pp): ")
nseg_ind = input("Number of segments per individual chain (nseg_ind): ")
a = input("Bead radius (a): ")
a_sph = input("Core sphere radius (a_sph): ")

#calculate chain properties
nseg = nseg_ind * nchain_pp
nbead_ind = nseg_ind+1
nbead = nbead_ind * nchain_pp

#printing out the user inputs
print '\n---------------------'
print 'Initial condition:'
print 'mode = ', mode
print 'nchain =    ',nchain
print 'nchain_pp = ', nchain_pp
print 'nseg =      ', nseg
print 'nbead =     ', nbead
print 'nseg_ind =  ', nseg_ind
print 'nbead_ind = ', nbead_ind
print 'a =         ', a
print 'a_sph =     ', a_sph
print '---------------------'

#opening the output files
qst_file = open("q.st.dat","w")
CoMst_file = open("CoM.st.dat","w")
rfinst_file = open("rfin.st.dat","w")

#initialize numpy arrays
qstart = np.zeros((3,nseg_ind,nchain,nchain_pp))
rcmstart = np.zeros((3,nchain,nchain_pp))
rf_in = np.zeros((3,nchain_pp))
rf_in_unit = np.zeros((3,nchain_pp))
r_curr = np.zeros((3))

#first generate random points on the sphere of radius a+a_sph
if mode == 1:
    for ichain_pp in xrange(0,nchain_pp):
        x = np.zeros((3))
        #to ensure that no divisions by a number that is v small will occur
        #while ((np.linalg.norm(x) < .0001) or (x[2]>=0)):
        while np.linalg.norm(x) < .0001:
            x = np.random.normal(0,1,3)
        rf_in_unit[:,ichain_pp] = (x/np.linalg.norm(x))
        rf_in[:,ichain_pp] = rf_in_unit[:,ichain_pp] * (a+a_sph)
        np.savetxt(rfinst_file,rf_in[:,ichain_pp],delimiter=' ',newline=' ', fmt='%1.10f')
        rfinst_file.write('\n')
elif mode == 2:
    gr = (1+math.sqrt(5))/2
    if shape == 1:
        rf_in_unit[:,0] = np.array([1,1,1])
        rf_in_unit[:,1] = np.array([1, -1, -1])
        rf_in_unit[:,2] = np.array([-1, 1, -1])
        rf_in_unit[:,3] = np.array([-1, -1, 1])
    elif shape == 2:
        rf_in_unit[:,0] = np.array([1, 0, 0])
        rf_in_unit[:,1] = np.array([-1, 0, 0])
        rf_in_unit[:,2] = np.array([0, 1, 0])
        rf_in_unit[:,3] = np.array([0, -1, 0])
        rf_in_unit[:,4] = np.array([0, 0, 1])
        rf_in_unit[:,5] = np.array([0, 0, -1])
    elif shape == 3:
        rf_in_unit[:,0] = np.array([1, 1, 1])
        rf_in_unit[:,1] = np.array([1, 1, -1])
        rf_in_unit[:,2] = np.array([1, -1, 1])
        rf_in_unit[:,3] = np.array([-1, 1, 1])
        rf_in_unit[:,4] = np.array([-1, -1, 1])
        rf_in_unit[:,5] = np.array([-1, 1, -1])
        rf_in_unit[:,6] = np.array([1, -1, -1])
        rf_in_unit[:,7] = np.array([-1, -1, -1])
    elif shape == 4:
        rf_in_unit[:,0] = np.array([0, 1, gr])
        rf_in_unit[:,1] = np.array([0, -1, gr])
        rf_in_unit[:,2] = np.array([0, 1, -gr])
        rf_in_unit[:,3] = np.array([0, -1, -gr])
        rf_in_unit[:,4] = np.array([1, gr, 0])
        rf_in_unit[:,5] = np.array([-1, gr, 0])
        rf_in_unit[:,6] = np.array([1, -gr, 0])
        rf_in_unit[:,7] = np.array([-1, -gr, 0])
        rf_in_unit[:,8] = np.array([gr, 0, 1])
        rf_in_unit[:,9] = np.array([gr, 0, -1])
        rf_in_unit[:,10] = np.array([-gr, 0, 1])
        rf_in_unit[:,11] = np.array([-gr, 0, -1])
    elif shape == 5:
        rf_in_unit[:,0] = np.array([1, 1, 1])
        rf_in_unit[:,1] = np.array([1, 1, -1])
        rf_in_unit[:,2] = np.array([1, -1, 1])
        rf_in_unit[:,3] = np.array([-1, 1, 1])
        rf_in_unit[:,4] = np.array([-1, -1, 1])
        rf_in_unit[:,5] = np.array([-1, 1, -1])
        rf_in_unit[:,6] = np.array([1, -1, -1])
        rf_in_unit[:,7] = np.array([-1, -1, -1])
        rf_in_unit[:,8] = np.array([0, 1/gr, gr])
        rf_in_unit[:,9] = np.array([0, -1/gr, gr])
        rf_in_unit[:,10] = np.array([0, 1/gr, -gr])
        rf_in_unit[:,11] = np.array([0, -1/gr, -gr])
        rf_in_unit[:,12] = np.array([1/gr, gr, 0])
        rf_in_unit[:,13] = np.array([-1/gr, gr, 0])
        rf_in_unit[:,14] = np.array([1/gr, -gr, 0])
        rf_in_unit[:,15] = np.array([-1/gr, -gr, 0])
        rf_in_unit[:,16] = np.array([gr, 0, 1/gr])
        rf_in_unit[:,17] = np.array([-gr, 0, 1/gr])
        rf_in_unit[:,18] = np.array([gr, 0, -1/gr])
        rf_in_unit[:,19] = np.array([-gr, 0, -1/gr])

    for ichain_pp in xrange(0,nchain_pp):
        rf_in_unit[:,ichain_pp] = (rf_in_unit[:,ichain_pp]/np.linalg.norm(rf_in_unit[:,ichain_pp]))
        rf_in[:,ichain_pp] = rf_in_unit[:,ichain_pp] * (a+a_sph)
        np.savetxt(rfinst_file,rf_in[:,ichain_pp],delimiter=' ',newline=' ', fmt='%1.10f')
        rfinst_file.write('\n')



#print rf_in_unit
#print rf_in
#sys.exit()


#generate spring vector and center of mass
for ichain in xrange(0,nchain):
    for ichain_pp in xrange(0,nchain_pp):

        r_curr = rf_in[:,ichain_pp].copy()
        rcmstart[:,ichain,ichain_pp] = r_curr[:].copy()

        for iseg in xrange(0,nseg_ind):
            #change to different magnitude of initial spring length if desired
            qstart[:,iseg,ichain,ichain_pp] = 1.0*rf_in_unit[:,ichain_pp].copy()

            r_curr[:] = r_curr[:] + qstart[:,iseg,ichain,ichain_pp]
            rcmstart[:,ichain,ichain_pp] = rcmstart[:,ichain,ichain_pp]+ r_curr[:]

            np.savetxt(qst_file,qstart[:,iseg,ichain,ichain_pp],delimiter=' ',newline=' ', fmt='%1.10f')
            qst_file.write('\n')

        rcmstart[:,ichain,ichain_pp] = rcmstart[:,ichain,ichain_pp]/nbead_ind
        np.savetxt(CoMst_file,rcmstart[:,ichain,ichain_pp],delimiter=' ',newline=' ', fmt='%1.10f')
        CoMst_file.write('\n')

#closing the output files
qst_file.close()
CoMst_file.close()
rfinst_file.close()
