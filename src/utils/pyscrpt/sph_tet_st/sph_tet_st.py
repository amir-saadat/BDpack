"""
Script to generate the initial configuration of polymers tethered to a sphere
for use in BDpack

Output files: q.st.dat, CoM.st.dat, rfin.st.dat

Written by: Tiras Y. Lin
Date: 2/19/18
"""

import numpy as np

#reading in user inputs
nchain = input("Number of particles (nchain): ")
nchain_pp = input("Number of chains per particle (nchain_pp): ")
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

#first generate random points on the sphere of radius a+a_sph
for ichain_pp in xrange(0,nchain_pp):
    x = np.zeros((3))
    #to ensure that no divisions by a number that is v small will occur
    while np.linalg.norm(x) < .0001:
        x = np.random.normal(0,1,3)
    rf_in_unit[:,ichain_pp] = (x/np.linalg.norm(x))
    rf_in[:,ichain_pp] = rf_in_unit[:,ichain_pp] * (a+a_sph)
    np.savetxt(rfinst_file,rf_in[:,ichain_pp],delimiter=' ',newline=' ', fmt='%1.8f')
    rfinst_file.write('\n')

#generate spring vector and center of mass
for ichain in xrange(0,nchain):
    for ichain_pp in xrange(0,nchain_pp):
        #initialize center of mass to tether point
        rcmstart[:,ichain,ichain_pp] = rf_in[:,ichain_pp]

        for iseg in xrange(0,nseg_ind):
            qstart[:,iseg,ichain,ichain_pp] = rf_in_unit[:,ichain_pp]
            rcmstart[:,ichain,ichain_pp] = rcmstart[:,ichain,ichain_pp] + qstart[:,iseg,ichain,ichain_pp]

            np.savetxt(qst_file,qstart[:,iseg,ichain,ichain_pp],delimiter=' ',newline=' ', fmt='%1.8f')
            qst_file.write('\n')

        rcmstart[:,ichain,ichain_pp] = rcmstart[:,ichain,ichain_pp]/nbead_ind
        np.savetxt(CoMst_file,rcmstart[:,ichain,ichain_pp],delimiter=' ',newline=' ', fmt='%1.8f')
        CoMst_file.write('\n')

#closing the output files
qst_file.close()
CoMst_file.close()
rfinst_file.close()
