import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt, pi, exp, linspace, loadtxt
from scipy.optimize import curve_fit

"""
User Inputs
"""
try:
  file_name=int(raw_input('File name [trXSqRel.dat]?:'))
except:
  file_name='trXSqRel.dat'
try:
  rlx_i=int(raw_input('Initial guess for relaxation time [500]?:'))
except:
  rlx_i=500

print 'Use Inputs:'
print 'file: ',file_name
print 'rlx initial guess: ',rlx_i


"""
Definitions
"""
def exp_func(t, a, rlx, b):
  return a*np.exp(-t/rlx)+b


"""
Extracting data
"""

my_data = loadtxt(file_name)
time = my_data[:,2]
Xsq = my_data[:,3]

"""
first loop counts and second loop populates
"""

count=0
it_fltr=0
it_ul=len(Xsq)

for itr in range(0,2):
  for it in range(0,it_ul,1):
    if itr == 0:
      if Xsq[it] < 0.09:
        count=count+1
    else:
      if it == 0:
        time_fltr=np.zeros(count)
        Xsq_fltr =np.zeros(count)
      if Xsq[it] < 0.09:
        time_fltr[it_fltr]=time[it]
        Xsq_fltr [it_fltr]=Xsq [it]
        it_fltr=it_fltr+1


"""
Fitting data
"""

popt, pcov = curve_fit(exp_func,time_fltr,Xsq_fltr,[0., rlx_i, 0.])


"""
Plot your data
"""

plt.plot(time_fltr,Xsq_fltr,'ro',label="Original Data")
plt.plot(time_fltr,exp_func(time_fltr,*popt), 'r-', label="Fitted Curve")

print 'Function: a*exp(-t/rlx)+b'
print 'a: ',popt[0]
print 'rlx: ',popt[1]
print 'c: ',popt[2]


plt.ylabel('$X^2/L_c^2$')
plt.xlabel('$t$')
plt.legend()
plt.show()
