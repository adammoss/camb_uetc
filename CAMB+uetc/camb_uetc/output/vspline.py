import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt('data/vsp.dat')
ax = plt.plot(A[:,0],A[:,1],'x')
plt.xscale('log')
plt.show()
