import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si


t = np.loadtxt('data/comp/tauval.dat')
v = np.loadtxt('data/comp/vval.dat')
xi = np.loadtxt('data/comp/xival.dat')

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','figure.figsize':[4.5,2.5],'legend.fontsize':6,'frameon':False,'font.size':6}
plt.rcParams.update(params)

fig = plt.figure()
ax = fig.add_subplot(121)
ax1 = fig.add_subplot(122)

ax.plot(t,v)
ax1.plot(t,xi)

plt.setp(ax,xlabel='$\eta$',ylabel='$v$',xscale='log')#,xlim=[ti,14373.2])
plt.setp(ax1,xlabel='$\eta$',ylabel='$\\xi$',xscale='log')#,xlim=[ti,14373.2])

plt.tight_layout()
plt.savefig('plots/VOS4.pdf')
plt.show()
