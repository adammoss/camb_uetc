import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si

a1 = np.loadtxt('data/a4.dat')
adot1 = np.loadtxt('data/adot4.dat')
time1 = np.loadtxt('data/time4.dat')
v1 = np.loadtxt('data/v4.dat')
xi1 = np.loadtxt('data/xi4.dat')/time1

n = 100
ti = np.min(time1)
tf = np.max(time1)
dt = np.log(tf/ti)/(n-1)
time = np.zeros(n)
for i in range(1,n):
	time[i] = ti*np.exp((i-1)*dt)

a = si.UnivariateSpline(time1,a1,k=3)
adot = si.UnivariateSpline(time1,adot1,k=3)
v = si.UnivariateSpline(time1,v1,k=3)
xi = si.UnivariateSpline(time1,xi1,k=3)

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','figure.figsize':[4.5,2.5],'legend.fontsize':6,'frameon':False,'font.size':6}
plt.rcParams.update(params)

fig = plt.figure()
ax = fig.add_subplot(121)
ax1 = fig.add_subplot(122)

ax.plot(time,v(time))
ax1.plot(time,xi(time))

plt.setp(ax,xlabel='$\eta$',ylabel='$v$',xscale='log')#,xlim=[ti,14373.2])
plt.setp(ax1,xlabel='$\eta$',ylabel='$\\xi$',xscale='log')#,xlim=[ti,14373.2])

plt.tight_layout()
plt.savefig('plots/VOS4.pdf')
plt.show()
