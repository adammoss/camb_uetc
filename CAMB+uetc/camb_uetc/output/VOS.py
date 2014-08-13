import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si

a1 = np.loadtxt('data/a.dat')
adot1 = np.loadtxt('data/adot.dat')
alpha1 = np.loadtxt('data/alpha.dat')
time1 = np.loadtxt('data/time.dat')
v1 = np.loadtxt('data/v.dat')
xi1 = np.loadtxt('data/xi.dat')/time1

n = 100
ti = np.min(time1)
tf = np.max(time1)
dt = np.log(tf/ti)/(n-1)
time = np.zeros(n)
for i in range(1,n):
	time[i] = ti*np.exp((i-1)*dt)

a = si.UnivariateSpline(time1,a1,k=3)
adot = si.UnivariateSpline(time1,adot1,k=3)
alpha = si.UnivariateSpline(time1,alpha1,k=3)
v = si.UnivariateSpline(time1,v1,k=3)
xi = si.UnivariateSpline(time1,xi1,k=3)

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','figure.figsize':[6,2],'legend.fontsize':6,'frameon':False,'font.size':6}
plt.rcParams.update(params)

fig = plt.figure()
ax = fig.add_subplot(131)
ax1 = fig.add_subplot(132)
ax2 = fig.add_subplot(133)
#ax3 = fig.add_subplot(234)
#ax4 = fig.add_subplot(235)

#ax.plot(time,a(time))
#ax1.plot(time,adot(time))
ax.plot(time,alpha(time))
ax1.plot(time,v(time))
ax2.plot(time,xi(time))

#plt.setp(ax,xlabel='$\eta$',ylabel='$a$',xscale='log',xlim=[ti,14373.2],ylim=[np.min(a(time)),3])
#plt.setp(ax1,xlabel='$\eta$',ylabel='$\dot{a}$',xscale='log',xlim=[ti,14373.2],ylim=[np.min(adot(time)),0.005])
plt.setp(ax,xlabel='$\eta$',ylabel='$\\alpha$',xscale='log')#,xlim=[ti,14373.2])
plt.setp(ax1,xlabel='$\eta$',ylabel='$v$',xscale='log')#,xlim=[ti,14373.2])
plt.setp(ax2,xlabel='$\eta$',ylabel='$\\xi$',xscale='log')#,xlim=[ti,14373.2])

plt.tight_layout()
plt.savefig('plots/VOS.pdf')
plt.show()
