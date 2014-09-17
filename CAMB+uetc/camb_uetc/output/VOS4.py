import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si

def xia(tau):
	if tau<t[0]:
		return xi1[0]/t[0]
	elif tau>t[len(t)-1]:
		return xi1[len(xi1)-1]/t[len(t)-1]
	else:
		return xis(tau)/tau

def va(tau):
	if tau<t[0]:
		return v1[0]
	elif tau>t[len(t)-1]:
		return v1[len(v1)-1]
	else:
		return vs(tau)

ktau = np.loadtxt('data/comp/VOS4_512_uetc_ktau.dat')
t = np.loadtxt('data/comp/time4.dat')
v1 = np.loadtxt('data/comp/v4.dat')
xi1 = np.loadtxt('data/comp/xi4.dat')
k = np.loadtxt('data/comp/k')

vs = si.UnivariateSpline(t,v1,k=5)
xis = si.UnivariateSpline(t,xi1,k=5)

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','figure.figsize':[4.5,2.5],'legend.fontsize':6,'frameon':False,'font.size':6}
plt.rcParams.update(params)

xib = []
vb = []
for x in ktau/k[16]:
	xib.append(xia(x))
	vb.append(va(x))

fig = plt.figure()
ax = fig.add_subplot(121)
ax1 = fig.add_subplot(122)

ax.plot(ktau/k[16],vb)
ax1.plot(ktau/k[16],xib)

plt.setp(ax,xlabel='$\eta$',ylabel='$v$',xscale='log')#,xlim=[ti,14373.2])
plt.setp(ax1,xlabel='$\eta$',ylabel='$\\xi$',xscale='log')#,xlim=[ti,14373.2])

plt.tight_layout()
plt.savefig('plots/VOS4.pdf')
plt.show()
