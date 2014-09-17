import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as s

def xia(tau):
	if tau<t[0]:
		return xi[0]/t[0]
	elif tau>t[len(t)-1]:
		return xi[len(xi)-1]/t[len(t)-1]
	else:
		return xis(tau)/tau


t = np.loadtxt('data/comp/time4.dat')
xi = np.loadtxt('data/comp/xi4.dat')
xis = s.UnivariateSpline(t,xi,k=5)
ktau = np.loadtxt('data/comp/uetc_uetc_ktau.dat')
k = 7.396588024304393E-006
xmin = 0.11
xapr = 2.5
etcmin = 0.001

xib = []
for i in ktau/k:
	xib.append(xia(i))

ind = []
j = 0
for i in ktau:
	if i < xmin/xib[j]:
		ind.append(j)
	j += 1

inde = []
inda = []
k = 0
for i in ktau:
	l = 0
	for j in ktau:
		if abs(i*xib[k]-j*xib[l]) < etcmin:
			inde.append([k,l])
		elif abs(i*xib[k]-j*xib[l]) > xapr:
			inda.append([k,l])
		elif j*xib[l]-i*xib[k] > xapr:
			indai.append([k,l])

		l += 1
	k += 1



Z=np.zeros([len(ktau),len(ktau)])

for i in ind:
	for j in ind:
		Z[i,j] +=1

for i in inde:
	Z[i[0],i[1]] += 1

for i in inda:
	Z[i[0],i[1]] += 1

plt.contourf(np.log(ktau),np.log(ktau),Z,alpha=0.1)

plt.show()
