import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib
import scipy.interpolate as s

def xia(tau):
	if tau<t[0]:
		return xi[0]/t[0]
	elif tau>t[len(t)-1]:
		return xi[len(xi)-1]/t[len(t)-1]
	else:
		return xis(tau)/tau

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','figure.figsize':[8,4.5],'legend.fontsize':6,'frameon':False,'font.size':6}
plt.rcParams.update(params)

mask = False

n = 343
j = 0
#filename = 'data/comp/VOS4_512_uetc_'
filename = '/home/ppxtc/Documents/cosmicstrings/camb/data/uetc/test_uetc_'

t = np.loadtxt('data/comp/time4.dat')
xi = np.loadtxt('data/comp/xi4.dat')
xis = s.UnivariateSpline(t,xi,k=5)
k = np.loadtxt('data/comp/k')
xmin = 0.15
xapr = 2.5
etcmin = 0.001

#for i in range(1,2):
for i in range(1,n):
	print i
	fig = plt.figure(j)
	ax = fig.add_subplot(231)
	ax1 = fig.add_subplot(232,sharey=ax)
	ax2 = fig.add_subplot(233,sharey=ax1)
	ax3 = fig.add_subplot(234,sharex=ax)
	ax4 = fig.add_subplot(235,sharex=ax1,sharey=ax3)
	ktau = np.loadtxt(filename+'ktau.dat')
	logktau = np.log(ktau)
	ss00 = np.loadtxt(filename+str(i)+'_ss00.dat')
	ss = np.loadtxt(filename+str(i)+'_ss.dat')
	sscross = np.loadtxt(filename+str(i)+'_sscross.dat')
	vv = np.loadtxt(filename+str(i)+'_vv.dat')
	tt = np.loadtxt(filename+str(i)+'_tt.dat')

	tick_locator = ticker.MaxNLocator(nbins=5)

	a = ax.contourf(logktau,logktau,ss00,levels=np.linspace(np.min(ss00),np.max(ss00),50))#levels=np.linspace(-2,13.5,50))#
	ax.set_ylabel('$\log k\eta$',fontsize=8,labelpad=0)
	c=plt.colorbar(a,ax=ax,pad=0)
	c.locator = tick_locator
	c.update_ticks()
	plt.setp(ax.get_xticklabels(),visible=False)
	a1 = ax1.contourf(logktau,logktau,ss,levels=np.linspace(np.min(ss),np.max(ss),50))#levels=np.linspace(-0.549,2.001,50))#
	plt.setp(ax1.get_xticklabels(),visible=False)
	plt.setp(ax1.get_yticklabels(),visible=False)
	c1 = plt.colorbar(a1,ax=ax1,pad=0)
	c1.locator = tick_locator
	c1.update_ticks()
	a2 = ax2.contourf(logktau,logktau,sscross,levels=np.linspace(np.min(sscross),np.max(sscross),50))#levels=np.linspace(-1.5452,1.481,50))#
	plt.setp(ax2.get_yticklabels(),visible=False)
	c2 = plt.colorbar(a2,ax=ax2,pad=0)
	c2.locator = tick_locator
	c2.update_ticks()
	plt.title('k = '+str(k[i-1]))
	a3 = ax3.contourf(logktau,logktau,vv,levels=np.linspace(np.min(vv),np.max(vv),50))#levels=np.linspace(-0.152,0.6961,50))#
	ax3.set_xlabel('$\log k\eta$',fontsize=8,labelpad=0)
	ax3.set_ylabel('$\log k\eta$',fontsize=8,labelpad=0)
	c3 = plt.colorbar(a3,ax=ax3,pad=0)
	c3.locator = tick_locator
	c3.update_ticks()
	a4 = ax4.contourf(logktau,logktau,tt,levels=np.linspace(np.min(tt),np.max(tt),50))#levels=np.linspace(-0.0491,0.6668,50))#
	ax4.set_xlabel('$\log k\eta$',fontsize=8,labelpad=0)
	plt.setp(ax4.get_yticklabels(),visible=False)
	c4 = plt.colorbar(a4,ax=ax4,pad=0)
	c4.locator = tick_locator
	c4.update_ticks()
	plt.tight_layout()
	j += 1
	if mask:
		xib = []
		for x in ktau/k[i-1]:
			xib.append(xia(x))

		ind = []
		y = 0
		for x in ktau:
			if x < xmin/xib[y]:
				ind.append(y)
			y += 1

		inde = []
		inda = []
		m = 0
		for x in ktau:
			l = 0
			for y in ktau:
				if abs(x*xib[m]-y*xib[l]) < etcmin:
					inde.append([m,l])
				elif abs(x*xib[m]-y*xib[l]) > xapr:
					inda.append([m,l])
				l += 1
			m += 1

		Z=np.zeros([len(ktau),len(ktau)])

		for x in ind:
			for y in ind:
				Z[x,y] +=1

		for x in inde:
			Z[x[0],x[1]] += 1

		for x in inda:
			Z[x[0],x[1]] += 1

		ax.contourf(np.log(ktau),np.log(ktau),Z,alpha=0.1)
		ax1.contourf(np.log(ktau),np.log(ktau),Z,alpha=0.1)
		ax2.contourf(np.log(ktau),np.log(ktau),Z,alpha=0.1)
		ax3.contourf(np.log(ktau),np.log(ktau),Z,alpha=0.1)
		ax4.contourf(np.log(ktau),np.log(ktau),Z,alpha=0.1)
	plt.savefig('/home/ppxtc/Documents/cosmicstrings/camb/plots/uetc/test/no_mask/uetc_'+str(i)+'.jpg')
	plt.close()

plt.show()


