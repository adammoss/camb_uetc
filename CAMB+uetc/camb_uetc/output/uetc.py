import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','figure.figsize':[8,4.5],'legend.fontsize':6,'frameon':False,'font.size':6}
plt.rcParams.update(params)

fig = plt.figure()
ax = fig.add_subplot(231)
ax1 = fig.add_subplot(232,sharey=ax)
ax2 = fig.add_subplot(233,sharey=ax1)
ax3 = fig.add_subplot(234,sharex=ax)
ax4 = fig.add_subplot(235,sharex=ax1,sharey=ax3)

filename = 'data/test_VOS4_uetc_'

ktau = np.loadtxt(filename+'ktau.dat')
logktau = np.log(ktau)
ss00 = np.loadtxt(filename+'ss00.dat')
ss = np.loadtxt(filename+'ss.dat')
sscross = np.loadtxt(filename+'sscross.dat')
vv = np.loadtxt(filename+'vv.dat')
tt = np.loadtxt(filename+'tt.dat')

tick_locator = ticker.MaxNLocator(nbins=5)

a = ax.contourf(logktau,logktau,ss00,levels=np.linspace(np.min(ss00),np.max(ss00),50))
ax.set_ylabel('$\log k\eta$',fontsize=8,labelpad=0)
c=plt.colorbar(a,ax=ax,pad=0)
c.locator = tick_locator
c.update_ticks()
plt.setp(ax.get_xticklabels(),visible=False)
a1 = ax1.contourf(logktau,logktau,ss,levels=np.linspace(np.min(ss),np.max(ss),50))
plt.setp(ax1.get_xticklabels(),visible=False)
plt.setp(ax1.get_yticklabels(),visible=False)
c1 = plt.colorbar(a1,ax=ax1,pad=0)
c1.locator = tick_locator
c1.update_ticks()
a2 = ax2.contourf(logktau,logktau,sscross,levels=np.linspace(np.min(sscross),np.max(sscross),50))
plt.setp(ax2.get_yticklabels(),visible=False)
c2 = plt.colorbar(a2,ax=ax2,pad=0)
c2.locator = tick_locator
c2.update_ticks()
a3 = ax3.contourf(logktau,logktau,vv,levels=np.linspace(np.min(vv),np.max(vv),50))
ax3.set_xlabel('$\log k\eta$',fontsize=8,labelpad=0)
ax3.set_ylabel('$\log k\eta$',fontsize=8,labelpad=0)
c3 = plt.colorbar(a3,ax=ax3,pad=0)
c3.locator = tick_locator
c3.update_ticks()
a4 = ax4.contourf(logktau,logktau,tt,levels=np.linspace(np.min(tt),np.max(tt),50))
ax4.set_xlabel('$\log k\eta$',fontsize=8,labelpad=0)
plt.setp(ax4.get_yticklabels(),visible=False)
c4 = plt.colorbar(a4,ax=ax4,pad=0)
c4.locator = tick_locator
c4.update_ticks()

plt.tight_layout()
plt.show()


