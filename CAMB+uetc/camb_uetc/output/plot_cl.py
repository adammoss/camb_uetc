import numpy as np
import matplotlib.pyplot as plt


plotact = True

gmu = 2E-7

filename = ['data/test_init','data/test_VOS','data/test_VOS4']
linestyles = ['-',':','-.','--']
colors = ['red','blue','green','black']

l_range = [2,3000]
tt = [[0,200],[0.0,100],[0.0,6]]
ee = [[0.0001,10],[0.00005,0.1],[0.00001,0.02]]
bb = [[0,1],[0.0001,1],[0.00001,0.01]]
te = [[-3,0.5],[0.0,0.5],[-0.02,0.015]]

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','figure.figsize':[6,4],'legend.fontsize':6,'frameon':False,'font.size':6}
plt.rcParams.update(params)

fig = plt.figure()
h=0

for k in filename:
	data = np.loadtxt(k+'_scalCls.dat')
	data[:,1] = data[:,1]*gmu**2
	data[:,2] = data[:,2]*gmu**2
	data[:,3] = data[:,3]*gmu**2
	j = 1 
	ax = fig.add_subplot(3,4,j)
	ax1 = fig.add_subplot(3,4,j+1)
	ax2 = fig.add_subplot(3,4,j+2)
	ax3 = fig.add_subplot(3,4,j+3)
	ax.plot(data[:,0],data[:,1],color=colors[h],linestyle=linestyles[h])
	ax1.plot(data[:,0],data[:,2],color=colors[h],linestyle=linestyles[h])
	ax3.plot(data[:,0],data[:,3],color=colors[h],linestyle=linestyles[h])
	plt.setp(ax,xscale='log',xlim=l_range,ylim=tt[j-1])
	ax.set_xlabel('$l$',labelpad=0)
	ax.set_ylabel('$C^{TT}_l\,l(l+1)/2\pi\,[\mu K^2]$',labelpad=0)
	plt.setp(ax1,xscale='log',yscale='log',xlim=l_range,ylim=ee[j-1])
	ax1.set_xlabel('$l$',labelpad=0)
	ax1.set_ylabel('$C^{EE}_l\,l(l+1)/2\pi\,[\mu K^2]$',labelpad=0)
	plt.setp(ax2,xscale='log',xlim=l_range,ylim=bb[j-1])
	ax2.set_xlabel('$l$',labelpad=0)
	ax2.set_ylabel('$C^{BB}_l\,l(l+1)/2\pi\,[\mu K^2]$',labelpad=0)
	plt.setp(ax3,xscale='log',xlim=l_range,ylim=te[j-1])
	ax3.set_xlabel('$l$',labelpad=0)
	ax3.set_ylabel('$C^{TE}_l\,l(l+1)/2\pi\,[\mu K^2]$',labelpad=0)

	j=1
	for i in ['_vecCls.dat','_tensCls.dat']:
		data = np.loadtxt(k+i)
		data[:,1] = data[:,1]*gmu**2
		data[:,2] = data[:,2]*gmu**2
		data[:,3] = data[:,3]*gmu**2
		data[:,4] = data[:,4]*gmu**2
		ax = fig.add_subplot(3,4,j*4+1)
		ax1 = fig.add_subplot(3,4,j*4+2)
		ax2 = fig.add_subplot(3,4,j*4+3)
		ax3 = fig.add_subplot(3,4,j*4+4)
		ax.plot(data[:,0],data[:,1],color=colors[h],linestyle=linestyles[h])
		ax1.plot(data[:,0],data[:,2],color=colors[h],linestyle=linestyles[h])
		ax2.plot(data[:,0],data[:,3],color=colors[h],linestyle=linestyles[h])
		ax3.plot(data[:,0],data[:,4],color=colors[h],linestyle=linestyles[h])
		plt.setp(ax,xscale='log',xlim=l_range,ylim=tt[j])
		ax.set_xlabel('$l$',labelpad=0)
		ax.set_ylabel('$C^{TT}_l\,l(l+1)/2\pi\,[\mu K^2]$',labelpad=0)
		plt.setp(ax1,xscale='log',yscale='log',xlim=l_range,ylim=ee[j])
		ax1.set_xlabel('$l$',labelpad=0)
		ax1.set_ylabel('$C^{EE}_l\,l(l+1)/2\pi\,[\mu K^2]$',labelpad=0)
		plt.setp(ax2,xscale='log',yscale='log',xlim=l_range,ylim=bb[j])
		ax2.set_xlabel('$l$',labelpad=0)
		ax2.set_ylabel('$C^{BB}_l\,l(l+1)/2\pi\,[\mu K^2]$',labelpad=0)
		plt.setp(ax3,xscale='log',xlim=l_range,ylim=te[j])
		ax3.set_xlabel('$l$',labelpad=0)
		ax3.set_ylabel('$C^{TE}_l\,l(l+1)/2\pi\,[\mu K^2]$',labelpad=0)
		j += 1 

	h += 1

plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
plt.tight_layout()

if plotact:
	cmbact_fact = 1
	tt = np.loadtxt('../cmbact4/cl_tt_100.d')
	stt = tt[:,1]*(gmu/1.1E-6)**2*cmbact_fact
	vtt = tt[:,2]*(gmu/1.1E-6)**2*cmbact_fact
	ttt = tt[:,3]*(gmu/1.1E-6)**2*cmbact_fact
	ee = np.loadtxt('../cmbact4/cl_ee_100.d')
	see = ee[:,1]*(gmu/1.1E-6)**2*cmbact_fact
	vee = ee[:,2]*(gmu/1.1E-6)**2*cmbact_fact
	tee = ee[:,3]*(gmu/1.1E-6)**2*cmbact_fact
	bb = np.loadtxt('../cmbact4/cl_bb_100.d')
	sbb = np.zeros(len(bb[:,0]))
	vbb = bb[:,1]*(gmu/1.1E-6)**2*cmbact_fact
	tbb = bb[:,2]*(gmu/1.1E-6)**2*cmbact_fact
	te = np.loadtxt('../cmbact4/cl_te_100.d')
	ste = te[:,1]*(gmu/1.1E-6)**2*cmbact_fact
	vte = te[:,2]*(gmu/1.1E-6)**2*cmbact_fact
	tte = te[:,3]*(gmu/1.1E-6)**2*cmbact_fact
	ls = [tt[:,0],ee[:,0],bb[:,0],te[:,0],tt[:,0],ee[:,0],bb[:,0],te[:,0],tt[:,0],ee[:,0],bb[:,0],te[:,0],tt[:,0],ee[:,0],bb[:,0],te[:,0]]
        stuff = [stt,see,sbb,ste,vtt,vee,vbb,vte,ttt,tee,tbb,tte]
	for i in range(len(stuff)):
		ax = fig.add_subplot(3,4,i+1)
		ax.plot(ls[i],stuff[i],color='yellow')


plt.savefig('plots/cls.pdf')

plt.show()
	
