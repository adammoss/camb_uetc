import numpy as np
import matplotlib.pyplot as plt
import operator

def converter(x):
	try:
		z = float(x)
	except:
		z = 0
	return z

	

one = np.loadtxt('data/comp/k1_uetc_tt_evec.dat')
for i in range(1,14):
	filename = 'data/comp/k'+str(i)+'_uetc_'
	filename1 = 'data/comp/k'+str(i+1)+'_uetc_'
	fktau = 'ktau.dat'
	fende = 'tt_eval.dat'
	fendv = 'tt_evec.dat'

	ktau = np.loadtxt(filename+fktau)
	logktau = np.log(ktau)
	tt_eval = np.loadtxt(filename+fende,converters={0:converter})
	tt_evec = np.loadtxt(filename+fendv)

	plt.figure(1)
	plt.plot(logktau,tt_evec[0,:]*tt_eval[0])
	plt.figure(2)
	plt.plot(logktau,tt_evec[1,:]*tt_eval[1])
	plt.figure(3)
	plt.plot(logktau,tt_evec[2,:]*tt_eval[2])
'''	
	l = 0
	for j in range(0,len(one[:,0])):
		if len(set(one[1,:]) & set(tt_evec[j,:])) > 10:
			print len(set(one[0,:]) & set(tt_evec[j,:])),filename,(set(one[0,:]) & set(tt_evec[j,:]))
			plt.plot(logktau,tt_evec[j,:]*tt_eval[j])
			l += 1
'''
plt.show()
