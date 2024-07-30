import numpy as np
import matplotlib.pyplot as plt


cmd = '/home/weicheng/Desktop/seu_project/netSimulation/simShotContactClose/datafiles/simDER_shotAngle_1.000000_viscosity_0.000000.txt'
data = np.loadtxt(cmd)
plt.plot(data[:,0], data[:,1], '-')

cmd = '/home/weicheng/Desktop/seu_project/netSimulation/simShotContactClose/datafiles/simDER_shotAngle_1.000000_viscosity_0.001000.txt'
data = np.loadtxt(cmd)
plt.plot(data[:,0], data[:,1], '--')

cmd = '/home/weicheng/Desktop/seu_project/netSimulation/simShotContactClose/datafiles/simDER_shotAngle_1.000000_viscosity_0.010000.txt'
data = np.loadtxt(cmd)
plt.plot(data[:,0], data[:,1], '-.')

'''
gVector = np.logspace(3, 6, 40)
option = [1, 2, 3]
getData = np.zeros((40, 4))
'''

'''
epsilon = np.linspace(0, 0.5, 40)
option = [1]
getData = np.zeros((40, 2))
temp = 0
for e in epsilon:

	getData[temp, 0] = e

	
	cmd = '/home/weicheng/Desktop/seu_project/nonlinear_gridshell/hypderElasticRod/datafiles_rod/simDER_rod_optionType_1_epsilon_' + '%7.6f' % e + '.txt'
	data123 = np.loadtxt(cmd)
	getData[temp, 1] = data123

	cmd = '/home/weicheng/Desktop/seu_project/nonlinear_gridshell/hypderElasticRod/datafiles/simDER_rod_optionType_1_epsilon_' + '%7.6f' % e + '.txt'
	data123 = np.loadtxt(cmd)
	#getData[temp, 2] = data123

	cmd = '/home/weicheng/Desktop/seu_project/nonlinear_gridshell/hypderElasticRod/datafiles/simDER_rod_optionType_2_gVector_' + '%7.6f' % e + '.txt'
	data123 = np.loadtxt(cmd)
	getData[temp, 2] = data123
	
	cmd = '/home/weicheng/Desktop/seu_project/nonlinear_gridshell/hypderElasticRod/datafiles/simDER_rod_optionType_3_gVector_' + '%7.6f' % e + '.txt'
	data123 = np.loadtxt(cmd)
	getData[temp, 3] = data123

	
	temp = temp + 1

plt.plot(getData[:,0], getData[:,1], 'o-')
#plt.plot(getData[:,0], getData[:,2], 's-')
#plt.plot(getData[:,3], getData[:,0], '^')

#plt.xscale('log')
#plt.yscale('log')
'''
#np.savetxt('/home/weicheng/Desktop/seu_project/nonlinear_gridshell/hypderElasticRod/simDER',getData)

'''
gVector = np.linspace(0.0, 10, 40)

getData = np.zeros((40, 3))
temp = 0

for e in gVector:

	getData[temp, 0] = e

	
	cmd = '/home/weicheng/Desktop/seu_project/netSimulation/simDiscreteNet/datafiles/simDiscreteNet_gVector_' + '%7.6f' % e + '.txt'
	data123 = np.loadtxt(cmd)
	getData[temp, 1] = data123

	cmd = '/home/weicheng/Desktop/seu_project/narrowStrip/elastic_beam/datafiles/simDER_rod_gVector_' + '%7.6f' % e + '.txt'
	data123 = np.loadtxt(cmd)
	getData[temp, 2] = data123
	
	temp = temp + 1

plt.plot(getData[:,0], getData[:,1], '-')
plt.plot(getData[:,0], getData[:,2], 'o')

np.savetxt('/home/weicheng/Desktop/seu_project/netSimulation/beam.txt',getData)
'''
plt.show()