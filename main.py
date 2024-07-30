import numpy as np

f = open('Commands.txt','w')

baVector = np.linspace(0.1, 10, 100)

for b in baVector:
	cmdline = 'export OMP_NUM_THREADS=1; ./simDER' + ' option.txt ' + '--' + ' ba 0.0 ' + '%7.6f' % b + ' 0.0 ' + '\n'
	f.write(cmdline)
	
f.close()
