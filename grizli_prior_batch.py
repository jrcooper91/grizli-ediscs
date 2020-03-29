#run in batch with LDP spec prior
from grizli import multifit, utils, fitting
import numpy as np

f = open('j130132m1137_ldpz.txt', 'r')
lines = f.readlines()[1:]
f.close()
id = []
z_i = []

for line in lines: 
	a = line.split()
	id.append(int(a[0]))
	z_i.append(float(a[3]))

id = np.array(id)
z_i = np.array(z_i)

# With prior
z = np.arange(0,2,.001)
for id_i in id:
	for z_ii in z_i:
		sig = 0.007*(1+z_ii) # times (1+z)
		p_z = np.exp(-(z - z_ii)**2/(2*sig**2))/((2*np.pi)**0.5*sig) 
		p_z /= np.trapz(p_z, z)
		# Enforce finite probabilities
		p_z = np.maximum(p_z, 1.e-10)
	_res = fitting.run_all_parallel(id_i, prior=(z,p_z), zr=[0.05,2], args_file='fit_args.npy', verbose=True, group_name=str(id_i) +'withLDPPrior', get_output_data=True) 

#run in batch for photoz prior 
from grizli import multifit, utils, fitting
import numpy as np

f = open('j122756m1136pz.txt', 'r')
lines = f.readlines()[1:]
f.close()
id = []
z_i = []
z_u = []
z_d = []
for line in lines: 
	a = line.split()
	id.append(int(a[0]))
	z_i.append(float(a[3]))
	z_u.append(float(a[4]))
	z_d.append(float(a[5]))

id = np.array(id)
z_i = np.array(z_i)
z_u = np.array(z_u)
z_d = np.array(z_d)
# With prior
z = np.arange(0,2,.001)
for id_i in id:
	for z_ii in z_i:
		for z_ui in z_u:
			for z_di in z_d:
				sig = 0.5*(z_ui+z_di)*(1+z_ii) # times (1+z)
		p_z = np.exp(-(z - z_ii)**2/(2*sig**2))/((2*np.pi)**0.5*sig) 
		p_z /= np.trapz(p_z, z)
		# Enforce finite probabilities
		p_z = np.maximum(p_z, 1.e-10)
	_res = fitting.run_all_parallel(id_i, prior=(z,p_z), zr=[0.05,2], args_file='fit_args.npy', verbose=True, group_name=str(id_i) +'withpzPrior', get_output_data=True) 
