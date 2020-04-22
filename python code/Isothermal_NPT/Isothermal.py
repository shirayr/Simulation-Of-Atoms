'''
 Isothermal.py
 Script that read npt files
 calculate for each file the Isothermal Compression for the last 50K timestep 
 - when the volume has stabilized.
 the formula: <(v(t) - <v>)^2> / (<v> * kB)
 plot the gaph of the Isothermal values.
 Usage:
 python Isothermal.py 

'''
import glob
import re
from numpy import mean
import numpy as np
import collections
from matplotlib import pyplot as plt
Isoth = {} # the Isothermal for each npt-file
kB = np.array([1.38065e-23], dtype=np.longfloat) # Boltzman constant
kB = kB[0] 
for file_name in glob.glob("Npt_Files\*.txt"): # reading all the npt - files
	nvtTimestep = int(re.findall(r'\d+',file_name)[0])
	file = open(file_name, 'r')
	npt_lines = file.read().split('\n')
	npt_lines = npt_lines[9555:10056] # last 50K timestep when the volume ha stabilized
	volumes = []
	for line in npt_lines:
		res = line.split()
		volumes.append(float(res[-2]))
	volumes = np.array(volumes)
	file.close()
	avgV = mean(volumes) # <v>
	delta = (volumes - avgV) # v(t) - <v>
	delta2 = mean(delta * delta) 
	Isothermal =  delta2 / (kB * avgV) 
	Isoth[nvtTimestep] = Isothermal
od = collections.OrderedDict(sorted(Isoth.items()))
print(od)
axis_x = []
axis_y = []
for k, v in od.items():
	axis_x.append(k)
	axis_y.append(v)
	
plt.plot(axis_x,axis_y)
plt.xlabel("Timestep of NVT")
plt.ylabel("Bt")
plt.title("Isothermal Compression")
plt.savefig("Isothermal Compression.png")
plt.show()