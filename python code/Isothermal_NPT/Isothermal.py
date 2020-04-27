'''
 Isothermal.py
 Script that read npt files
 calculate for each file the Isothermal Compression for the last 50K timestep 
 - when the volume has stabilized.
 the formula: <(v(t) - <v>)^2> / (<v> * kB *T)
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
T = 300 # celvin temprature (const) of all the system
avg_volums = []
tziluv_percent = []
tziluv = -0.1
for file_name in glob.glob("Npt_Files\*.txt"): # reading all the npt - files
	nvtTimestep = int(re.findall(r'\d+',file_name)[0])
	file = open(file_name, 'r')
	npt_lines = file.read().split('\n')
	npt_lines = npt_lines[9555:10056] # last 50K timestep when the volume ha stabilized
	tziluv = tziluv + 0.1
	tziluv_percent.append(tziluv)
	volumes = []
	for line in npt_lines:
		res = line.split()
		volumes.append(float(res[-2]))
	volumes = np.array(volumes)
	file.close()
	avgV = mean(volumes) # <v>
	avg_volums.append(avgV)
	delta = (volumes - avgV) # v(t) - <v>
	delta2 = mean(delta * delta) 
	Isothermal =  delta2 / (kB * avgV * T) 
	Isoth[nvtTimestep] = Isothermal
od = collections.OrderedDict(sorted(Isoth.items()))
print(od)
axis_x = []
axis_y = []
print(len(avg_volums))
for (t, v) in zip(tziluv_percent, avg_volums):
	axis_x.append(t)
	axis_y.append(v)
plt.clf()
plt.plot(axis_x,axis_y)
plt.xlabel("Tziluv of NVT")
plt.ylabel("V")
plt.title("Avg_volume")
plt.savefig("Avg_volume.png")
plt.show()
axis_x = []
axis_y = []
for (t, (k, v)) in zip(tziluv_percent, od.items()):
	axis_x.append(t)
	axis_y.append(v)
plt.clf()
plt.plot(axis_x,axis_y)
plt.xlabel("Tziluv of NVT")
plt.ylabel("Bt")
plt.title("Isothermal Compression")
plt.savefig("Isothermal Compression.png")
plt.show()