'''
 Isothermal.py
 Script that read npt files
 calculate for each file the Isothermal Compression for the last 50K timestep 
 - when the volume has stabilized.
 the formula: <(v(t) - <v>)^2> / (<v> * kB *T)
 plot the gaph of the Isothermal values.
 in addition plot the Volume of % Tziluv, Avg Volume for each % Tziluv and volume module graphes
 Usage:
 python Isothermal.py 

'''
import glob
import re
from numpy import mean
import numpy as np
import collections
import matplotlib.pyplot as plt
Isoth = {} # the Isothermal for each npt-file
kB = np.array([1.38065e-23], dtype=np.longfloat) # Boltzman constant
kB = kB[0] 
T = 300 # kelvin temprature (const) of all the system
avg_volums = []
tziluv_percent = []
tziluv = -0.1
all_volume_files = []
all_densities_files = []
for file_name in glob.glob("Npt_Files\*.txt"): # reading all the npt - files
	nvtTimestep = int(re.findall(r'\d+',file_name)[0])
	file = open(file_name, 'r')
	npt_lines = file.read().split('\n')
	v_file = []
	d_file = []
	t_steps = []
	for line_npt in npt_lines[55:10056]:# lines of timestep - npt, with volume
		res = line_npt.split()
		v_file.append(float(res[-2]))
		d_file.append(float(res[-1]))
		t_steps.append(int(res[0]))
	all_volume_files.append(v_file)
	all_densities_files.append(d_file)
	npt_last_lines = npt_lines[9555:10056] # last 50K timestep when the volume ha stabilized
	tziluv = tziluv + 0.1
	tziluv_percent.append(tziluv)
	volumes = []
	for line in npt_last_lines:
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
    
    
###############################visualization of the data###################
######################## Density of % Tziluv
plt.figure()
for (d_f, p_tziluv) in zip(all_densities_files, range(len(all_densities_files))):	
	l = "{}% Tziluv".format( str(p_tziluv*10+20))
	plt.plot(t_steps, d_f, label = l) # Densities all files
plt.legend()
plt.xlabel("Timestep of NPT")
plt.ylabel("Density")
plt.title("Density of % Tziluv")
plt.savefig("Density of % Tziluv.png")
plt.show()

######################## Volume of % Tziluv
plt.figure()
for (v_f, p_tziluv) in zip(all_volume_files, range(len(all_volume_files))):	
	l = "{}% Tziluv".format( str(p_tziluv*10+20))
	plt.plot(t_steps, v_f, label = l) # volumes all files
plt.legend()
plt.xlabel("Timestep of NPT")
plt.ylabel("Volume")
plt.title("Volume of % Tziluv")
plt.savefig("Volume of % Tziluv.png")
plt.show()

######################## Avg_volume
od = collections.OrderedDict(sorted(Isoth.items()))
print(od)
axis_x = []
axis_y = []
print(len(avg_volums))
for (t, v) in zip(tziluv_percent, avg_volums):
	axis_x.append(t+0.20)
	axis_y.append(v)
plt.clf()
plt.plot(axis_x,axis_y,linestyle='--', marker='o', color='b')
plt.xlabel("Tziluv of NVT")
plt.ylabel("V")
plt.title("Avg_volume")
plt.savefig("Avg_volume.png")
plt.show()

######################## Isothermal Compression
axis_x = []
axis_y1 = []
axis_y2 = []
for (t, (k, v)) in zip(tziluv_percent, od.items()):
	axis_x.append(t+0.20)
	axis_y1.append(v)
	axis_y2.append(1/v)
plt.clf()
plt.plot(axis_x,axis_y1,linestyle='--', marker='o', color='b')
plt.xlabel("Tziluv of NVT")
plt.ylabel("Bt")
plt.title("Isothermal Compression")
plt.savefig("Isothermal Compression.png")
plt.show()

######################### Volume Module
plt.clf()
plt.plot(axis_x,axis_y2,linestyle='--', marker='o', color='b')
plt.xlabel("Tziluv of NVT")
plt.ylabel("K")
plt.title("Volume Module")
plt.savefig("Volume Module")
plt.show()