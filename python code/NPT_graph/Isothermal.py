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
files = []
kB = np.array([1.38065e-23], dtype=np.longfloat) # Boltzman constant
kB = kB[0] 
T = 300 # kelvin temprature (const) of all the system
avg_volums = []
avg_densities = []
tziluv_percent = []
tziluv = -0.1
all_volume_files = []
all_densities_files = []
for file_name in glob.glob("Npt_Files_four_M\*.txt"): # reading all the npt - files
	nvtTimestep = int(re.findall(r'\d+',file_name)[0])
	files.append(nvtTimestep)
files.sort() # sorting the file by the % Tziluv
for nvtTimestep in files: # reading all the npt - files
	file_name = "Npt_Files_four_M\d{}_lmps.dat_log.lammps.txt".format(nvtTimestep)
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
	densities = []
	for line in npt_last_lines:
		res = line.split()
		volumes.append(float(res[-2]))
		densities.append(float(res[-1]))
	volumes = np.array(volumes)
	densities = np.array(densities)
	file.close()
	avgV = mean(volumes) # <v>
	avg_volums.append(avgV)
	avg_densities.append(mean(densities))
	delta = (volumes - avgV) # v(t) - <v>
	delta2 = mean(delta * delta) 
	Isothermal =  delta2 / (kB * avgV * T) 
	Isoth[nvtTimestep] = Isothermal
###############################visualization of the data###################
######################## Density of % Tziluv
plt.figure()
for (d_f, p_tziluv) in zip(all_densities_files, range(len(all_densities_files))):	
	l = "{}% Tziluv".format( str(p_tziluv*10))
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
	l = "{}% Tziluv".format( str(p_tziluv*10))
	plt.plot(t_steps, v_f, label = l) # volumes all files
plt.legend()
plt.xlabel("Timestep of NPT")
plt.ylabel("Volume")
plt.title("Volume of % Tziluv")
plt.savefig("Volume of % Tziluv.png")
plt.show()
######################## Avg_volume
axis_x = []
axis_y = []
for (t, v) in zip(tziluv_percent, avg_volums):
	axis_x.append(t)
	axis_y.append(v)
plt.clf()
plt.plot(axis_x,axis_y,linestyle='--', marker='o', color='b')
plt.xlabel("Tziluv of NVT")
plt.ylabel("V")
plt.title("Avg_volume")
plt.savefig("Avg_volume.png")
plt.show()
######################## Avg_Density
axis_x = []
axis_y = []
for (t, d) in zip(tziluv_percent, avg_densities):
	axis_x.append(t)
	axis_y.append(d)
plt.clf()
plt.plot(axis_x,axis_y,linestyle='--', marker='o', color='b')
plt.xlabel("Tziluv of NVT")
plt.ylabel("D")
plt.title("Avg_density")
plt.savefig("Avg_density.png")
plt.show()
######################## Isothermal Compression
od = collections.OrderedDict(sorted(Isoth.items()))
axis_x = []
axis_y1 = []
axis_y2 = []
for (t, (k, v)) in zip(tziluv_percent, od.items()):
	axis_x.append(t)
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