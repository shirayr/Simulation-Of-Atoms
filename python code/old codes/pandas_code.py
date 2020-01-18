import sys, os, getopt, re 
import shutil
import pandas as pd
def pandas_func():
	run_dir = 'all_results_for_f1_f2 _154'
	f1_vals = range(50,151,50)
	f2 = 0.75
	with open('{}/species50_50_0_50.txt'.format(run_dir), 'r') as f:
		lines = f.read().splitlines()
		second_line = lines[1]
		delta = re.match(r"([0-9]+)[ ]+([0-9]+)[ ]+([0-9]+)", second_line)
		delta = int(delta.group(1))
		last_line = lines[-1]
		temp = re.match(r"([0-9]+)[ ]+([0-9]+)[ ]+([0-9]+)", last_line)
		Timestep = int(temp.group(1))
		firstT = int(Timestep*4/5)
		print(firstT)
		firstT = round(firstT,-len(str(delta))+1)
		step = int((Timestep - firstT)/10)
		step = round(step, -len(str(step))+1)
		print(firstT)
		print(step)
	directory = '{}/CSV'.format(run_dir)
	if not os.path.exists(directory):
		os.makedirs(directory)
	good = 0
	prams = ''
	k=1
	for f11 in f1_vals:
		for f12 in f1_vals:
			#for f13 in f1_vals:
			for f14 in f1_vals:
				
				nextStep = firstT
				with open('{}/species{}_{}_0_{}.txt'.format(run_dir, f11, f12,f14), 'r') as file:	
					
					f1_f2_list = []
					f1_f2 = [] #[num_Moles, num_Specs, num_Timestep]
					num_Moles = 0 
					num_Specs = 0
					num_Timestep = 0
					csv_file = open('{}/CSV/species{}_{}_0_{}.csv'.format(run_dir, f11, f12,f14),"r+")
					headers = "Step,f11,f12,f13,f14,D,E,DE,D2E,DE2,DE3,D2E3,D0E0"
					csv_file.write(headers)
					line = file.readline()
					while line:
						# an array of arrays that contains the full details required for each f1_f2
						if re.match(r"#", line):
							headers = line
							line = file.readline()
						if re.match(r"([0-9]+)", line):
							temp = re.match(r"([0-9]+)[ ]+([0-9]+)[ ]+([0-9]+)", line)
							if temp:
								if int(temp.group(1)) == nextStep:
									print(nextStep)
									k = k + 1
									dict = {"step":nextStep,"f11":f11,"f12":f12,"f13":0,"f14":f14,"C11H18N2": 0, "C19H20O4": 0, "C30H38O4N2": 0, "C41H56O4N4": 0, "C49H58O8N2": 0, "C68H78O12N2": 0,"C79H96O12N4": 0, "Other": 0}
									headers = re.sub(' +', ' ', headers)
									headers = re.sub('\t+', ' ', headers)
									headers = headers.split()
									headers = headers[4:]
									info = line
									info = re.sub(' +', ' ', info)
									info = re.sub('\t+', ' ', info)
									info = re.sub(' +', ' ', info)
									info = info.split()
									info = info[3:]
									count_others = 0
									for i, (mole_name, no_moles) in enumerate(zip(headers, info)):
										if mole_name in dict.keys():
											dict[mole_name] = no_moles
										else:
											count_others = count_others + int(no_moles)
									dict["Other"] = count_others
									list_keys = []
									list_values = []
									for key, value in dict.items():
										list_keys.append(key)
										list_values.append(value)
									
									res = [int(i) for i in list(dict.values())]
									#print(str(res)[1:-1])
									csv_file.write('\n'+str(res)[1:-1])
									nextStep = nextStep +step
						line = file.readline()
				##############pandas
				df_csv = '{}/CSV/species{}_{}_0_{}.csv'.format(run_dir, f11, f12,f14)
				df = pd.read_csv(df_csv)
				df_runs = pd.read_csv(df_csv,  usecols=['DE','D2E','DE2','DE3','D2E3','D0E0'])
				medians =[]
				tziluv = [1.0,2.0,2.0,3.0,4.0,0.0]
				for column in df_runs:
					medians.append(df_runs[column].sort_values(ascending=True).median())
					#print((df_runs[column]))
				if medians[-1] == 0.0:
					# print(df_csv)
					# print(df_runs)
					# print(sum([a*b for a,b in zip(medians,tziluv)]))
					goodness = sum([a*b for a,b in zip(medians,tziluv)])
					if(goodness>good):
						good = goodness
						prams = 'f11={}_f12={}_f13=0_f14={}__f2=0.75'.format(f11, f12,f14)
	print(prams)
					
				

if __name__ == "__main__":
	pandas_func()
