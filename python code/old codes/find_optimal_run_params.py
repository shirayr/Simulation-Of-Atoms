'''
 find_optimal_run_params.py
 Script that read species files
 calculate for each run the amount of tziluv by meadian method (by 11 timestemp)
 print the good runs and the best run with the params
 Usage:
 python find_optimal_run_params.py 

'''
import sys, os, getopt
import shutil
import pandas as pd
import operator
def pandas_func(run_dir,f1_vals):
	f2 = 0.75
	dict_of_res =  {}
	with open('{}/species50_50_0_50.txt'.format(run_dir), 'r') as f:
		species_line = f.read().split('\n')
		species_value = species_line[1].split()
		delta = int(species_value[0]) # size of step of species
		last_line = species_line[-2].split()
		Timestep = int(last_line[0]) # last timestep
		firstT = int(Timestep*4/5) # firt timestep to check
		firstT = round(firstT,-len(str(delta))+1)
		step = int((Timestep - firstT)/10)
		step = round(step, -len(str(step))+1) # size of step to check
		print("firstT ",firstT)
		print("Step ",step)
	directory = '{}/CSV'.format(run_dir)
	if os.path.exists(directory):
		shutil.rmtree(directory, ignore_errors=True)
	os.makedirs(directory)
	for f11 in f1_vals:
		for f12 in f1_vals:
			#for f13 in f1_vals:
			for f14 in f1_vals:
				
				nextStep = firstT
				with open('{}/species{}_{}_0_{}.txt'.format(run_dir, f11, f12,f14), 'r') as file:	

					with open('{}/CSV/species{}_{}_0_{}.csv'.format(run_dir, f11, f12,f14),'a') as csv_file:
						headers_csv = "Step,f11,f12,f13,f14,D,E,DE,D2E,DE2,DE3,D2E3,D0E0"
						csv_file.write(headers_csv)
						species_line = file.read().split('\n')
						species_headlines = [species_line[i].split()[1:] for i in range(0,len(species_line)-1,2)]
						species_values = [species_line[i].split() for i in range(1,len(species_line),2)]
						for headers, info in zip(species_headlines, species_values):
							if int(info[0]) == nextStep: # line to check - saving the info of the step
								dict = {"step":nextStep,"f11":f11,"f12":f12,"f13":0,"f14":f14,"C11H18N2": 0, "C19H20O4": 0, "C30H38O4N2": 0, "C41H56O4N4": 0, "C49H58O8N2": 0, "C68H78O12N2": 0,"C79H96O12N4": 0, "Other": 0}
								count_others = 0
								for i, (mole_name, no_moles) in enumerate(zip(headers[3:], info[3:])):
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
								csv_file.write('\n'+str(res)[1:-1])
								nextStep = nextStep +step
				##############pandas - for each mole in each run calculate the meadian value of the 11 step
				df_csv = '{}/CSV/species{}_{}_0_{}.csv'.format(run_dir, f11, f12,f14)
				df = pd.read_csv(df_csv)
				df_runs = pd.read_csv(df_csv,  usecols=['DE','D2E','DE2','DE3','D2E3','D0E0'])
				medians =[]
				tziluv = [1,2,2,3,4,0]
				for column in df_runs:
					medians.append(df_runs[column].sort_values(ascending=True).median())
				if medians[-1] == 0.0: # others mole shuld be - 0 for good run
					goodness = sum([a*b for a,b in zip(medians,tziluv)])
					if goodness > 0: # saving the runs that tziluv was happend
						dict_of_res['f11={}_f12={}_f13=0_f14={}_f2=0.75'.format(f11, f12,f14)] = goodness
					
	print("\nruns with good parameters:")
	for keys,values in dict_of_res.items():
		print(keys.replace("_"," "))
		print("Tziluv: ",values)
	print("\nThe best parameters:")
	maxRun = max(dict_of_res.items(), key=operator.itemgetter(1))[0].replace("_","\n")
	print(maxRun)
	print("Tziluv: ",dict_of_res[max(dict_of_res.items(), key=operator.itemgetter(1))[0]])
					
				

if __name__ == "__main__":
	pandas_func('all_results_for_f1_f2 _506',range(50,151,50))
