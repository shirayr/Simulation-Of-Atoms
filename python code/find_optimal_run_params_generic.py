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
import collections

'''
get a molecula
return: how many D, E, name of the mole and Tziluv
'''
def recognize_mole(mole_name):
	Epon = [19, 20, 4, 0] # C, H, O, N
	Detda = [11, 18, 0, 2] # C, H, O, N
	delta = [5, 5, 2, 0] # can be missed in a mole
	mole_name_org = mole_name
	mole_name = mole_name.replace('C', '').replace('H', ' ').replace('O', ' ').replace('N', ' ')
	amounts = mole_name.split()
	amounts = [int(i) for i in amounts]
	if mole_name_org == "H2O": # this a a water, its not consider a tear
		return 0, 0, "H2O", 0
	elif 'N' not in mole_name_org and len(amounts) == 3: # Epon not contain N
		for i in (0, 1, 2): # checking how many atoms for each base atom
			if (Epon[i] - delta[i]) <= amounts[i] <= (Epon[i] + delta[i]):
				continue
			else: # its a tear
				return 0, 0, "Other", 0			
		return 0, 1, "E1", 0
	elif 'O' not in mole_name_org and len(amounts) == 3: # Detda not contain O
		amounts.append(amounts[2]) # 0
		for i in (0, 1, 3): # checking how many atoms for each base atom
			if (Detda[i] - delta[i]) <= amounts[i] <= (Detda[i] + delta[i]):
				continue
			else: # its a tear
				return 0, 0, "Other", 0
		return 1, 0, "D1", 0
	elif len(amounts) == 4: # Detda and Epon - probably - Tziluv
		count_D = round(amounts[3]/Detda[3]) # because only D has N atom
		count_E = round(amounts[2]/Epon[2]) # because only E has O atom
		for i in (0, 1 ,2 ,3): # checking the amount of the staemed atoms
			if abs(amounts[i] - count_D*Detda[i] - count_E*Epon[i])  <= delta[i]:
				continue
			else: #there were too many missing atoms
				return 0, 0, "Other", 0					
		str_mol = "D" + str(count_D) + "E" + str(count_E) 
		Tziluv = count_D +count_E -1
		return count_D, count_E, str_mol, Tziluv
	return 0, 0, "Other", 0

'''
get the species file and the f valuesand return the amount of the Tziluv 
by the meadian method of the 11 sample of in the end of the run
'''
def analyze_func(run_dir,f1_vals):
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
			for f14 in f1_vals:
				nextStep = firstT
				headers_csv = ''
				with open('{}/species{}_{}_0_{}.txt'.format(run_dir, f11, f12,f14), 'r') as file:	
					#this section predict all the name of moles can be in the 11 timestep in order to write the header of the CSV file
					species_line = file.read().split('\n')
					species_headlines = [species_line[i].split()[1:] for i in range(0,len(species_line)-1,2)]
					species_values = [species_line[i].split() for i in range(1,len(species_line),2)]
					cols_dict = {}
					for headers, info in zip(species_headlines, species_values):
						if int(info[0]) == nextStep: # line to check - saving the info of the step
							for i, (mole_name, no_moles) in enumerate(zip(headers[3:], info[3:])):
								count_D, count_E, str_mol, Tziluv = recognize_mole(mole_name)
								if str_mol not in cols_dict:
									cols_dict[str_mol] = Tziluv
							nextStep = nextStep +step
					cols_dict = collections.OrderedDict(cols_dict)
					cols_dict = sorted(cols_dict.items(), key = lambda kv:(kv[1], kv[0]))
					
					nextStep = firstT
					with open('{}/CSV/species{}_{}_0_{}.csv'.format(run_dir, f11, f12,f14),'a') as csv_file:
						headers_csv = "Step,f11,f12,f13,f14"
						tziluv = []# saving the tziluv of each mole
						for (m_name, m_tziluv) in cols_dict: # init
							if m_name != 'H2O' and m_name != 'Other':
								headers_csv = headers_csv + ',' + m_name
								tziluv.append(m_tziluv)
						headers_csv = headers_csv + ",H2O,Other"
						tziluv.append(0)# H2O = 0
						tziluv.append(0)# Otheres = 0
						csv_file.write(headers_csv)
						# checking the last 11 timesetps 
						for headers, info in zip(species_headlines, species_values):
							if int(info[0]) == nextStep: # line to check - saving the info of the step
								dict = {"step":nextStep,"f11":f11,"f12":f12,"f13":0,"f14":f14}
								for (m_name, m_tziluv) in cols_dict:
									if m_name != 'H2O' and m_name != 'Other':
										dict[m_name] = 0
								dict['H2O'] = 0 
								dict['Other'] = 0 
								# calculate the moles of the run
								for i, (mole_name, no_moles) in enumerate(zip(headers[3:], info[3:])):
									count_D, count_E, str_mol, c_Tziluv = recognize_mole(mole_name)
									if str_mol in dict.keys():
										dict[str_mol] = no_moles

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
				u_cols = headers_csv.split(",")[5:]	
				df_runs = pd.read_csv(df_csv,  usecols = u_cols)
				medians =[]
				if( f11 == 50 and f12==100 and f14 == 150):
					print(tziluv)
				for column in df_runs:
					medians.append(df_runs[column].sort_values(ascending=True).median())
				if len(medians) >0 and medians[-1] == 0.0: # others mole shuld be - 0 for good run
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
	analyze_func('all_results_for_f1_f2 _27',range(50,151,50))
