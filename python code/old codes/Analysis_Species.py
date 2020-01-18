'''
Analtsis_Species.py
Script that Goes through the species.out file and
analyzes how many time steps have been for each num_Moles&num_Specs (the atom types and their amount)
Usage:
python Analtsis_Species.py
The #'s are for exeute on the colleage computere
'''
import re
import os
import csv

f1 = 0.5
f2 = 0.25
f1_f2_list = []
f1_f2 = [] #[num_Moles, num_Specs, num_Timestep]
num_Moles = 0 
num_Specs = 0
num_Timestep = 0
inputFile = os.path.join(os.getcwd(),"species.out")
file = open(inputFile, 'r')
#inputFile_t = str('{}/nvt_BB_real/species.out'.format(run_dir))
#file = open(inputFile_t, 'r')
line = file.readline()
headers = ''
while line:
	# creating the reviews_list - an array of arrays that contains the full details required for each f1_f2	   
	if re.match(r"#", line):
		headers = line
		line = file.readline()
	if re.match(r"([0-9]+)", line):
		temp = re.match(r"([0-9]+)[ ]+([0-9]+)[ ]+([0-9]+)", line)
		if temp:
			if int(temp.group(1)) == 152000:
				dict = {"f1":f1,"f2":f2,"C11H18N2": 0, "C19H20O4": 0, "C30H38O4N2": 0, "C41H56O4N4": 0, "C49H58O8N2": 0, "C68H78O12N2": 0,"C79H96O12N4": 0, "Other": 0}
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

				with open('output.csv', 'w') as output:
					writer = csv.writer(output)
					for key, value in dict.items():
						writer.writerow([key, str(value)])
				
			if num_Moles == 0 and num_Specs == 0 and num_Timestep == 0:
				num_Moles = int(temp.group(2))
				num_Specs = int(temp.group(3))
				num_Timestep = 1
				
				
			elif num_Moles == int(temp.group(2)) and num_Specs == int(temp.group(3)):
				num_Timestep = num_Timestep + 1
				
			elif num_Moles != int(temp.group(2)) or num_Specs != int(temp.group(3)):
				f1_f2.extend([num_Moles, num_Specs,num_Timestep])	
				f1_f2_list.append(f1_f2)
				f1_f2 = []
				num_Moles = int(temp.group(2))
				num_Specs = int(temp.group(3))
				num_Timestep = 1
	line = file.readline()
f1_f2.extend([num_Moles, num_Specs,num_Timestep])	
f1_f2_list.append(f1_f2)
print(f1_f2_list)

file.close()

f1_f2_info_file = open(os.getcwd() + "\\result_f1=" + str(f1)+ "_f2=" + str(f2) + ".txt", 'w')
#f1_f2_info_file_t = str('{}/nvt_BB_real/all_results_for_f1_f2/result_f1='.format(run_dir)) + str(f1)+ "_f2=" + str(f2) +".txt"
#f1_f2_info_file = open(f1_f2_info_file_t, 'w')
str_forces = "f1 = " + str(f1) + " f2 = " + str(f2) + "\n"
f1_f2_info_file.write(str_forces + str(f1_f2_list))
f1_f2_info_file.close()
