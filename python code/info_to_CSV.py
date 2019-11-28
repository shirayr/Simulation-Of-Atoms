# info_to_CSV.py
# Script that present the result of the program find_optimal_params 
# by CSV file and define the types of the Molecules
# Input: 
# resultCSV.txt
# Output:
# 4_forces_last_step_Timestep.csv	
# Usage:
# python info_to_CSV.py 

import re
import os
import csv

inputFile = os.path.join(os.getcwd(),"resultCSV.txt")
file = open(inputFile, 'r')
line = file.readline()
csv_file = open(os.getcwd() + "\\4_forces_last_step_506,000.csv", 'w')
headers = str(line)[1:-2]
headers = headers.replace("'","")
headers = headers.replace(" ","")
headers = headers.split(',')
for i in range(len(headers)):
	if headers[i] == 'C11H18N2':
		headers[i] = 'D' #'C11H18N2 - D'
	if headers[i] == 'C19H20O4':
		headers[i] = 'E' #'C19H20O4 - E'
	if headers[i] == 'C30H38O4N2':
		headers[i] = 'DE' #'C30H38O4N2 - DE'
	if headers[i] == 'C41H56O4N4':
		headers[i] = 'D2E' #'C41H56O4N4 - D2E'
	if headers[i] == 'C49H58O8N2':
		headers[i] = 'DE2' #'C49H58O8N2 - DE2'
	if headers[i] == 'C68H78O12N2':
		headers[i] = 'DE3' #'C68H78O12N2 - DE3'
	if headers[i] == 'C79H96O12N4':
		headers[i] = 'D2E3' #'C79H96O12N4 - D2E3'
	if headers[i] == 'Other':
		headers[i] = 'D0E0' #'Other - D0E0'
headers = str(headers)[1:-2]
headers = headers.replace("'","")


csv_file.write(headers+"\n")
line = file.readline()
while line:
	data = str(line)[1:-2]
	data = data.replace("'","")
	csv_file.write(data+"\n")
	line = file.readline()
	line = file.readline()

csv_file.close()
