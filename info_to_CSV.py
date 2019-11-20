import re
import os
import csv

inputFile = os.path.join(os.getcwd(),"resultCSV.txt")
file = open(inputFile, 'r')
line = file.readline()
csv_file = open(os.getcwd() + "\\4_forces_last_step_500,000.csv", 'w')
headers = str(line)[1:-2]
headers = headers.replace("'","")
headers = headers.replace(" ","")
headers = headers.split(',')
for i in range(len(headers)):
	if headers[i] == 'C11H18N2':
		headers[i] = 'C11H18N2 - D'
	if headers[i] == 'C19H20O4':
		headers[i] = 'C19H20O4 - E'
	if headers[i] == 'C30H38O4N2':
		headers[i] = 'C30H38O4N2 - DE'
	if headers[i] == 'C41H56O4N4':
		headers[i] = 'C41H56O4N4 - D2E'
	if headers[i] == 'C49H58O8N2':
		headers[i] = 'C49H58O8N2 - DE2'
	if headers[i] == 'C68H78O12N2':
		headers[i] = 'C68H78O12N2 - DE3'
	if headers[i] == 'C79H96O12N4':
		headers[i] = 'C79H96O12N4 - D2E3'
	if headers[i] == 'Other':
		headers[i] = 'Other - D0E0'
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
