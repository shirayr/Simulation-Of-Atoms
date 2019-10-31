import re
import os
import csv

inputFile = os.path.join(os.getcwd(),"info.txt")
file = open(inputFile, 'r')
line = file.readline()

csv_file = open(os.getcwd() + "\\result_.csv", 'w')
headers = str(line)[1:-2]
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
