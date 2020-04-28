'''
 read_res_nvt.py
 Script that read the res_nvt.txt file (output of nvt_BB_real run, very big file)
 and extract from the file only the relevant rows
 save those rows in res_nvt_count_fourths.txt
 Usage:
 python read_res_nvt.py 

'''
file = open("res_nvt.txt", 'r')
file_count = open("res_nvt_count_fourths.txt", 'w')
res_lines = file.read().split('\n')
for l in res_lines:# lines of timestep - npt, with volume
	if 'The crossover attempt' in l:
		file_count.write(l+'\n')
file.close()
file_count.close()
		
