# createNPTruns.py
# Script that compiles and executes .cpp files
# Usage:
# python createNPTruns.py -i <filename> (without .cpp extension)

import sys, os, getopt, re 
import shutil
def main():
	run_dir =  '/home/student/Desktop/Shira_Michal/before_ delete/8_16'
	# 'd400000_lmps.dat', 'd800000_lmps.dat', 
	input_files = [ 'd1200000_lmps.dat', 'd1600000_lmps.dat', 'd2000000_lmps.dat']
	timeStep = 1000000

	for input_f in input_files:
		print ('----------------------------', input_f, '---',timeStep )
		############################# update the input_file, timeStep params in the file #################################
		in_file = open('{}/npt/in.npt'.format(run_dir), 'r')
		lines = in_file.readlines()
		in_file.close()
		in_file = open('{}/npt/in.npt'.format(run_dir), 'w')
		data = lines[2].split()
		data[1] = input_f
		lines[2] = ' '.join(data) + '\n'
		data = lines[-1].split()
		data[1] = str(timeStep)
		lines[-1] = ' '.join(data) + '\n'
		in_file.writelines(lines)
		in_file.close()
		###################################################################################################
		curr_dir = '{}/npt'.format(run_dir)
		cmd = '/home/student/lammps/src/lmp_serial -sf omp < in.npt'

		out_file = 'res_npt.txt'
		run(curr_dir, cmd, out_file) # do cd + run in.suffix
		###################################################################################################
		################################ save results of npt run ######################################
		inputFile_t = str('{}/npt/log.lammps'.format(run_dir))
		all_results = str('{}/npt/all_results_for_npt/'.format(run_dir, input_f))+input_f+'_log.lammps.txt'
		shutil.copy2(inputFile_t, all_results)
		for i in range(800000, 1000001, 50000):
			inputFile_t = str('{}/npt/dump/dump.reax.npt.'.format(run_dir))+str(i)
			all_results = str('{}/npt/all_results_for_npt/'.format(run_dir, i))+input_f+'_dump'+str(i)
			shutil.copy2(inputFile_t, all_results)
		break

		###################################################################################################
	print("yayy")

def run(dir, cmd, out_file):
    os.chdir(dir)
    os.system("echo do cd to " + os.getcwd())
    os.system("echo Running " + cmd)
    os.system('{} > {}'.format(cmd, out_file))
    os.system("echo -------------------")

if __name__=='__main__':
    main()
