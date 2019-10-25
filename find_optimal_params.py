# find_optimal_params.py
# Script that compiles and executes .cpp files
# Usage:
# python find_optimal_params.py -i <filename> (without .cpp extension)

import sys, os, getopt, re 

def main():
    run_dir = '/home/student/Desktop/Shira_Michal/level3_run_noPBC_2_4/run/for_yehuda_send'
    f1_vals = range(60, 200, 50)
    f2_vals = [0.25, 0.5, 1.0]
    suffixes = {'min': 'min',
                'nvt_1': 'nvt',
                'nvt_BB_real': 'nvt'}
    mult_arr = [1, 2, 0, 2]
    for f1 in f1_vals:
        for f2 in f2_vals:
            print ('------------------------------ f1, f2:', f1, f2)
            ############################# update the f1, f2 params in the file #################################
            params_file = open('{}/nvt_BB_real/Extra_Potential_Parameters.txt'.format(run_dir), 'r')
            lines = params_file.readlines()
            params_file.close()
            params_file = open('{}/nvt_BB_real/Extra_Potential_Parameters.txt'.format(run_dir), 'w')
            for line_num, (line, mult) in enumerate(zip(lines[3:], mult_arr)):
                data = line.split()
                data[2] = str(f1 * mult)
                data[3] = str(f2 * mult)
                lines[line_num+3] = ' '.join(data) + '\n'
            params_file.writelines(lines)
            params_file.close()
            ###################################################################################################
            ########################################## runs 3 levels ##########################################
            for suffix, exe_file in suffixes.items():
                curr_dir = '{}/{}'.format(run_dir, suffix)
                cmd = 'OMP_NUM_THREADS=1 /home/student/Desktop/original_lammps/lammps/src/lmp_omp -sf omp < in.{}'.format(exe_file)
                out_file = 'res_{}.txt'.format(exe_file)
                run(curr_dir, cmd, out_file) # do cd + run in.suffix
            ###################################################################################################
            ################################ save results of nvt__BB run ######################################
            # add code
            f1_f2_list = []
            f1_f2 = [] #[num_Moles, num_Specs, num_Timestep]
            num_Moles = 0 
            num_Specs = 0
            num_Timestep = 0
            inputFile_t = str('{}/nvt_BB_real/species.out'.format(run_dir))

                        
            #inputFile = open('{}/nvt_BB_real/species.out'.format(run_dir), 'r')
            file = open(inputFile_t, 'r')
            line = file.readline()
            while line:
	            # creating the reviews_list - an array of arrays that contains the full details required for each f1_f2	   
	            if re.match(r"([0-9]+)", line):
		            temp = re.match(r"([0-9])+[ ]+([0-9])+[ ]+([0-9])+", line)
		            if temp:
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
            print(f1_f2_list)

            file.close()

            #f1_f2_info_file = open(os.getcwd() + "\\result_f1=" + str(f1)+ "_f2=" + str(f2) + ".txt", 'w') 
            f1_f2_info_file_t = str('{}/nvt_BB_real/all_results_for_f1_f2/result_f1='.format(run_dir)) + str(f1)+ "_f2=" + str(f2) +".txt"
            f1_f2_info_file = open(f1_f2_info_file_t, 'w')
            str_forces = "f1 = " + str(f1) + " f2 = " + str(f2) + "\n"
            f1_f2_info_file.write(str_forces + str(f1_f2_list))
            f1_f2_info_file.close()
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
