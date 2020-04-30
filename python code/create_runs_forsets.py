# create_runs_forsets.py
# Script that compiles and executes .cpp files
# Usage:
# python create_runs_forsets.py (without .cpp extension)

import sys, os, getopt, re 
import shutil
def main():
    run_dir =  '/home/student/Desktop/Shira_Michal/before_ delete/8_16_forsets'
    params_file = open('{}/nvt_BB_real/Extra_Potential_Parameters.txt'.format(run_dir), 'r')
    lines = params_file.readlines()
    f1_vals = lines[19]
    f1_vals = f1_vals.split(",")
    f2_vals = lines[21]
    f2_vals = f2_vals.split(",")
    params_file.close()
    f1_vals = range(int(f1_vals[0]), int(f1_vals[1]), int(f1_vals[2]))# range(50, 151, 50)
    f2_vals = float(f2_vals[0])# range(0.75,1,1)
	# 3 steps of run
    suffixes = [('min', 'min'),
           ('nvt_1', 'nvt'),
            ('nvt_BB_real', 'nvt')]
    for f11 in f1_vals:
        for f12 in f1_vals:
            #for f13 in f1_vals:
            for f14 in f1_vals:
                if True: #(f11 == 50 and f12 == 100 and f14 == 150): # for case of run only one combination
                    print ('------------------------------ f11, f12, f14:', f11, f12,0, f14)
                    ############################# update the f1, f2 params in the file #################################
                    params_file = open('{}/nvt_BB_real/Extra_Potential_Parameters.txt'.format(run_dir), 'r')
                    lines = params_file.readlines()
                    params_file.close()
                    params_file = open('{}/nvt_BB_real/Extra_Potential_Parameters.txt'.format(run_dir), 'w')
                    mult = []
                    mult.append(f11)
                    mult.append(f12)
                    mult.append(0)#f13)
                    mult.append(f14)
                    for line_num, (line, f1) in enumerate(zip(lines[14:], mult)):
                        data = line.split()
                        data[2] = str(f1)
                        data[3] = str(f2_vals)
                        if line_num == 2:
                            data[3] = str(0)
                        lines[line_num+14] = ' '.join(data) + '\n'
                    params_file.writelines(lines)
                    params_file.close()
                    ###################################################################################################
                    ########################################## runs 3 levels ##########################################
                    for suffix, exe_file in suffixes:
                        if True: # suffix == 'nvt_BB_real': #just for the first time you run the program on this size of run(2-4, 8-18, 16-32...) you need the True else use the 'nvt_BB_real' condition
                            curr_dir = '{}/{}'.format(run_dir, suffix)
                            # cmd = '/home/student/lammps_forsets/src/lmp_serial < in.{}'.format(exe_file) #for serial run, first compile the serial orders
                            cmd = 'OMP_NUM_THREADS=4 /home/student/lammps_forsets/src/lmp_omp -sf omp < in.{}'.format(exe_file)
                            out_file = 'res_{}.txt'.format(exe_file)
                            run(curr_dir, cmd, out_file) # do cd + run in.suffix
                    ###################################################################################################
                    ################################ save results of nvt__BB run ######################################
                    inputFile_t = str('{}/nvt_BB_real/species.out'.format(run_dir))
                    all_species = str('{}/nvt_BB_real/all_results_for_f1_f2/'.format(run_dir))+'species'+str(f11)+"_"+ str(f12)+ "_"+str(0)+ "_"+str(f14)+'.txt'
                    shutil.copy2(inputFile_t, all_species)
                    inputFile_t = str('{}/nvt_BB_real/bonds.reax'.format(run_dir))
                    all_species = str('{}/nvt_BB_real/all_results_for_f1_f2/'.format(run_dir)) + 'bonds' + str(f11) + "_" + str(f12) + "_" + str(0) + "_" + str(f14) + '.reax'
                    shutil.copy2(inputFile_t, all_species)

                    ###################################################################################################
    print("complited")

def run(dir, cmd, out_file):
	# this func execute cpp functions
    os.chdir(dir)
    os.system("echo do cd to " + os.getcwd())
    os.system("echo Running " + cmd)
    os.system('{} > {}'.format(cmd, out_file))
    os.system("echo -------------------")

if __name__=='__main__':
    main()
