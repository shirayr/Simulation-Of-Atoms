# find_optimal_params.py
# Script that compiles and executes .cpp files
# Usage:
# python find_optimal_params.py -i <filename> (without .cpp extension)

import sys, os, getopt, re 

def main():
    run_dir = '/home/student/Desktop/Shira_Michal/level3_run_PBC_2_4/run/for_yehuda_send'
    f1_vals = range(50, 151, 50)# range(60, 301, 30)
    f2_vals = 0.75#[0.25, 0.5, 1.0]
    suffixes = {'min': 'min',
                'nvt_1': 'nvt',
                'nvt_BB_real': 'nvt'}
    mult_arr = [1, 2, 0, 2]
    for f11 in f1_vals:
        for f12 in f1_vals:
            #for f13 in f1_vals:
            for f14 in f1_vals:
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
                for suffix, exe_file in suffixes.items():
                    curr_dir = '{}/{}'.format(run_dir, suffix)
                    cmd = '/home/student/Desktop/original_lammps/lammps/src/lmp_serial < in.{}'.format(exe_file)
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
                    if re.match(r"#", line):
                        headers = line
                        line = file.readline()
                    if re.match(r"([0-9]+)", line):
                        temp = re.match(r"([0-9]+)[ ]+([0-9]+)[ ]+([0-9]+)", line)
                        if temp:
                            if int(temp.group(1)) == 152000:
                                dict = {"f11":f11,"f12":f12,"f13":0,"f14":f14,"C11H18N2": 0, "C19H20O4": 0, "C30H38O4N2": 0, "C41H56O4N4": 0, "C49H58O8N2": 0, "C68H78O12N2": 0,"C79H96O12N4": 0, "Other": 0}
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

                                all_info_file_t = str('{}/nvt_BB_real/all_results_for_f1_f2/resultCSV.txt'.format(run_dir))
                                all_info_file = open(all_info_file_t, 'a+')
                                list_keys = []
                                list_values = []
                                for key, value in dict.items():
                                    list_keys.append(key)
                                    list_values.append(value)
                                all_info_file.write(str(list_keys)+"\n")
                                all_info_file.write(str(list_values)+"\n")

                                all_info_file.close()


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

                #f1_f2_info_file = open(os.getcwd() + "\\result_f1=" + str(f11)+"_"+ str(f12)+ "_"+str(f13)+ "_f2=" + str(f2) + ".txt", 'w')
                f1_f2_info_file_t = str('{}/nvt_BB_real/all_results_for_f1_f2/result_f1_f2.txt'.format(run_dir))
                f1_f2_info_file = open(f1_f2_info_file_t, 'a+')
                str_forces = "f1 = "+ str(f11)+"_"+ str(f12)+ "_"+str(0)+ "_"+str(f14)+ " f2 = " + str(f2_vals) + "\n"
                f1_f2_info_file.write(str_forces + str(f1_f2_list)+ "\n")
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
