'''
 cal_cross_number.py
 This script reads species files,
 find and shows gragh of the crossover number for some points during the run. 
 Usage:
 python cal_cross_number.py 

'''

import matplotlib.pyplot as plt


def make_graph(axis_x, axis_y, xlabel="x", ylabel="y", title="Graph", name="Graph.png"):
    plt.plot([axis_x], [axis_y], marker='o', markersize=3, color="red")
    plt.plot(axis_x, axis_y)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(name)
    plt.show()
    

'''
get a molecula
return what the Tziluv for this molecula
'''
def recognize_mole(mole_name):
	Epon = [19, 20, 4, 0] # C, H, O, N
	Detda = [11, 18, 0, 2] # C, H, O, N
	delta = [5, 5, 2, 0] # can be missed in a mole
	mole_name_org = mole_name
	mole_name = mole_name.replace('C', '').replace('H', ' ').replace('O', ' ').replace('N', ' ')
	amounts = mole_name.split()
	amounts = [int(i) for i in amounts]
	if mole_name_org == "H2O": # this a a water, its not consider a rupture
		return 0
	elif 'N' not in mole_name_org and len(amounts) == 3: # Epon not contain N
		for i in (0, 1, 2): # checking how many atoms for each base atom
			if (Epon[i] - delta[i]) <= amounts[i] <= (Epon[i] + delta[i]):
				continue
			else: # its a rupture
				return 0	
		return 0
	elif 'O' not in mole_name_org and len(amounts) == 3: # Detda not contain O
		amounts.append(amounts[2]) # 0
		for i in (0, 1, 3): # checking how many atoms for each base atom
			if (Detda[i] - delta[i]) <= amounts[i] <= (Detda[i] + delta[i]):
				continue
			else: # its a rupture
				return 0
		return 0
	elif len(amounts) == 4: # Detda and Epon - probably - Tziluv
		count_D = round(amounts[3]/Detda[3]) # because only D has N atom
		count_E = round(amounts[2]/Epon[2]) # because only E has O atom
		for i in (0, 1 ,2 ,3): # checking the amount of the staemed atoms
			if abs(amounts[i] - count_D*Detda[i] - count_E*Epon[i])  <= delta[i]:
				continue
			else: #there were too many missing atoms
				return 0					
		str_mol = "D" + str(count_D) + "E" + str(count_E) 
		Tziluv = count_D +count_E -1
		return Tziluv
	return 0


def cross_percent_graph():
    species_file = open('species.out', 'r')
    species_line = species_file.read().split('\n')
    run_size = int(species_line[-2].split()[0])  # get the number of timestep from the last line
    step = int(species_line[1].split()[0])  # get the step intervals that dumped to the file

    sample_size = int(run_size * 0.00025)  # we took the sample size to be 0.025 percent of the run size
    partial_sample = [int(run_size / sample_size) * i for i in range(1, sample_size+1)]
    indexes = [int(i/step * 2 - 1) for i in partial_sample]
    for i in range(len(indexes)):
        if indexes[i]%2 == 0:
            indexes[i] = indexes[i]+1
    species_headlines = [species_line[i - 1].split()[1:] for i in indexes]
    species_values = [species_line[i].split() for i in indexes]
    print('indexes', indexes)
    print('partial_sample', partial_sample)
    print('species_headlines', species_headlines)
    print('species_values', species_values)

    cross_percent = {}
    for species_headline, species_value in zip(species_headlines, species_values):
        cross_percent[species_value[0]] = 0
        for mole_type, num_mole in zip(species_headline[3:], species_value[3:]):
            cross_percent[species_value[0]] += recognize_mole(mole_type) * int(num_mole)

    cross_percent = list(cross_percent.values())
    print(cross_percent)
    make_graph(partial_sample, cross_percent, 'Time Step', 'Cross Number', 'Cross Number as function of time', 'cross_number.png')


if __name__ == '__main__':
    cross_percent_graph()
