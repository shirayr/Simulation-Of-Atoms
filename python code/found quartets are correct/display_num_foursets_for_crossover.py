'''
 display_num_foursets_for_crossover.py
 Script that read the res_nvt_count_fourths.txt file (how many foursets found correct for applay the potential)
 plot the graph- Number of foursets for crossover, as function of time
 Usage:
 python display_num_foursets_for_crossover.py 

'''
import os
from matplotlib import pyplot as plt


def make_graph(axis_x, axis_y):
    plt.plot([axis_x], [axis_y], marker='o', markersize=3, color="red")
    plt.plot(axis_x, axis_y)
    plt.xlabel('timeStep')
    plt.ylabel('num of foursets')
    plt.title('Number of foursets for crossover, as function of time')
    plt.savefig('num_foursets_for_crossover.png')
    plt.show()


def find_fourset_output(file_path, str):
    timesteps = []
    num_foursets = []
    file = open(file_path, 'r')
    lines = file.read().split('\n')
    for line in lines:
        if str in line:
            print(line)
            words = line.split(' ')
            timesteps.append(int(words[2][:-1]))
            num_foursets.append(int(words[9]))

    return timesteps, num_foursets


if __name__ == '__main__':
    file_path = 'res_nvt_count_fourths.txt'
    str = 'The crossover attempt'
    timesteps, num_foursets = find_fourset_output(file_path, str)
    make_graph(timesteps, num_foursets)

