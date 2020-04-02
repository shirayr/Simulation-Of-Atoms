import os


def find_fourset_output(file_path, str):
    file = open(file_path, 'r')
    lines = file.read().split('\n')
    for line in lines:
        if str in line:
            print(line)


if __name__ == '__main__':
    file_path = 'res_nvt.txt'
    str = 'The crossover attempt'
    # str = 'fourset #'
    find_fourset_output(file_path, str)

