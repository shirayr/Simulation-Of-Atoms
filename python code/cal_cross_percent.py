import matplotlib.pyplot as plt


def make_graph(axis_x, axis_y, xlabel="x", ylabel="y", title="Graph"):
    plt.plot([axis_x], [axis_y], marker='o', markersize=3, color="red")
    plt.plot(axis_x, axis_y)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()


def cross_percent_graph():
    species_file = open('species.out', 'r')
    species_line = species_file.read().split('\n')
    run_size = int(species_line[-2].split()[0])  # get the number of timestep from the last line
    step = int(species_line[1].split()[0])  # get the step intervals that dumped to the file

    sample_size = int(run_size * 0.00025)  # we took the sample size to be 0.025 percent of the run size
    partial_sample = [int(run_size / sample_size) * i for i in range(1, sample_size+1)]
    indexes = [int(i/step * 2 - 1) for i in partial_sample]
    species_headlines = [species_line[i - 1].split()[1:] for i in indexes]
    species_values = [species_line[i].split() for i in indexes]
    print(partial_sample)
    print(species_headlines)
    print(species_values)

    mole_type_to_num_of_cross = {'C11H18N2': 0,
                                 'C19H20O4': 0,
                                 'C30H38O4N2': 1,
                                 'C41H56O4N4': 2,
                                 'C49H58O8N2': 2,
                                 'C68H78O12N2': 3,
                                 'C79H96O12N4': 4}

    cross_percent = {}
    for species_headline, species_value in zip(species_headlines, species_values):
        cross_percent[species_value[0]] = 0
        for mole_type, num_mole in zip(species_headline[3:], species_value[3:]):
            cross_percent[species_value[0]] += mole_type_to_num_of_cross[mole_type] * int(num_mole) if mole_type in mole_type_to_num_of_cross else 0

    cross_percent = list(cross_percent.values())
    print(cross_percent)
    make_graph(partial_sample, cross_percent, 'Time Step', 'Cross Percent', 'Cross Percent as a function of time')


if __name__ == '__main__':
    cross_percent_graph()
