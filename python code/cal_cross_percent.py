import matplotlib.pyplot as plt


def make_graph(axis_x, axis_y, xlabel="x", ylabel="y", title="Graph"):
    plt.plot([axis_x], [axis_y], marker='o', markersize=5, color="red")
    plt.plot(axis_x, axis_y)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()


def cross_percent_graph(run_size):
    species_file = open('species.out', 'r')
    species_line = species_file.read().split('\n')

    five_samples = [int(run_size / 5) * i for i in range(1, 6)]
    indexes = [int(i/100 * 2 - 1) for i in five_samples]
    species_headlines = [species_line[i - 1].split()[1:] for i in indexes]
    species_values = [species_line[i].split() for i in indexes]
    print(five_samples)
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
    make_graph(five_samples, cross_percent, 'Time Step', 'Cross Percent', 'Cross Percent as a function of time')


if __name__ == '__main__':
    cross_percent_graph(2000000)
