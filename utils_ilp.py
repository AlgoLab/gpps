import re


def read_matrix_tab(file_path):
    def parse_line(array):
        return [int(x) for x in array]

    if '.sim' in file_path:
        with open(file_path, 'r') as f_in:
            first_line = f_in.readline()[:-2]
            names = list(first_line.split('\t'))

            m = []
            for line in f_in:
                row = list(line[:-1].split('\t'))
                m.append(map(lambda x: float(x), row))

            return m, names

    else:
        with open(file_path, 'r') as fo:
            lines = fo.readlines()
            matrix = [parse_line(re.split("\s", line.rstrip("\n")))
                      for line in lines]

            return matrix


def expand_name(s, max_gains, max_losses):
    positive_names = [s + '+' + str(i) for i in range(max_gains)]
    negative_names = [s + '-' + str(i) for i in range(max_losses)]
    return positive_names + negative_names
