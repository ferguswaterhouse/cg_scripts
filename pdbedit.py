import argparse

def read_file(filename):
    with open(filename) as f:
        return f.readlines()


def replace_atom_type_with_water(lines, atom_types):
    new_pdb_lines = []
    last_w_index = 0
    number_of_changes = 0
    for line in lines:
        line_data = line.split()
        
        if len(line_data) == 11 and line_data[3] == 'W':
            last_w_index = int(line_data[4])
        
        if len(line_data) == 11 and line_data[3] in atom_types:
            last_w_index += 1
            res_number = line_data[1]
            x, y, z = line_data[5], line_data[6], line_data[7]
            new_data = ['ATOM', res_number, 'W', 'W', str(last_w_index), x, y, z, '1.00', '0.00', 'WAT']
            new_pdb_lines.append(new_data)
            number_of_changes += 1
        else:
            new_pdb_lines.append(line_data)
    return new_pdb_lines, number_of_changes


def write_new_pdb(new_lines):
    with open('modified.pdb', 'w') as f:
        for line in new_lines:
            if len(line) == 11:
                terms = []
                terms.append(line[0])
                terms.append((7-len(line[1])) * ' ' + line[1])
                terms.append((5-len(line[2])) * ' ' + line[2])
                terms.append((5-len(line[3])) * ' ' + line[3])
                terms.append((5-len(line[4])) * ' ' + line[4])
                terms.append((12-len(line[5])) * ' ' + line[5])
                terms.append((8-len(line[6])) * ' ' + line[6])
                terms.append((8-len(line[7])) * ' ' + line[7])
                terms.append((6-len(line[8])) * ' ' + line[8])
                terms.append((6-len(line[9])) * ' ' + line[9])
                terms.append(6 * ' ')
                terms.append(line[10])
                f.write(''.join(terms))
                f.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f')
    args = parser.parse_args()
    file_lines = read_file(args.f)
    new_lines, number_of_changes = replace_atom_type_with_water(file_lines, ('CA', 'NA', 'CL'))
    write_new_pdb(new_lines)
    print('Modified pdb modified.pdb created with {n} changes made'.format(n=number_of_changes))

