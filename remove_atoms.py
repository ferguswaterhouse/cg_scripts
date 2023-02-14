import argparse


def read_pdb(file):
    contents = []
    with open(file, 'r') as f:
        for line in f.readlines()[3:-2]:
            contents.append(line.split())
    return contents


def create_new_contents(old_contents, molecule, atoms_to_remove):
    new_contents = []
    atom_id = 1
    for i, content in enumerate(old_contents):
        print(i, len(old_contents))
        if content[3] == molecule and content[2] in atoms_to_remove:
            pass
        else:
            content[1] = str(atom_id)
            atom_id += 1
            new_contents.append(content)
    return new_contents


def write_new_contents(file, new_contents, new_file):
    with open(file, 'r') as f:
        lines = f.readlines()
        top = lines[:3]
        tail = lines[-1]
    with open(new_file, 'w') as f:
        f.write(''.join(top))
        for content in new_contents:
            spacing = [4, 11, 16, 21, 26, 38, 46, 54, 60, 66, 76]
            line = ''
            for i, data in enumerate(content):
                new_line = line + (spacing[i] - len(line) - len(data)) * ' ' + data
                line = new_line
            f.write(line+'\n')
        f.write(''.join(tail))
    print('Written to {}'.format(new_file))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb')
    parser.add_argument('-mol')
    parser.add_argument('-atom', nargs='+')
    parser.add_argument('-o')
    args = parser.parse_args()
    contents = read_pdb(args.pdb)
    new_contents = create_new_contents(contents, args.mol, args.atom)
    write_new_contents(args.pdb, new_contents, args.o)

