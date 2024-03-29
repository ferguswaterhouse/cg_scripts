import subprocess
import argparse


def read_gro(gro):
    contents = []
    with open(gro, 'r') as f:
        for line in f.readlines()[2:-1]:
            content = line.split()
            content[0] = content[0].replace('1', '')
            if len(content[1]) > 3:
                name = content[1][:3]
                number = content[1][3:]
                content[1:2] = (name, number)
            contents.append(content)
    return contents


def divide_contents_into_molecule(contents):
    molecules = {}
    for content in contents:
        if content[0] not in molecules.keys():
            molecules[content[0]] = [content]
        else:
            molecules[content[0]].append(content)
    return molecules


def split_contents_into_residues(molecule_contents):
    residues = []

    first_atom_name = molecule_contents[0][1]
    molecule_length = 1
    while molecule_contents[molecule_length][1] != first_atom_name:
        molecule_length += 1
    
    residue_number, residue_length = 0, 0
    for content in molecule_contents:
        if residue_length >= molecule_length:
            residue_number += 1
            residue_length = 0
        residues.append(residue_number + 1)
        residue_length += 1

    return residues


def write_gro_file(gro_file, residues, molecule, out):
    with open(gro_file, 'r') as f:
        lines = f.readlines()

    with open(out, 'w') as w:
        number_of_molecules_parsed = 0
        for line in lines:
            if line.split()[0].replace('1', '') == molecule:
                new_line = (5-len(str(residues[number_of_molecules_parsed]))) * ' ' +  str(residues[number_of_molecules_parsed]) + line[5:]
                number_of_molecules_parsed += 1
            else:
                new_line = line
            w.write(new_line)
        

if __name__ =='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-gro')
    parser.add_argument('-mol')
    parser.add_argument('-out')
    args = parser.parse_args()
    contents = read_gro(args.gro)
    molecules = divide_contents_into_molecule(contents)
    residues = split_contents_into_residues(molecules[args.mol])
    write_gro_file(args.gro, residues, args.mol, args.out)