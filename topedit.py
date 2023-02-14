import argparse
from dataclasses import dataclass


@dataclass
class AtomType:
    name: str
    mass: float
    charge: float
    ptype: str
    c6: float
    c12: float
    user: bool


@dataclass
class NonBondParam:
    i: str
    j: str
    funda: int
    level: int
    ringed: bool
    user: bool
    lj_change: bool = False # for change all lj so that the change does not happen twixe between to modified atom types


lj_levels = {
    '0.76824E-00':-1,
    '0.24145E-00': 0,
    '0.21558E-00': 1,
    '0.19402E-00': 2,
    '0.17246E-00': 3,
    '0.15091E-00': 4,
    '0.13366E-00': 5,
    '0.11642E-00': 6,
    '0.99167E-01': 7,
    '0.86233E-01': 8,
    '0.45440E-00': 9,
    '0.10620E-00': 0,
    '0.94820E-01': 1,
    '0.85338E-01': 2,
    '0.75856E-01': 3,
    '0.66375E-01': 4,
    '0.58789E-01': 5,
    '0.51203E-01': 6,
    '0.43617E-01': 7,
    '0.37928E-01': 8
}

ringed_c6 = ['0.10619E-00', '0.94820E-01', '0.85338E-01', '0.75856E-01', '0.66375E-01', '0.58789E-01', '0.51203E-01', '0.43617E-01', '0.37928E-01']

#epsilon, sigma
bond_terms = {
   -1: (5.6, 0.57),
    0: (5.6, 0.47),    
    1: (5.0, 0.47),
    2: (4.5, 0.47),
    3: (4.0, 0.47),
    4: (3.5, 0.47),
    5: (3.1, 0.47),
    6: (2.7, 0.47),
    7: (2.3, 0.47),
    8: (2.0, 0.47),
    9: (2.0, 0.62),
}

roman_num = {
    -1: 'B',
    0: 'O',
    1: 'I',
    2: 'II',
    3: 'III',
    4: 'IV',
    5: 'V',
    6: 'VI',
    7: 'VII',
    8: 'VIII',
    9: 'IX'
}


def calculate_lj_constants(nonbond):
    eps, s = bond_terms[nonbond.level]
    r = 1
    if nonbond.ringed: r, s = 0.75, 0.43
    a = 4 * r * eps * s ** 12
    b = 4 * r * eps * s ** 6
    return a, b


def standard_form(decimal):
    sciform = "{:.4E}".format(decimal)
    exp = str(int(sciform[-1])-1)
    return('0.' + sciform.replace('.', '')[:-1] + exp)


def get_bond_params(c6):
    level = lj_levels[c6]
    ringed = False
    if c6 in ringed_c6:
        ringed = True
    return level, ringed


def process_atom_type_content(atom_type_content):
    return AtomType(
        name=atom_type_content[0],
        mass=float(atom_type_content[1]),
        charge=float(atom_type_content[2]),
        ptype=atom_type_content[3],
        c6=float(atom_type_content[4]),
        c12=float(atom_type_content[5]),
        user=False
    )


def process_nonbond_param_content(nonbond_param_content):
    level, ringed = get_bond_params(nonbond_param_content[3])
    return NonBondParam(
        i=nonbond_param_content[0],
        j=nonbond_param_content[1],
        funda=int(nonbond_param_content[2]),
        level=level,
        ringed=ringed, 
        user=False
    )


def read_itp(file):
    contents = []
    with open(file, 'r') as f:
        for line in f.readlines():
            line_content = line.split(';')[0].split()
            if len(line_content) > 0: contents.append(line_content)
    return contents


def process_itp_contents(itp_contents):

    # SPLITS CONTENTS INTO SECTIONS BY LINE INDEX
    headers = {}
    for i, line_content in enumerate(itp_contents):
        if line_content[0] == '[':
            header = line_content[1]
            if header in headers.keys():
                header += '_2'
            headers[header] = i
    # ATOM TYPES
    atom_types_content = itp_contents[
        headers['atomtypes']+1 : headers[list(headers.keys())[list(headers.keys()).index('atomtypes')+1]]
    ]
    atom_types = [process_atom_type_content(content) for content in atom_types_content]
    # NON BOND PARAMS
    nonbond_params_content = itp_contents[
        headers['nonbond_params']+1 : headers[list(headers.keys())[list(headers.keys()).index('nonbond_params')+1]]
    ]
    nonbond_params = [process_nonbond_param_content(content) for content in nonbond_params_content]
    
    return atom_types, nonbond_params


def copy_atom_type(target_name, copy_name, atom_types, nonbond_params):
    copied = False
    for atom_type in atom_types:
        if atom_type.name == target_name:
            new_atom_type = AtomType(
                name=copy_name,
                mass=atom_type.mass,
                charge=atom_type.charge,
                ptype=atom_type.ptype,
                c6=atom_type.c6,
                c12=atom_type.c12,
                user=True
            )
            atom_types.append(new_atom_type)
            for nonbond in nonbond_params:
                if nonbond.i == target_name or nonbond.j == target_name:
                    i_name, j_name = nonbond.i, nonbond.j
                    if i_name == target_name: i_name = copy_name
                    if j_name == target_name: j_name = copy_name
                    nonbond_params.append(NonBondParam(
                        i=i_name,
                        j=j_name,
                        funda=nonbond.funda,
                        level=nonbond.level,
                        ringed=nonbond.ringed,
                        user=True
                    ))
            # for cross term with unmodified self
            for nonbond in nonbond_params:
                if nonbond.i == target_name and nonbond.j == target_name:
                    nonbond_params.append(NonBondParam(
                        i=copy_name,
                        j=nonbond.j,
                        funda=nonbond.funda,
                        level=nonbond.level,
                        ringed=nonbond.ringed,
                        user=True
                    ))
            copied = True
    if copied == True:
        print(' >>> SUCCESS: atom type {tg} copied to {cp}'.format(tg=target_name, cp=copy_name))
    else:
        print(' >>> ERROR: atom type {tg} does not exist'.format(tg=target_name))
    return atom_types, nonbond_params


def change_lj(target_i, target_j, new_level, nonbond_params):
    success = False
    for nonbond in nonbond_params:
        if nonbond.i == target_i and nonbond.j == target_j:
            old_level = nonbond.level
            nonbond.level = new_level
            success = True
            print(' >>> SUCCESS: LJ level of {tgi} {tgj} interaction changed from {lvlo} to {lvln}'.format(tgi=target_i, tgj=target_j, lvlo=roman_num[old_level], lvln=roman_num[new_level]))
    if success == False:
        print(' >>> ERROR: nonbond interaction between {tgi} and {tgj} does not exist'.format(tgi=target_i, tgj=target_j))
    return nonbond_params


def change_all_lj(target, change, nonbond_params):
    success = False
    for nonbond in nonbond_params:
        if nonbond.i == target or nonbond.j == target:
            if nonbond.lj_change == False:
                new_level = nonbond.level + change
                if new_level < 0:
                    new_level = 0
                elif new_level > 9:
                    new_level = 9
                nonbond.level = new_level
                success = True
                nonbond.lj_change = True
            else:
                print('     ({tgi} {tgj} interactions has already been edited)'.format(tgi=nonbond.i, tgj=nonbond.j))
    if success == False:
        print(' >>> ERROR: atom {tg} does not exist'.format(tg=target))
    else:
        print(' >>> SUCCESS: LJ level of all {tg} interactions changed by {chg}'.format(tg=target, chg=str(change)))


def write_atom_type(atom_type, file):
    file.write(' '.join([atom_type.name, str(atom_type.mass), str(atom_type.charge), atom_type.ptype, str(atom_type.c6), str(atom_type.c12)]) + '\n')


def write_nonbond_params(nonbond, file):
    a, b = calculate_lj_constants(nonbond)
    ringed = ''
    if nonbond.ringed:
        ringed = 'R'
    file.write('\t'.join([nonbond.i, nonbond.j, str(nonbond.funda), standard_form(b), standard_form(a), '; ' + ringed + roman_num[nonbond.level]]) + '\n')


def write(atom_types, nonbond_params, outdir):

    custom_atom_types = []
    default_atom_types = []
    custom_self_nonbond_params = []
    default_self_nonbond_params = []
    custom_nonbond_params = []
    default_nonbond_params = []
    for atom_type in atom_types:
        if atom_type.user == True:
            custom_atom_types.append(atom_type)
        else:
            default_atom_types.append(atom_type)
    for nonbond in nonbond_params:
        if nonbond.i == nonbond.j:
            if nonbond.user == True:
                custom_self_nonbond_params.append(nonbond)
            else:
                default_self_nonbond_params.append(nonbond)
        else:
            if nonbond.user == True:
                custom_nonbond_params.append(nonbond)
            else:
                default_nonbond_params.append(nonbond)
    
    with open(outdir, 'w') as f:

        f.write("; MODIFIED MARTINI 2.2 FORCEFIELD CREATED BY TOPEDIT.PY (FERGUS WATERHOUSE)\n")
        f.write("""


; O     - supra attractive:     (eps=5.6, s=0.47)
; I     - attractive:           (eps=5.0, s=0.47)
; II    - almost attractive:    (eps=4.5, s=0.47)
; III   - semi attractive:      (eps=4.0, s=0.47)
; IV    - intermediate:         (eps=3.5, s=0.47)
; V     - almost intermediate:  (eps=3.1, s=0.47)
; VI    - semi repulsive:       (eps=2.7, s=0.47)
; VII   - almost repulsive:     (eps=2.3, s=0.47)
; VIII  - repulsive:            (eps=2.0, s=0.47)
; IX    - super repulsive:      (eps=2.0, s=0.62)
;
; RINGS: for ring-ring interactions eps is reduced to 75%, sigma=0.43.


        """)
        f.write("\n\n[ defaults ]\n")
        f.write("1 1\n")

        # ====== ATOM TYPES ======
        f.write("\n\n[ atomtypes ]\n")
        # CUSTOM
        f.write("\n; ====== CUSTOM ATOM TYPES ======\n")
        for atom_type in custom_atom_types:
            write_atom_type(atom_type, f)
        # DEFAULT
        f.write("\n; ====== default atom types ======\n")
        for atom_type in default_atom_types:
            write_atom_type(atom_type, f)
        # ====== NONBOND PARAMS ======
        # CUSTOM
        f.write("\n\n[ nonbond_params ]\n")
        f.write("\n; ====== CUSTOM SELF TERMS ======\n")
        for nonbond in custom_self_nonbond_params:
            write_nonbond_params(nonbond, f)
        f.write("\n; ====== default self terms ======\n")
        for nonbond in default_self_nonbond_params:
            write_nonbond_params(nonbond, f)
        f.write("\n; ====== CUSTOM CROSS TERMS ======\n")
        for nonbond in custom_nonbond_params:
            write_nonbond_params(nonbond, f)
        f.write("\n; ====== default cross terms ======\n")
        for nonbond in default_nonbond_params:
            write_nonbond_params(nonbond, f)

        # ====== SOLVENT TOPOLOGY ======
        f.write("""
; WATER
[ moleculetype ]
; molname       nrexcl
  W             1
[ atoms ]
;id     type    resnr   residu  atom    cgnr    charge
 1      P4      1       W       W       1       0

; ANTIFREEZE
[ moleculetype ]
; molname        nrexcl
  WF             1
[ atoms ]
;id     type    resnr   residu  atom    cgnr    charge
 1      BP4     1       WF      WF      1       0
""")

    print(" >>> SUCCESS: modified topology file written to {fnm}".format(fnm=outdir))


def help():
    print("""
        COMMANDS =>
            COPY:           c <atom> <newatom>
            ATOMS:          a
            NONBONDS:       b (<atom> to get atom specific terms)
            EPS, S:         params
            CHANGE LJ:      lj <atomi> <atomj> <new level>
            WRITE:          w <filename>
            HELP:           h
            """)


def run(itp_file):
    contents = read_itp(itp_file)
    atom_types, nonbond_params = process_itp_contents(contents)
    quit = False
    while quit == False:
        cmd = input(' >> ').split()
        # QUIT
        if len(cmd) == 1 and cmd[0] == 'q':
            quit = True
        # BEADS
        elif len(cmd) == 1 and cmd[0] == 'a':
            for i, atom_type in enumerate(atom_types):
                print('\t'.join(['\t' + str(i)+'>', atom_type.name]))
        # BONDS
        elif cmd[0] == 'b':
            if len(cmd) == 1:
                for i, bond in enumerate(nonbond_params):
                    print('\t'.join(['\t' + str(i)+'>', bond.i, bond.j, roman_num[bond.level]]))
            elif len(cmd) == 2:
                for i, bond in enumerate(nonbond_params):
                    if bond.i == cmd[1] or bond.j == cmd[1]:
                        print('\t'.join(['\t' + str(i)+'>', bond.i, bond.j, roman_num[bond.level]]))
        # COPY
        elif len(cmd) == 3 and cmd[0] == 'c':
            atom_types, nonbond_params = copy_atom_type(cmd[1], cmd[2], atom_types, nonbond_params)
        # LJ PARAMS
        elif len(cmd) == 1 and cmd[0] == 'params':
            for nonbond in nonbond_params:
                a, b = calculate_lj_constants(nonbond)
                print('\t'.join([nonbond.i, nonbond.j, str(b), str(a)]))
        # MODIFY LJ
        elif len(cmd) == 4 and cmd[0] == 'lj':
            change_lj(cmd[1], cmd[2], int(cmd[3]), nonbond_params)
        # MODIFY ALL LJ
        elif len(cmd) == 3 and cmd[0] == 'all_lj':
            change_all_lj(cmd[1],  int(cmd[2]), nonbond_params)
        # WRITE
        elif len(cmd) == 2 and cmd[0] == 'w':
            write(atom_types, nonbond_params, cmd[1])
        # HELP
        elif len(cmd) == 1 and cmd[0] == 'h':
            help()
        # OTHERWISE
        else:
            help()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f')
    args = parser.parse_args()
    run(args.f)
