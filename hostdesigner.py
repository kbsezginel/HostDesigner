import os

from tabulate import tabulate


def hdo_read(hdo_dir, num='all', table=True, hdo='out_1.hdo'):
    """
    Read hdo output files, tabulate results and return structure information.
    """
    hdo_path = os.path.join(hdo_dir, hdo)
    hdo_file = open(hdo_path, 'r')
    hdo_lines = hdo_file.readlines()
    hdo_file.close()

    num_structures = 0
    hdo_structures = {'rmsd': [], 'energy': [], 'num_atoms': [], 'atom': [], 'coor': [],
                      'info': [], 'xyz': [], 'host': [], 'linker': [], 'index': []}
    for line_index, line in enumerate(hdo_lines):
        if line.replace(' ', '').strip('\n').isdigit():
            hdo_structures['num_atoms'].append(int(line.replace(' ', '').strip('\n')))
        if 'RMSD' in line:
            # hdo_structures['rmsd'].append(float(line.split(',')[0].split('=')[-1]))
            hdo_structures['rmsd'].append(line.split(',')[0].split('=')[-1])
            hdo_structures['energy'].append(float(line.split(',')[2].split('_')[1]))
            hdo_structures['info'].append(line.strip('_'))
            hdo_structures['linker'].append(line.strip('_').split('_')[2])
            hdo_structures['index'].append(num_structures)
            num_structures += 1

    if num == 'all':
        num = num_structures

    start_line = 2
    for structure in range(min(num, num_structures)):
        host_lines = hdo_structures['info'][structure]
        host_lines += ' ' + str(hdo_structures['num_atoms'][structure]) + '\t1\n'

        xyz_lines = str(hdo_structures['num_atoms'][structure]) + '\n'
        xyz_lines += str(structure + 1) + '_' + str(hdo_structures['rmsd'][structure]) + '\n'

        hdo_structures['atom'].append([])
        hdo_structures['coor'].append([])

        end_line = start_line + hdo_structures['num_atoms'][structure]
        for line_index, line in enumerate(hdo_lines[start_line:end_line]):
            atom_name = line.split()[0]
            new_line = '  {:2} {:3}{}'.format(atom_name, line_index + 1, line[3:])
            host_lines += new_line

            name, x, y, z = line.split()[0], line.split()[1], line.split()[2], line.split()[3]
            xyz_lines += name + '\t' + x + '\t' + y + '\t' + z + '\n'

            hdo_structures['atom'][structure].append(name)
            hdo_structures['coor'][structure].append([float(x), float(y), float(z)])

        hdo_structures['host'].append(host_lines)
        hdo_structures['xyz'].append(xyz_lines)
        start_line += hdo_structures['num_atoms'][structure] + 2

    if table:
        print(tabulate({'Structure': hdo_structures['index'][:num],
                        'Linker': hdo_structures['linker'][:num],
                        'Num Atoms': hdo_structures['num_atoms'][:num],
                        'RMSD': hdo_structures['rmsd'][:num],
                        'Energy': hdo_structures['energy'][:num]}, headers="keys"))

    return hdo_structures


def hdo_linker(hdo_structures, linkers, table=True, export=[False, os.getcwd(), 'xyz']):
    linker_structures = {'index': [], 'linker': [], 'rmsd': [], 'energy': [], 'num_atoms': [], 'atom': [], 'coor': []}
    for structure_index, hdo_linker in enumerate(hdo_structures['linker']):
        for linker in linkers:
            if hdo_linker == linker:
                linker_structures['index'].append(hdo_structures['index'][structure_index])
                linker_structures['rmsd'].append(hdo_structures['rmsd'][structure_index])
                linker_structures['energy'].append(hdo_structures['energy'][structure_index])
                linker_structures['num_atoms'].append(hdo_structures['num_atoms'][structure_index])
                linker_structures['linker'].append(hdo_structures['linker'][structure_index])
                linker_structures['atom'].append(hdo_structures['atom'][structure_index])
                linker_structures['coor'].append(hdo_structures['coor'][structure_index])
                if export[0]:
                    hdo_export(hdo_structures, structure_index, export[1], name=linker, exp_format=export[2])

    if table:
        print(tabulate({'Structure': linker_structures['index'],
                        'Linker': linker_structures['linker'],
                        'Num Atoms': linker_structures['num_atoms'],
                        'RMSD': linker_structures['rmsd'],
                        'Energy': linker_structures['energy']}, headers="keys"))

    return linker_structures


def drive_read(drive_dir, num='all', table=True, drive='out_testa.hdo'):
    """
    Reads test_drive output files, tabulate results and return structure information.
    """
    drive_path = os.path.join(drive_dir, drive)
    drive_file = open(drive_path, 'r')
    drive_lines = drive_file.readlines()
    drive_file.close()

    # Read number of atoms and drive information
    num_structures = 0
    drive_structures = {'drive': [], 'num_atoms': [], 'xyz': [], 'host': [],
                        'index': [], 'atom': [], 'coor': []}
    for line_index, line in enumerate(drive_lines):
        if line.replace(' ', '').strip('\n').isdigit():
            drive_structures['num_atoms'].append(int(line.replace(' ', '').strip('\n')))
        if 'Drive' in line:
            drive_structures['drive'].append(line)
            drive_structures['index'].append(num_structures)
            num_structures += 1

    if num == 'all':
        num = num_structures

    # Write xyz informationm
    start_line = 2
    for structure in range(min(num, num_structures)):
        host_lines = drive_structures['drive'][structure]
        host_lines += ' ' + str(drive_structures['num_atoms'][structure]) + '   1\n'

        xyz_lines = str(drive_structures['num_atoms'][structure]) + '\n'
        xyz_lines += str(structure + 1) + '_' + drive_structures['drive'][structure]

        drive_structures['atom'].append([])
        drive_structures['coor'].append([])

        end_line = start_line + drive_structures['num_atoms'][structure]
        for line_index, line in enumerate(drive_lines[start_line:end_line]):
            atom_name = line.split()[0]
            new_line = '  {:2} {:3}{}'.format(atom_name, line_index + 1, line[3:])
            host_lines += new_line

            name, x, y, z = line.split()[0], line.split()[1], line.split()[2], line.split()[3]
            xyz_lines += name + '\t' + x + '\t' + y + '\t' + z + '\n'

            drive_structures['atom'][structure].append(name)
            drive_structures['coor'][structure].append([float(x), float(y), float(z)])

        start_line += drive_structures['num_atoms'][structure] + 2

        drive_structures['host'].append(host_lines)
        drive_structures['xyz'].append(xyz_lines)

    if table:
        print(tabulate({'Structure': drive_structures['index'][:num],
                        'Drive': drive_structures['drive'][:num],
                        'Num Atoms': drive_structures['num_atoms'][:num]}, headers="keys"))

    return drive_structures


def host_export(host_dir, export_dir, host='hosta'):
    """
    Convert host structure file to xyz format.
    """
    host_path = os.path.join(host_dir, host)
    host_file = open(host_path, 'r')
    host_lines = host_file.readlines()
    host_file.close()
    num_atoms = host_lines[1].split()[0]

    xyz_lines = ''
    for line in host_lines[2:2 + int(num_atoms)]:
        name = line.split()[0]
        x, y, z = line.split()[2], line.split()[3], line.split()[4]
        xyz_lines += name + '\t' + x + '\t' + y + '\t' + z + '\n'

    host_name = os.path.basename(host_path)
    xyz_path = os.path.join(export_dir, host_name + '.xyz')
    xyz_file = open(xyz_path, 'w')
    xyz_file.write(num_atoms + '\n')
    xyz_file.write(host_name + '\n')
    xyz_file.write(xyz_lines)
    xyz_file.close()


def linker_export(linker_path, export_dir, num=5, xyz_name=''):
    """
    Export a given number of structures with lowest RMSD values from linker.out file.
    """
    linker_file = open(linker_path, 'r')
    linker_lines = linker_file.readlines()
    linker_file.close()

    num_structures = 0
    linker_structures = {'rmsd': [], 'energy': [], 'num_atoms': []}
    for line_index, line in enumerate(linker_lines):
        if line.replace(' ', '').strip('\n').isdigit():
            linker_structures['num_atoms'].append(int(line.replace(' ', '').strip('\n')))
        if 'RMSD' in line:
            linker_structures['rmsd'].append(float(line.split(',')[0].split('=')[-1]))
            linker_structures['energy'].append(float(line.split(',')[2].split('_')[1]))
            num_structures += 1

    start_line = 2
    for structure in range(min(num, num_structures)):
        xyz_lines = str(linker_structures['num_atoms'][structure]) + '\n'
        xyz_lines += str(structure + 1) + '_' + str(linker_structures['rmsd'][structure]) + '\n'

        end_line = start_line + linker_structures['num_atoms'][structure]
        for line in linker_lines[start_line:end_line]:
            name, x, y, z = line.split()[0], line.split()[1], line.split()[2], line.split()[3]
            xyz_lines += name + '\t' + x + '\t' + y + '\t' + z + '\n'

        xyz_path = os.path.join(export_dir, str(structure + 1) + xyz_name + '.xyz')
        xyz_file = open(xyz_path, 'w')
        xyz_file.write(xyz_lines)
        xyz_file.close()

        start_line += linker_structures['num_atoms'][structure] + 2


# Change to read_drive similar to read_hdo
def drive_export(drive_path, export_dir, num=5, xyz_name=''):
    """
    Export a given number of structures with lowest RMSD values from drive.out file.
    """
    drive_file = open(drive_path, 'r')
    drive_lines = drive_file.readlines()
    drive_file.close()

    num_structures = 0
    drive_structures = {'drive': [], 'num_atoms': []}
    for line_index, line in enumerate(drive_lines):
        if line.replace(' ', '').strip('\n').isdigit():
            drive_structures['num_atoms'].append(int(line.replace(' ', '').strip('\n')))
        if 'Drive' in line:
            drive_structures['drive'].append(line)
            num_structures += 1

    start_line = 2
    for structure in range(min(num, num_structures)):
        xyz_lines = str(drive_structures['num_atoms'][structure]) + '\n'
        xyz_lines += str(structure + 1) + '_' + drive_structures['drive'][structure]

        end_line = start_line + drive_structures['num_atoms'][structure]
        for line in drive_lines[start_line:end_line]:
            name, x, y, z = line.split()[0], line.split()[1], line.split()[2], line.split()[3]
            xyz_lines += name + '\t' + x + '\t' + y + '\t' + z + '\n'

        xyz_path = os.path.join(export_dir, str(structure + 1) + xyz_name + '.xyz')
        xyz_file = open(xyz_path, 'w')
        xyz_file.write(xyz_lines)
        xyz_file.close()

        start_line += drive_structures['num_atoms'][structure] + 2


def hdo_export(hdo_structures, hdo_index, export_dir, name='hosta', exp_format='xyz'):

    if exp_format == 'xyz':
        xyz_path = os.path.join(export_dir, str(hdo_index) + name + '.xyz')
        xyz_file = open(xyz_path, 'w')
        xyz_file.write(hdo_structures['xyz'][hdo_index])
        xyz_file.close()

    if exp_format == 'host':
        host_path = os.path.join(export_dir, str(hdo_index) + name)
        host_file = open(host_path, 'w')
        host_file.write(hdo_structures['host'][hdo_index])
        host_file.close()


# Required variables:
# - atachment info
def host_import(host_path):
    """
    Convert structure file to host format.
    """
    structure_file = open(structure_path, 'r')
    structure_lines = structure_file.readlines()
    structure_file.close()

    host_file = open(host_path, 'w')
    host_file.write(host_name + '\n')
    host_file.write(' ' + num_atoms + '\t1\n')

    for atom_info in atom_list:
        host_file.write()

    host_file.write(' ' + num_attachment + '\n')

    for attachment_info in attachment_list:
        host_file.write()

    # Add drive info and additional host info

    host_file.close()


def angle_between(C1, C0, C2):
    """ Calculates bond angle for three given atom positions. => C1-CO-C2 """
    import numpy as np
    from math import acos, degrees
    coor1 = []
    coor2 = []
    for c0, c1, c2 in zip(C0, C1, C2):
        coor1.append(c1 - c0)
        coor2.append(c2 - c0)

    angle = np.dot(coor1, coor2) / (np.linalg.norm(coor1) * np.linalg.norm(coor2))

    return degrees(acos(angle))
