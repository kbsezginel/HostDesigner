import os


def host_export(host_path, export_dir):
    """
    Convert host structure file to xyz format.
    """
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
