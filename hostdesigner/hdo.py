# Reading and writing hostecule files for HostDesigner
# Author: Kutay B Sezginel
# Date: February 2017
import os
import copy
from tabulate import tabulate
from hostdesigner.host import Host
from hostdesigner.visualize import show, write_pdb


hdo_export_dir = os.path.join(os.getcwd(), 'results', 'hdo')


class Hdo:
    """
    HostDesigner output file object.
    """
    def __init__(self, hdo_path, drive=False):
        self.export_dir = hdo_export_dir
        self.path = hdo_path
        if drive:
            self.read_drive()
        else:
            self.read()
            self.read_structures()

    def read(self):
        """
        Read .hdo output files from a HostDesigner run with given path.
        """
        with open(self.path, 'r') as hdo:
            self.lines = hdo.readlines()

        self.n_structures = 0
        self.structures = {'rmsd': [], 'energy': [], 'n_atoms': [], 'atom': [], 'coor': [],
                           'info': [], 'xyz': [], 'host': [], 'linker': [], 'index': []}
        for line_index, line in enumerate(self.lines):
            if line.replace(' ', '').strip('\n').isdigit():
                self.structures['n_atoms'].append(int(line.replace(' ', '').strip('\n')))
            if 'RMSD' in line:
                self.structures['rmsd'].append(line.split(',')[0].split('=')[-1])
                self.structures['energy'].append(float(line.split(',')[2].split('_')[1]))
                self.structures['info'].append(line.strip('_'))
                self.structures['linker'].append(line.strip('_').split('_')[2])
                self.structures['index'].append(self.n_structures)
                self.n_structures += 1

    def tabulate(self, structures=None):
        """
        Print hdo results in table format.
        """
        if type(structures) == int:
            num = structures
        else:
            num = self.n_structures
        print(tabulate({'Structure': self.structures['index'][:num],
                        'Linker': self.structures['linker'][:num],
                        'N_atoms': self.structures['n_atoms'][:num],
                        'RMSD': self.structures['rmsd'][:num],
                        'Energy': self.structures['energy'][:num]}, headers="keys"))

    def sort(self, var='n_atoms', structures=None, table=True):
        """
        Sort hdo results according to given variable and print in table format.
        Available variables:
            - n_atoms
            - energy
            - RMSD (default)
        """
        i = self.structures['index']
        v = self.structures[var]
        sorted_var = sorted(zip(v, i))
        sorted_indices = [i[1] for i in sorted_var]
        sorted_structures = {'rmsd': [], 'energy': [], 'n_atoms': [], 'atom': [], 'coor': [],
                             'info': [], 'xyz': [], 'host': [], 'linker': [], 'index': []}
        for ni, si in enumerate(sorted_indices):
            for k in self.structures:
                sorted_structures[k].append(self.structures[k][si])
        if table:
            if type(structures) == int:
                num = structures
            else:
                num = self.n_structures
            print(tabulate({'Structure': sorted_structures['index'][:num],
                            'Linker': sorted_structures['linker'][:num],
                            'N_atoms': sorted_structures['n_atoms'][:num],
                            'RMSD': sorted_structures['rmsd'][:num],
                            'Energy': sorted_structures['energy'][:num]}, headers="keys"))
        new_hdo = copy.deepcopy(self)
        new_hdo.structures = sorted_structures
        return new_hdo

    def read_structures(self, structures=None):
        """
        Read atom names and coordinates and convert to xyz and host formats.
        """
        if type(structures) == int:
            num = structures
        else:
            num = self.n_structures
        start_line = 2
        for structure in range(min(num, self.n_structures)):
            host_lines = self.structures['info'][structure]
            host_lines += ' ' + str(self.structures['n_atoms'][structure]) + '\t1\n'

            xyz_lines = str(self.structures['n_atoms'][structure]) + '\n'
            xyz_lines += str(structure + 1) + '_' + str(self.structures['rmsd'][structure]) + '\n'

            self.structures['atom'].append([])
            self.structures['coor'].append([])

            end_line = start_line + self.structures['n_atoms'][structure]
            for line_index, line in enumerate(self.lines[start_line:end_line]):
                atom_name = line.split()[0]
                new_line = '  {:2} {:3}{}'.format(atom_name, line_index + 1, line[3:])
                host_lines += new_line

                name, x, y, z = line.split()[0], line.split()[1], line.split()[2], line.split()[3]
                xyz_lines += name + '\t' + x + '\t' + y + '\t' + z + '\n'

                self.structures['atom'][structure].append(name)
                self.structures['coor'][structure].append([float(x), float(y), float(z)])

            self.structures['host'].append(host_lines)
            self.structures['xyz'].append(xyz_lines)
            start_line += self.structures['n_atoms'][structure] + 2

    def read_drive(self):
        """
        Read test drive output
        """
        with open(self.path, 'r') as drive:
            self.lines = drive.readlines()

        self.n_structures = 0
        self.structures = {'rmsd': [], 'energy': [], 'n_atoms': [], 'atom': [], 'coor': [],
                           'info': [], 'xyz': [], 'host': [], 'linker': [], 'index': []}
        xyz_indices = []
        for line_index, line in enumerate(self.lines):
            if 'Drive' in line:
                xyz_indices.append(line_index - 1)

        for i in range(1, len(xyz_indices) - 1):
            start = xyz_indices[i]
            end = xyz_indices[i + 1]
            self.structures['n_atoms'].append(int(self.lines[start].strip()))
            self.structures['xyz'].append(self.lines[start:end])
            self.structures['atom'].append([line.strip().split()[0] for line in self.lines[start + 2:end]])
            self.structures['coor'].append([[float(i) for i in line.strip().split()[1:]] for line in self.lines[start + 2:end]])
            self.structures['info'].append(self.lines[start + 1].strip())
            self.structures['index'].append(i)
            self.structures['host'].append('---')
            self.structures['linker'].append('---')
            self.structures['energy'].append('---')
            self.structures['rmsd'].append('---')

        self.n_structures = len(xyz_indices)

    def show(self, structures=3, start=0, color=False,
             camera='perspective', move='auto', div=5, distance=(-10, 10), axis=0, caps=True, save=None, group=True, rotate=None):
        """
        Show hdo results output structures
        - structures: number of structures to visualize
        - start: start showing structures from given index
        - color: color Dummy atoms separetely
        other arguments are directly used with show method in visualize library.
        """
        hosts = []
        for i in range(structures):
            h = Host()
            h.read_hdo(self.structures, idx=i + start)
            if color:
                hosts.append(h.color())
            else:
                hosts.append(h)
        return show(*hosts, camera=camera, move=move, div=div, distance=distance, axis=axis, caps=caps, save=save, group=group, rotate=rotate)
