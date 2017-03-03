# Reading and writing hostecule files for HostDesigner
# Author: Kutay B Sezginel
# Date: February 2017
from tabulate import tabulate
import os
from hostdesigner.host import Host
from hostdesigner.visualize import show, write_pdb


hdo_export_dir = os.path.join(os.getcwd(), 'results', 'hdo')


class Hdo:
    """
    HostDesigner output file object.
    """
    def __init__(self, path=None, run_type=None):
        self.run_type = run_type
        self.export_dir = hdo_export_dir
        if type(path) is str:
            self.path = path
            self.read(path)
            self.read_structures()

    def read(self, hdo_path=None):
        """
        Read .hdo output files from a HostDesigner run with given path.
        """
        hdo_path = self.path if type(hdo_path) is not str else hdo_path
        with open(hdo_path, 'r') as hdo:
            self.lines = hdo.readlines()

        self.structures = []
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
        new_hdo = Hdo()
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

    def show(self, structures=3, start=0, move='auto', div=5, distance=(-10, 10), axis=0, rename='F'):
        """
        Show hdo results output structures
        - structures: number of structures to visualize
        """
        hosts = []
        for i in range(structures):
            h = Host()
            h.read_hdo(self.structures, idx=i + start)
            hosts.append(h)
        return show(*hosts, move=move, div=div, distance=distance, axis=axis, rename=rename)
