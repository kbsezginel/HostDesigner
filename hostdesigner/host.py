# Reading and writing hostecule files for HostDesigner
# Author: Kutay B Sezginel
# Date: February 2017
import os
import yaml
from hostdesigner.visualize import write_pdb


class Host:
    """
    Host class.
    """
    def __init__(self, host_input=None, run_type=None):
        self.run_type = run_type
        self.export_dir = os.getcwd()
        if type(host_input) is str:
            self.read(host_input)
        elif type(host_input) is dict:
            self.read_dict(host_input)
        else:
            self.name = ''
            self.atom_names = []
            self.atom_coors = []

    def __repr__(self):
        return "<Host object: %s>" % self.name

    def read(self, host_path):
        """ Read host file """
        with open(host_path, 'r') as host:
            self.lines = host.readlines()

        self.name = self.lines[0].strip()
        self.n_atoms, self.n_guest = [int(i) for i in self.lines[1].split()]

        self.read_structure()
        if self.run_type is not None:
            self.read_attachments()
            self.read_drive()

    def read_structure(self):
        """ Read structure information """
        self.atom_names = []
        self.atom_coors = []
        self.atom_types = []
        self.bonds = []
        for line in self.lines[2:self.n_atoms + 2]:
            name, idx, x, y, z, a_type = line.split()[:6]
            bonds = line.split()[6:]

            self.atom_names.append(name)
            self.atom_coors.append([float(i) for i in [x, y, z]])
            self.atom_types.append(int(a_type))
            self.bonds.append(bonds)

    def read_attachments(self):
        """ Read attachment information """
        self.n_attachments = None
        self.attachments = []
        if self.run_type == 'LINKER':
            self.n_attachments = int(self.lines[self.n_atoms + 2].strip())
            for line in self.lines[self.n_atoms + 3:self.n_atoms + self.n_attachments + 3]:
                self.attachments.append(line.split())
            self.attachment_list = [int(i[0]) for i in self.attachments]
        elif self.run_type == 'OVERLAY':
            attachment_line = self.lines[self.n_atoms + 2].strip().split()
            self.n_attachments = int(attachment_line[0])
            if len(attachment_line) > 1:
                self.n_sym = int(attachment_line[1])
                self.n_attachments = self.n_attachments * (self.n_sym + 1)
            else:
                self.n_sym = 0
            for line in self.lines[self.n_atoms + 3:self.n_atoms + self.n_attachments + 3]:
                self.attachments.append(line.split())
            self.attachment_list = [int(i[0]) for i in self.attachments]
        else:
            print('No run type found for %s!' % self.name)

    def read_drive(self):
        """ Read drive information if exists """
        self.drive = []
        if len(self.lines) > (self.n_atoms + self.n_attachments + 3):
            drive_start = self.n_atoms + self.n_attachments + 4
            self.n_drive = int(self.lines[drive_start - 1].split()[1])
            drive_end = drive_start + self.n_drive
            for line in self.lines[drive_start:drive_end]:
                self.drive.append(line.split())
        else:
            print('No drive information found for %s!' % self.name)

    def read_hdo(self, hdo_structures, idx=0):
        """ Generate host object from hdo results file. """
        self.atom_coors = hdo_structures['coor'][idx]
        self.atom_names = hdo_structures['atom'][idx]
        self.name = hdo_structures['linker'][idx]

    def remove_atoms(self, atom=['X']):
        """ Remove all atoms with given name from host object """
        new_coors = []
        new_names = []
        for name, coor in zip(self.atom_names, self.atom_coors):
            for rmv_atom in atom:
                if rmv_atom != name:
                    new_names.append(name)
                    new_coors.append(coor)
        self.atom_names = new_names
        self.atom_coors = new_coors

    def export(self, file_format='pdb', export_dir=None, name=None, atom_type=None):
        """ Export host object """
        if export_dir is None:
            export_dir = self.export_dir
        if name is None:
            name = self.name
        if file_format is 'pdb':
            pdb_path = os.path.join(export_dir, '%s.pdb' % name)
            with open(pdb_path, 'w') as pdb_file:
                write_pdb(pdb_file, self.atom_names, self.atom_coors, header=name)
        elif file_format is 'dict':
            self.to_dict()
        elif file_format is 'host':
            host_path = os.path.join(export_dir, '%s' % name)
            self.write_host(host_path, atom_type=atom_type)
        elif file_format is 'yaml':
            yaml_path = os.path.join(export_dir, '%s.yaml' % name)
            yaml.dump(self, open(yaml_path, 'w'))

    def write_host(self, host_path, atom_type=None):
        """ Export host object to HostDesigner input structures (host) format """
        with open(host_path, 'w') as h:
            h.write('%s\n' % self.name)
            h.write(' %i %i\n' % (self.n_atoms, self.n_guest))
            i = 0
            for atom, coor in zip(self.atom_names, self.atom_coors):
                x, y, z = coor
                bonds = '    '.join(self.bonds[i])
                if atom_type is None:
                    # Control file must include 'notype' keyword to read this format properly
                    h.write('  %s    %i    %f    %f   %f    %s\n' % (atom, i + 1, x, y, z, bonds))
                else:
                    h.write('  %s    %i    %f    %f   %f    %s    %s\n' % (atom, i + 1, x, y, z, atom_type, bonds))
                i += 1
            # Attachments:
            h.write(' %i\n' % self.n_attachments)
            for a in self.attachments:
                h.write('  %s\n' % '  '.join(a))

    def to_dict(self):
        """ Export host information to dictionary """
        return dict(name=self.name, n_atoms=self.n_atoms, n_guest=self.n_guest,
                    atom_names=self.atom_names, atom_coors=self.atom_names, bonds=self.bonds,
                    n_attachments=self.n_attachments, attachments=self.attachments,
                    n_drive=self.n_drive, drive=self.drive)
