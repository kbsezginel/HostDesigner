# Visualizing molecule files for HostDesigner
# Author: Kutay B Sezginel
# Date: February 2017
import os
import math
import tempfile
import nglview


def show(*args, camera='perspective', move='auto', div=5, distance=(-10, -10), axis=0, caps=True, save=None, group=True, rotate=None):
    """
    Show given structure(s) using nglview
        - save: file_path -> save pdb file to given path
        - camera: 'perspective' / 'orthographic'
        - move: 'auto' / 'single' / None -> separate multiple structures equally distant from each other
        - distance: tuple (x, y) -> separation distance
        - axis: 0 / 1 / 2 -> separation direction (0:x, 1:y, 2:z)
        - caps: bool -> capitalize atom names so they show true colors in nglview+
        - group: bool -> group molecules in pdb file (to show bonds properly in nglview)
        - rotate: 0 / 1 / 2 -> rotate 90 degrees in (0:x, 1:y, 2:z) axis (view in xy plane)
    """
    if move is 'auto':
        translation_vectors = arrange_structure_positions(len(args), div=div, distance=distance)
    elif move is 'single':
        translation_vectors = axis_translation(len(args), distance=distance[0], axis=axis)
    else:
        translation_vectors = [[0, 0, 0]] * len(args)

    atom_names = []
    atom_coors = []
    group_numbers = []
    for mol_index, (molecule, vec) in enumerate(zip(args, translation_vectors), start=1):
        atom_names += molecule.atom_names
        coordinates = molecule.atom_coors
        if rotate is not None:
            coordinates = rotate_90(coordinates, axis=rotate)
        coordinates = translate(coordinates, vector=vec)
        atom_coors += coordinates
        group_numbers += [mol_index] * len(molecule.atom_names)

    # nglview require atom names in all caps to color them properly
    if caps:
        atom_names = [name.upper() for name in atom_names]

    # nglview requires molecules to be grouped separately to show proper bonding
    if not group:
        group_numbers = [1] * len(atom_names)
    temp_pdb_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb')
    write_pdb(temp_pdb_file, atom_names, atom_coors, group=group_numbers)

    view = nglview.show_structure_file(temp_pdb_file.name)
    view.camera = camera
    temp_pdb_file.close()

    if save is not None:
        with open(save, 'w') as save_file:
            write_pdb(save_file, atom_names, atom_coors, group=group_numbers)

    return view


def arrange_structure_positions(n_structures, div=5, distance=(10, 10)):
    """
    Arrange structure positions according to number of structures given.
    """
    n_structures_lateral = div
    split = math.ceil(n_structures / n_structures_lateral)
    vertical_positions = axis_translation(split, distance=distance[1], axis=1)
    translation_vectors = []
    for s in range(split):
        if s < 2:
            n_split = math.ceil(n_structures / split)
        else:
            n_split = math.floor(n_structures / split)
        new_vectors = axis_translation(n_split, distance=distance[0], axis=0)
        new_vectors = translate(new_vectors, vector=vertical_positions[s])
        translation_vectors += new_vectors
    return translation_vectors


def translate(atom_coors, vector=[-10, 0, 0]):
    """ Translate given coordinates with given vector """
    translated_coors = []
    x, y, z = vector
    for coor in atom_coors:
        new_coor = [coor[0] + x, coor[1] + y, coor[2] + z]
        translated_coors.append(new_coor)
    return translated_coors


def axis_translation(n_structures, distance=-10, axis=0):
    """
    Automatically adjust structure positions equally distant from each other in given axis
        - distance: distance between each structure
        - axis: axis selection for translation (0: x-axis, 1: y-axis, 2: z-axis)
    """
    translation_vectors = []
    lim = (n_structures - 1) * distance / 2
    for i in range(n_structures):
        vec = [0, 0, 0]
        vec[axis] = -lim + i * distance
        translation_vectors.append(vec)
    return translation_vectors


def rotate_90(atom_coors, axis=0):
    """
    Rotate coordinates 90 degrees in x, y or z axis.
    """
    if axis == 0:             # Rotate 90 degrees on x-axis
        new_coors = [[c[0], c[2], c[1]] for c in atom_coors]
    elif axis == 1:           # Rotate 90 degrees on y-axis
        new_coors = [[c[2], c[1], c[0]] for c in atom_coors]
    elif axis == 2:           # Rotate 90 degrees on z-axis
        new_coors = [[c[1], c[0], c[2]] for c in atom_coors]

    return new_coors


def write_pdb(pdb_file, names, coors, group=None, header='Host'):
    """ Write given atomic coordinates to file object in pdb format """
    pdb_file.write('HEADER    ' + header + '\n')
    format = 'HETATM%5d%3s  M%4i %3i     %8.3f%8.3f%8.3f  1.00  0.00          %2s\n'
    if group is None:
        group = [1] * len(names)
    for atom_index, (atom_name, atom_coor) in enumerate(zip(names, coors), start=1):
        x, y, z = atom_coor
        residue_no = group[atom_index - 1]
        pdb_file.write(format % (atom_index, atom_name, residue_no, residue_no, x, y, z, atom_name.rjust(2)))
    pdb_file.write('END\n')
    pdb_file.flush()
