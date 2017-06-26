# Control file for HostDesigner
# Author: Kutay B Sezginel
# Date: February 2017
import os

""" Control file keywords """
sample = dict(
    run_type=None,      # Type of HostDesigner run (LINK or OVER)
    hosta=None,         # File name for OVERLAY input or LINKER input fragment 1 (default: hosta`)
    hostb=None,         # File name for LINKER input fragment 2 (default: hostb)
    notype=False,       # If true eliminates the need to include atom type in input
    out=None,           # Prefix for the output files (default: out) -> out.summ, out_1.hdo, out_2.hdo
    numview=None,       # Number of structures to print (def:100, max:100)
    xyz=False,          # Limits file output to atom name and coordinates only
    linklib=None,       # Path for alternative library (ex: linklib=/.../LIBRARY_syn)
    mirroraoff=False,   # If true turns off using mirror image of hosta for building structures
    mirrorboff=False,   # If true turns off using mirror image of hostb for building structures
    drivea=False,       # Apply geometry drive for hosta
    driveb=False,       # Apply geometry drive for hostb
    testdrive=False,    # Run test for geometry drives only and export to: out_testa.hdo and out_testb.hdo
    metshape=None,      # Check binding site topology (options: TETR, SQPL, TBPY, SQPY, OCTA | default: NONE)
    useclass=False,     # If true only first structure from each class in the link library is used
    nochiral=False,     # If true all chiral links are discarded
    noprochiral=False,  # If true all prochiral links are discarded
    noasym=False,       # All links that would give a different connectivity when the attachment points are switched are discarded
    minconn=None,       # Min. number of atoms for linker btw. attachment points (default: 0)
    maxconn=None,       # Max. number of atoms for linker btw. attachment points (default: 20)
    maxconfe=None,      # Max. conformational energy for linker (default: 100.0 kcal/mol)
    maxnrot=None,       # Max. num. of rotatable bonds for linker + input (default: 20)
    maxrmsd=None,       # Max. RMSD for superposition (default: 100.0 Angstrom)
    tightness=None)     # Max. torsion angle deviation in OVERLAY (default: 1.0 corresponds to 20° -> 20°/tightness)


def clean(control):
    """ Remove unused keywords in control dictionary """
    cont = dict(control)
    for kwd in control:
        if cont[kwd] is False or cont[kwd] is None:
            del cont[kwd]
    return cont


def generate(control, export_dir, space=20):
    """ Generate control file with specified keywords """
    cont = clean(control)
    if len(cont) == 1:
        cont_lines = ['%s\n' % cont['run_type']]
    else:
        cont_lines = ['%s%sAND\n' % (cont['run_type'], (space - len(cont['run_type'])) * ' ')]
        del cont['run_type']
        for i, kwd in enumerate(cont):
            if type(cont[kwd]) is bool:
                cl = len(kwd)
                end = '%sAND\n' % (' ' * (space - cl)) if i < len(cont) - 1 else '\n'
                cont_lines.append('%s%s' % (kwd, end))
            else:
                cl = len(kwd) + len(str(cont[kwd]))
                end = '%sAND\n' % (' ' * (space - cl - 1)) if i < len(cont) - 1 else '\n'
                cont_lines.append('%s=%s%s' % (kwd, str(cont[kwd]), end))
    control_path = os.path.join(export_dir, 'control')
    with open(control_path, 'w') as c:
        for line in cont_lines:
            c.write(line)
    return cont


def overlay(hosta=None, export_dir=None):
    """ Get a generic OVERLAY control file """
    cont = dict()
    if hosta is not None:
        cont['hosta'] = '%s  ' % hosta
    cont['run_type'] = 'OVER'
    cont['notype'] = True
    cont['maxconn'] = 50
    if export_dir is not None:
        generate(cont, export_dir)
    else:
        return cont
