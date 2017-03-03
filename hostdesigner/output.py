# Read output files of HostDesigner
# Author: Kutay B Sezginel
# Date: February 2017
import os
from tabulate import tabulate
from hostdesigner.hdo import Hdo


def find_results(results_dir):
    """
    Find HostDesigner results in a given directory
    - out_1.hdo: sorted according to RMSD
    - out_2.hdo: sorted according to energy
    """
    results = dict(hdo=[], summ=[])
    for f in os.listdir(results_dir):
        if os.path.splitext(f)[1] == '.hdo':
            results['hdo'].append(os.path.join(results_dir, f))
        elif os.path.splitext(f)[1] == '.summ':
            results['summ'].append(os.path.join(results_dir, f))
    return results


def read_hdo(hdo_path):
    """
    Read .hdo output files from a HostDesigner run with given path.
    """
    return Hdo(hdo_path)
