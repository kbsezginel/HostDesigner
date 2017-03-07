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


def read_summary(summary_path, text=False):
    """
    Read summary output file.
    """
    with open(summary_path, 'r') as s:
        summary_lines = s.readlines()
    summary = dict(date=None, hosts=[], time=None, links_read=None, links_used=None, structures=None)
    for line in summary_lines:
        if 'Date of run' in line:
            summary['date'] = line.split('Date of run:')[1].strip()
        if 'Name of Host' in line:
            summary['hosts'].append(' '.join(line.split()[4:]).strip())
        if 'Number of links read' in line:
            summary['links_read'] = int(line.split()[-1])
        if 'Number of links used' in line:
            summary['links_used'] = int(line.split()[-1])
        if 'Total number of structures examined' in line:
            summary['structures'] = int(line.split()[-1])
        if 'Total time (sec)' in line:
            summary['time'] = float(line.split()[-1])
    if text:
        return print_summary(summary)
    else:
        return summary


def print_summary(summary):
    """ Format summary output for printing """
    summary_out = """
    Date: %s
    Time (s): %f
    Number of links read: %i
    Number of links used: %i
    Total number of structures examined: %i
    Hosts(s): %s
    """ % (summary['date'], summary['time'], summary['links_read'], summary['links_used'],
           summary['structures'], ('\n%s' % (' ' * 14)).join(summary['hosts']))
    return summary_out
