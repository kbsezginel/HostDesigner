import os
import shutil
import subprocess
from datetime import datetime
from hostdesigner.output import read_summary
from hostdesigner.control import overlay, generate


def hd_run(run_dir, control=None, verbose=1, archive=False):
    """
    Make a HostDesigner run in given directory
    - run_dir: Directory of the simulation
    - control: dict -> If a control object is given it will be saved to run directory
    - verbose: 1 / 2 -> Verbosity
    - archive: bool -> Create archive of run environment
    """
    summary_file = 'out'
    if control is not None:
        run_control = generate(control, run_dir)
        if 'out' in run_control:
            summary_file = run_control['out']
        summary_path = os.path.join(run_dir, 'out.summ')
        print('Control file generated')
    print('Starting HostDesigner run in %s' % run_dir) if verbose <= 2 else None
    hd = subprocess.run(['hd3.0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=run_dir)
    stdout = hd.stdout.decode()
    stderr = hd.stderr.decode()
    if verbose == 2:
        summary_path = os.path.join(run_dir, '%s.summ' % summary_file)
        summary = read_summary(summary_path, text=True)
        print("Stdout:\n\n%s\nStderr:\n%s\n\nSummary:\n%s\n" % (stdout, stderr, summary))
    if archive:
        record_name = record_environment(run_dir)
        print('Environment archived to: %s' % record_name)


def initialize(run_dir, hosta, hostb=None, run_type=None):
    """ initialize host designer run into given directory """
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    hosta.export(file_format='host', export_dir=run_dir)
    if hostb is not None:
        hosta_path = os.path.join(run_dir, '%s' % hosta.name)
        hostb_path = os.path.join(run_dir, '%s' % hostb.name)
        shutil.copy(hosta_path, hostb_path)
    if run_type is 'OVERLAY':
        overlay(hosta=hosta.name, export_dir=run_dir)


def record_environment(run_dir, record_name=None, archive_format='tar'):
    """
    Record HostDesigner run environment
    - record_name: output file name (if not given run directory and timestamp is used)
    - archive: archive format
    """
    if record_name is None:
        record_name = '%s_%f' % (os.path.basename(run_dir), datetime.now().timestamp())
        record_name = os.path.join(run_dir, record_name)
    shutil.make_archive(record_name, archive_format, run_dir)
    return '%s.%s' % (record_name, archive_format)
