import os
import subprocess
from hostdesigner.output import read_summary


def hd_run(run_dir, verbose=False):
    """ Make a HostDesigner run in given directory """
    hd = subprocess.run(['hd3.0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=run_dir)
    stdout = hd.stdout.decode()
    stderr = hd.stderr.decode()
    if verbose:
        summary_path = os.path.join(run_dir, 'out.summ')
        summary = read_summary(summary_path, text=True)
        print("Stdout:\n\n%s\nStderr:\n%s\n\nSummary:\n%s\n" % (stdout, stderr, summary))
