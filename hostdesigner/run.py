import os
import subprocess


def hd_run(run_dir):
    """ Make a HostDesigner run in given directory """
    hd = subprocess.run(['hd3.0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=run_dir)
    hd_out = hd.stdout.decode().split('\n')
    return hd_out
