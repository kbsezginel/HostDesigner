import os
import sys
import subprocess


def hd_run(control):
    hd = subprocess.run(['hd3.0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    hd_out = hd.stdout.decode().split('\n')


class Control:
    def __init__(self, run_type):
        if run_type == 'LINKER':
            self.run = 'LINK'
        elif run_type = 'OVERLAY':
            self.run = 'OVER'
        else:
            print('Please enter run type [LINKER or OVERLAY]')

    def export():
        self.path = os.path.join(run_path, 'control')
        with open(self.path, 'w') as f:
            f.write()
