'''Create and print out a list of POT for each file in the input.

You can then trim this and add the lines together with something
like np.sum(np.loadtxt(file.txt)). This is a workaround since
somehow just getting the totgoodpot branch inside ROOT doesn't
work at all.

Usage:

  $ python pot.py "files*.root"

The input is a set of art ROOT files.

A. Mastbaum <mastbaum@uchicago.edu>, 2017/09/07
'''

import sys
import ROOT

import subprocess32
import tempfile
from decimal import *

def main(glob):
    tfil = tempfile.TemporaryFile()
    proc = subprocess32.Popen(['python', 'pot.py', glob], stdout=tfil)
    proc.wait()

    tfil.seek(0)
    s = Decimal(0)
    for line in tfil:
        info = line.split("*")
        try:
            s += Decimal(info[2].rstrip().lstrip())
        except: 
            pass
    print "_____FINAL SUM______: %f" % s

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage:', sys.argv[0], '"files*.root"'
        sys.exit(1)

    main( sys.argv[1] )

