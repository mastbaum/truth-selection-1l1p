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

def main(files):
    t = ROOT.TChain('tree')
    for fi in files:
        t.Add('%s/SubRuns' % fi.rstrip("\n"))

    t.SetScanField(0)
    t.Scan("sumdata::POTSummary_generator__GenieGen.obj.totgoodpot")
    

if __name__ == '__main__':
    from glob import glob

    if len(sys.argv) < 2:
        print 'Usage:', sys.argv[0], '"files*.root"'
        sys.exit(1)

    with open(sys.argv[1]) as f:
        main(f.readlines())

