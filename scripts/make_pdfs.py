'''Make dE/dx PDFs for the truth-based selection.

Usage:

  $ python make_pdfs.py "files*.root"

The input is a set of art ROOT files.

A. Mastbaum <mastbaum@uchicago.edu>, 2017/09
'''

import sys
from glob import glob
from ROOT import galleryfmwk

def process_files(files):
    # Create ana_processor
    my_proc = galleryfmwk.ana_processor()

    # Set input root file
    for _f in glob(files):
        # May check here if files open correctly
        my_proc.add_input_file(_f)
    
    # Set output ROOT file name
    my_proc.set_ana_output_file("ana_out.root")

    # Configure analysis module
    exampleModule = galleryfmwk.TSPDFGen()
    exampleModule.setMCTruthProducer("generator");
    exampleModule.setMCTrackProducer("mcreco");
    exampleModule.setMCShowerProducer("mcreco"); 
    exampleModule.setVerbose(True)
    
    # Attach analysis module
    my_proc.add_process(exampleModule)

    my_proc.run()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage:', sys.argv[0], '"input*.root"'
        sys.exit(1)

    process_files(sys.argv[1])

