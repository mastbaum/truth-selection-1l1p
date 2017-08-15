from glob import glob
import sys
from ROOT import gallery, galleryfmwk, TFile

def process_files(files):
    # Create ana_processor
    my_proc = galleryfmwk.ana_processor()

    # Set input root file
    print glob(files)
    for _f in glob(files):
        print 'FILE', _f
        ff = TFile(_f)
        ok = ff.IsOpen() and not ff.IsZombie()
        ff.Close()

        if ok:
            my_proc.add_input_file(_f)
    
    # Set output ROOT file name
    my_proc.set_ana_output_file("ana_out.root")

    # Configure analysis module
    exampleModule = galleryfmwk.make_dedx_pdfs()
    #exampleModule.setEventWeightProducer("eventweight")
    exampleModule.setMCTruthProducer("generator");
    #exampleModule.setMCFluxProducer("generator");
    exampleModule.setMCTrackProducer("mcreco");
    exampleModule.setMCShowerProducer("mcreco"); 
    exampleModule.setVerbose(True)
    
    # Attach analysis module
    my_proc.add_process(exampleModule)

    my_proc.run()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage:', sys.argv[0], 'input.root'
        sys.exit(1)

    process_files(sys.argv[1])

