import sys
from glob import glob
from ROOT import galleryfmwk

def process_files(outfile, dataset_id, files):
    ds_id = int(dataset_id)

    # Create ana_processor
    my_proc = galleryfmwk.ana_processor()

    # Set input root file
    for _f in glob(files):
        my_proc.add_input_file(_f)
    
    # Set output ROOT file name
    my_proc.set_ana_output_file(outfile)

    # Configure analysis module
    exampleModule = galleryfmwk.TSSelection()
    exampleModule.setFluxWeightProducer("eventweight")
    exampleModule.setEventWeightProducer("mcweight")
    exampleModule.setMCTruthProducer("generator");
    exampleModule.setMCTrackProducer("mcreco");
    exampleModule.setMCShowerProducer("mcreco"); 
    exampleModule.setVerbose(True)
    exampleModule.setDatasetID(ds_id)
    
    # Attach analysis module
    my_proc.add_process(exampleModule)

    my_proc.run()


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Usage:', sys.argv[0], 'output.root dataset_id "input*.root"'
        sys.exit(1)

    process_files(*sys.argv[1:])

