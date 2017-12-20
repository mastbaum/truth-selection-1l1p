import glob
import sys
import ROOT

file_list = glob.glob(sys.argv[1])
f_out = sys.argv[2]
with open(f_out, "w+") as f:
    for f_name in file_list:
        if not ROOT.TFile(f_name).IsZombie():
            f.write("%s\n" % f_name)
        else:
            print "ZOMBIE: %s" % f_name

