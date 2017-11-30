import glob
import sys

file_list = glob.glob(sys.argv[1])
f_out = sys.argv[2]
with open(f_out, "w+") as f:
    for f_name in file_list:
        f.write("%s\n" % f_name)

