import sys
import os

import ROOT

def loop_merge(source_dir, dir_name, temp_dir, min_selection, max_selection, append, final_name, max_dir_name):
    nselections = max_selection - min_selection+1
    n_events = 0
    data_tree_list_per_selection = [[] for i in range(nselections)]
    file_list_per_selection = [[] for i in range(nselections)]
    t_list_per_selection = [ROOT.TList() for i in range(nselections)]
    merge_count = 0

    n_good_1e1p = [0 for i in range(nselections)]
    n_true_1e1p = [0 for i in range(nselections)]
    n_miss_1e1p = [0 for i in range(nselections)]

    n_good_1m1p = [0 for i in range(nselections)]
    n_true_1m1p = [0 for i in range(nselections)]
    n_miss_1m1p = [0 for i in range(nselections)]

    while 1:
        if dir_name > max_dir_name:
            break
        for i in range(nselections):
            fname = source_dir + str(dir_name) + "/out%i.root" % (i + min_selection)
            f = open_root(fname)
            file_list_per_selection[i].append(f)
            if f is None:
                print "MISSING or CORRUPTED FILE %i in DIRECTORY %i" % (i + min_selection, dir_name)
                continue
            """
            if f is None:
                fname_2 = fname.replace("_try2","")
                f2 = open_root(fname_2)
                if f2 is None:
                    print "MISSING or CORRUPTED FILE %i in DIRECTORY %i" % (i + min_selection, dir_name)
                    continue
                else:
                    file_list_per_selection[i][-1] = f2
            """
            for j, entry in enumerate(f.Get("header")):
		n_good_1e1p[i] += entry.good_1e1p
		n_true_1e1p[i] += entry.true_1e1p
		n_miss_1e1p[i] += entry.miss_1e1p
		
		n_good_1m1p[i] += entry.good_1m1p
		n_true_1m1p[i] += entry.true_1m1p
		n_miss_1m1p[i] += entry.miss_1m1p
                assert(j==0)

            data_tree_list_per_selection[i].append( file_list_per_selection[i][-1].Get("data") )
            t_list_per_selection[i].Add(data_tree_list_per_selection[i][-1])

        merge_count += 1
        max_n_merge = 100 / nselections
        if (dir_name + 1) % max((max_n_merge / 10),1) == 0:
            print "Processed %i files" % ((dir_name + 1) * nselections)
            print "AT DIRECTORY: %i" % dir_name
        if merge_count % max(max_n_merge,1) == 0:
            break
        dir_name += 1

    if append:
        print "LOAD PREVIOUS AND MERGE AND WRITE" 
        for i in range(nselections):
	    temp_file_a = temp_dir + "temp%i" % (i+min_selection) + final_name
	    last_file = ROOT.TFile(temp_file_a)
	    data_tree = last_file.Get("data")
	    t_list_per_selection[i].Add(data_tree)
		
            header_file = file_list_per_selection[i][-1]
		
            n_good_1e1p[i] += last_file.Get("good_1e1p").GetUniqueID()
            n_true_1e1p[i] += last_file.Get("true_1e1p").GetUniqueID()
            n_miss_1e1p[i] += last_file.Get("miss_1e1p").GetUniqueID()
            n_good_1m1p[i] += last_file.Get("good_1m1p").GetUniqueID()
            n_true_1m1p[i] += last_file.Get("true_1m1p").GetUniqueID()
            n_miss_1e1p[i] += last_file.Get("miss_1m1p").GetUniqueID()
      

	    temp_file_b = temp_dir + "temp%iB" % (i+min_selection) + final_name
	    merge_and_write(t_list_per_selection[i], header_file, temp_file_b, 
                             n_good_1e1p[i], n_true_1e1p[i], n_miss_1e1p[i], 
                             n_good_1m1p[i], n_true_1m1p[i], n_miss_1m1p[i])
		
            last_file.Close()
    else:
        print "MERGE AND WRITE"
        for i in range(nselections):
            header_file = file_list_per_selection[i][-1]
	    temp_file_b = temp_dir + "temp%iB" % (i+min_selection) + final_name
	    merge_and_write(t_list_per_selection[i], header_file, temp_file_b, 
                             n_good_1e1p[i], n_true_1e1p[i], n_miss_1e1p[i], 
                             n_good_1m1p[i], n_true_1m1p[i], n_miss_1m1p[i])

    print "CLOSING FILES"
    [[f.Close() for f in lst if f is not None] for lst in file_list_per_selection]
    for i in range(nselections):
        temp_file_a = temp_dir + "temp%i" % (i+min_selection) + final_name
        temp_file_b = temp_dir + "temp%iB" % (i+min_selection) + final_name
        print temp_file_a
        print temp_file_b
        os.rename(temp_file_b, temp_file_a)

    return dir_name + 1

def open_root(fname):
    if not os.path.isfile(fname):
        return None
    
    f = ROOT.TFile(fname)
    # check if file is corrupted
    try:
        f.Get("truth").GetEntries()
    except AttributeError:
        return None
    return f

def merge_and_write(t_list, header_file, f_out_name, n_good_1e1p, n_true_1e1p, n_miss_1e1p, n_good_1m1p, n_true_1m1p, n_miss_1m1p):
    f_out = ROOT.TFile(f_out_name, "RECREATE")
    f_out.cd()
    merged = ROOT.TTree.MergeTrees(t_list)
    merged.SetName("data")
    header = header_file.Get("header")

    good_1e1p = ROOT.TNamed("good_1e1p", "good_1e1p")
    good_1e1p.SetUniqueID(n_good_1e1p)
    true_1e1p = ROOT.TNamed("true_1e1p", "true_1e1p")
    true_1e1p.SetUniqueID(n_true_1e1p)
    miss_1e1p = ROOT.TNamed("miss_1e1p", "miss_1e1p")
    miss_1e1p.SetUniqueID(n_miss_1e1p)

    good_1m1p = ROOT.TNamed("good_1m1p", "good_1m1p")
    good_1m1p.SetUniqueID(n_good_1m1p)
    true_1m1p = ROOT.TNamed("true_1m1p", "true_1m1p")
    true_1m1p.SetUniqueID(n_true_1m1p)
    miss_1m1p = ROOT.TNamed("miss_1m1p", "miss_1m1p")
    miss_1m1p.SetUniqueID(n_miss_1m1p)

    good_1e1p.Write()
    miss_1e1p.Write()
    true_1e1p.Write()
    good_1m1p.Write()
    miss_1m1p.Write()
    true_1m1p.Write()

    merged.Write()
    header.CloneTree().Write()
    f_out.Close()

def main(source_dir, output_dir, identifier, temp_dir, min_selection, max_selection, min_dir, max_dir):
    append = min_dir != 0
    dir_name = min_dir
    max_dir_name = max_dir
    while 1:
        dir_name = loop_merge(source_dir, dir_name, temp_dir, min_selection, max_selection, append, identifier, max_dir_name)
        append = True
        if dir_name > max_dir_name:
            break
    for i in range(min_selection,max_selection+1):
        os.rename(temp_dir + "temp%i" % i + identifier, output_dir + identifier + "%i.root" % i)
    print "FINISHED"

if __name__=="__main__":
    source_dir = sys.argv[1]
    output_dir = sys.argv[2]
    output_identifier = sys.argv[3]
    temp_dir = sys.argv[4]
    min_selection = int(sys.argv[5])
    max_selection = int(sys.argv[6])
    min_dir = int(sys.argv[7])
    max_dir = int(sys.argv[8])
   
    main(source_dir, output_dir, output_identifier, temp_dir, min_selection, max_selection, min_dir, max_dir) 

