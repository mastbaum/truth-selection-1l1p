import sys
import os

import ROOT

def loop_merge(source_dir, dir_name, temp_dir, n_selections, append, final_name, max_dir_name):
    n_events = 0
    data_tree_list_per_selection = [[] for i in range(n_selections)]
    file_list_per_selection = [[] for i in range(n_selections)]
    t_list_per_selection = [ROOT.TList() for i in range(n_selections)]
    merge_count = 0

    while 1:
        if dir_name > max_dir_name:
            break
        for i in range(n_selections):
            fname = source_dir + str(dir_name) + "/out%i.root" % i
            f = open_root(fname)
            file_list_per_selection[i].append(f)
            if f is None:
                fname_2 = fname.replace("_try2","")
                f2 = open_root(fname_2)
                if f2 is None:
                    print "MISSING or CORRUPTED FILE %i in DIRECTORY %i" % (i, dir_name)
                    continue
                else:
                    file_list_per_selection[i][-1] = f2

            data_tree_list_per_selection[i].append( file_list_per_selection[i][-1].Get("data") )
            t_list_per_selection[i].Add(data_tree_list_per_selection[i][-1])

        merge_count += 1
        max_n_merge = 10 / n_selections
        if (dir_name + 1) % max((max_n_merge / 10),1) == 0:
            print "Processed %i files" % ((dir_name + 1) * n_selections)
            print "AT DIRECTORY: %i" % dir_name
        if merge_count % max(max_n_merge,1) == 0:
            break
        dir_name += 1

    if append:
        print "LOAD PREVIOUS AND MERGE AND WRITE" 
        for i in range(n_selections):
            if file_list_per_selection[i][-1] is not None:
		temp_file_a = temp_dir + "temp%i" % i + final_name
		last_file = ROOT.TFile(temp_file_a)
		data_tree = last_file.Get("data")
		t_list_per_selection[i].Add(data_tree)
		
                header_file = file_list_per_selection[i][-1]
		
		temp_file_b = temp_dir + "temp%iB" % i + final_name
		merge_and_write(t_list_per_selection[i], header_file, temp_file_b)
		
		last_file.Close()
    else:
        print "MERGE AND WRITE"
        for i in range(n_selections):
            if file_list_per_selection[i][-1] is not None:
                header_file = file_list_per_selection[i][-1]
		temp_file_b = temp_dir + "temp%iB" % i + final_name
		merge_and_write(t_list_per_selection[i], header_file, temp_file_b)

    print "CLOSING FILES"
    [[f.Close() for f in lst if f is not None] for lst in file_list_per_selection]
    for i in range(n_selections):
        temp_file_a = temp_dir + "temp%i" % i + final_name
        temp_file_b = temp_dir + "temp%iB" % i + final_name
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

def merge_and_write(t_list, header_file, f_out_name):
    f_out = ROOT.TFile(f_out_name, "RECREATE")
    f_out.cd()
    merged = ROOT.TTree.MergeTrees(t_list)
    merged.SetName("data")
    header = header_file.Get("header")

    merged.Write()
    header.CloneTree().Write()
    f_out.Close()

def main(source_dir, output_dir, identifier, temp_dir, n_selections, min_dir, max_dir):
    append = min_dir != 0
    dir_name = min_dir
    max_dir_name = max_dir
    while 1:
        dir_name = loop_merge(source_dir, dir_name, temp_dir, n_selections, append, identifier, max_dir_name)
        append = True
        if dir_name > max_dir_name:
            break
    for i in range(n_selections):
        os.rename(temp_dir + "temp%i" % i + identifier, output_dir + identifier + "%i.root" % i)

if __name__=="__main__":
    source_dir = sys.argv[1]
    output_dir = sys.argv[2]
    output_identifier = sys.argv[3]
    temp_dir = sys.argv[4]
    n_selections = int(sys.argv[5])
    min_dir = int(sys.argv[6])
    max_dir = int(sys.argv[7])
   
    main(source_dir, output_dir, output_identifier, temp_dir, n_selections, min_dir, max_dir) 

