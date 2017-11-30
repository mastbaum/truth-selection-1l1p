import sys
import os

import ROOT

def loop_merge(dir_name, data_dir, temp_dir, append, n_data, final_name):
    last_failed = False
    n_events = 0
    data_tree_list = []
    file_list = []
    t_list = ROOT.TList()
    merge_count = 0

    while 1:
        if dir_name == 926:
            print "END OF DATA"
            last_failed = True
            break
        fname = datadir + str(dir_name) + "/out.root"
        if not os.path.isfile(fname):
            if last_failed:
                print "END OF DATA"
                break
            else:
                last_failed = True
                dir_name += 1
                continue
        else:
            if last_failed:
                print "MISSING DIRECTIORY %i " % (dir_name-1)
            last_failed = False

        file_list.append(ROOT.TFile(fname))
        data_tree_list.append( file_list[-1].Get("data") )
        try:
            n_events += file_list[-1].Get("truth").GetEntries()
        except:
            print "CORRUPTED FILE: %s" % fname
            dir_name += 1
            continue
        t_list.Add(data_tree_list[-1])
        n_data += data_tree_list[-1].GetEntries()

        merge_count += 1
        if (dir_name + 1) % 50 == 0:
            print "Processed %i files" % (dir_name + 1)
            print "N Events: %i" % n_events
        if merge_count % 400 == 0:
            break
        dir_name += 1

    if append:
        last_file = ROOT.TFile(temp_dir + "temp" + final_name)
        data_tree = last_file.Get("data")
        last_n_events = last_file.Get("NEvents").GetUniqueID()
        
        t_list.Add(data_tree)
        n_events += last_n_events

    merge_and_write(t_list, n_events, temp_dir + "temp2" + final_name)
    if append:
        last_file.Close()
    [f.Close() for f in file_list]
    os.rename(temp_dir + "temp2" + final_name, temp_dir + "temp" + final_name)
    if last_failed:
        return (0, False, n_data)
    else: 
        return (dir_name+1, True, n_data)

def merge_and_write(t_list, n_events, f_out_name):
    print "MERGING"
    f_out = ROOT.TFile(f_out_name, "RECREATE")
    f_out.cd()
    merged = ROOT.TTree.MergeTrees(t_list)
    merged.SetName("data")

    print "WRITING TO FILE"
    merged.Write()
    n_events_obj = ROOT.TNamed("NEvents", "NEvents")
    n_events_obj.SetUniqueID(n_events)
    n_events_obj.Write()
    f_out.Close()

def main(datadir, out_file_name, temp_dir):
    append = False
    dir_name = 0
    n_data = 0
    final_name = out_file_name.split("/")[-1]
    while 1:
        (dir_name, append, n_data) = loop_merge(dir_name, datadir, temp_dir, append, n_data, final_name)
        if not append:
            break

    os.rename(temp_dir + "temp" + final_name, out_file_name)
    print "N DATA: %i" % n_data

if __name__=="__main__":
    datadir = sys.argv[1]
    main(datadir, sys.argv[2], sys.argv[3])

