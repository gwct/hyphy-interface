############################################################
# Functions for the FEL model from Hyphy.
# 11.2020
############################################################

import sys, os, json, re, hpcore, hpseq, hptree as tp
from collections import defaultdict

############################################################

def generate(indir, tree_input, gt_opt, gt_extension, aln_id_delim, hyphy_path, outdir, logdir, outfile):
    if aln_id_delim:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : f.split(aln_id_delim)[0], "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    else:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : "NA", "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    # Read and sort the alignment file names

    i = 0;
    tree_skipped, stop_skipped = 0, 0;
    for aln in aligns:

        i += 1;
        print(str(i) + "\t" + aln);

        if gt_opt:
            if aln_id_delim:
                tree_dir = os.path.join(tree_input, aln);
                if os.path.isdir(tree_dir):
                    tree_dir_files = os.listdir(tree_dir);
                    tree_files = [ f for f in tree_dir_files if re.findall(aligns[aln]['id'] + '(.*)' + gt_extension, f) != [] and "rooted" not in f ];
                    if len(tree_files) != 1:
                        print(tree_files)
                        outfile.write(" # Multiple trees found. Skipping: " + aln + "\n");
                        tree_skipped += 1;
                        continue;
                    tree_file = "".join(tree_files);
                    if not tree_file:
                        continue;
                    tree_file = os.path.join(tree_dir, tree_file);
                else:
                    tree_file = False;
            # If we need to split the input alignment directory to get the tree file name
            else:
                tree_file = os.path.join(tree_input, aln, aln + gt_extension);
            # Get the tree file name

            if os.path.isfile(tree_file):
                aligns[aln]['tree'] = tree_file;
            # Assign the current tree file to the alignment
            if not aligns[aln]['tree']:
                print(aligns[aln]['tree'])
                outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
                tree_skipped += 1;
                continue;
            # Check the tree file.

        else:
            aligns[aln]['tree'] = tree_input;
        # If a single tree is input, we need to set the tree and targets to whatever was determined above in the
        # single tree block.        

        seq_dict = hpseq.fastaGetDict(aligns[aln]['aln-file']);
        prem_stop_flag = False
        for title in seq_dict:
            prem_stop, new_seq = hpseq.premStopCheck(seq_dict[title], allowlastcodon=True, rmlast=True);
            if prem_stop:
                prem_stop_flag = True;
            seq_dict[title] = new_seq;
        # Read the sequence and check for premature stop codons.

        if prem_stop_flag:
            outfile.write(" # Premature stop found. Skipping: " + aln + "\n");
            stop_skipped += 1;
            continue;
        # Check the alignment for premature stop codons (which PAML hangs on)

        # cur_outdir = os.path.join(outdir, aln);
        # if not os.path.isdir(cur_outdir):
        #     os.system("mkdir " + cur_outdir);
        # Make the output directory for this alignment

        cur_jsonfile = os.path.join(outdir, aln + ".json");
        #cur_outfile = os.path.join(cur_outdir, align + "-out.txt");
        cur_logfile = os.path.join(logdir, aln + ".log");
        # Get the control and output file names

        hyphy_cmd = "hyphy fel --alignment " + aligns[aln]['aln-file'] + " --tree " +  aligns[aln]['tree'] + " --output " + cur_jsonfile + " &> " + cur_logfile 
        outfile.write(hyphy_cmd + "\n");
        # Construct and write the hyphy command

    hpcore.PWS("# Num skipped because tree file not found     : " + str(tree_skipped), outfile);
    hpcore.PWS("# Num skipped because of premature stop codons: " + str(stop_skipped), outfile);

############################################################

def parse(indir, features, outfile, pad):

    outdir = os.path.join(indir, "csv");
    if not os.path.isdir(outdir):
        os.system("mkdir " + outdir);
    # Make the csv directory if it doesn't exist
    ##########################

    if features:
        summary_headers = ["file","id","chr","start","end","num.sites","num.sites.alpha.ne.beta","num.sites.ps"];
    else:
        summary_headers = ["file","num.sites","num.sites.alpha.ne.beta","num.sites.ps"];
    
    
    # Write the output headers 
    ##########################

    hyphy_files = os.listdir(indir);
    num_files = len(hyphy_files);
    num_files_str = str(num_files);
    num_files_len = len(num_files_str);

    # Read align file names from input directory
    ##########################

    num_unfinished = 0;
    # A count of the number of unfinished hyphy runs as determined by empty output files

    init_alpha = 0.01;
    # The initial significance level

    total_sites_tested = 0;

    print("Counting number of tests...");

    for f in os.listdir(indir):
        cur_json_file = os.path.join(indir, f);
        if not os.path.isfile(cur_json_file) or not cur_json_file.endswith(".json"):
            continue;
        if os.stat(cur_json_file).st_size == 0:
            continue;
        # Get the current json file.   

        with open(cur_json_file) as json_data:
            cur_data = json.load(json_data);
        # Read the Hyphy json file.

        for site in cur_data["MLE"]["content"]["0"]:
            if site[3] != 0:
                total_sites_tested += 1;

    adj_alpha = init_alpha / total_sites_tested;

    hpcore.PWS("# Initial alpha      : " + str(init_alpha), outfile);
    hpcore.PWS("# Total sites tested : " + str(total_sites_tested), outfile);
    hpcore.PWS("# Adjusted alpha     : " + str(adj_alpha), outfile);

    outfile.write(",".join(summary_headers) + "\n");

    ##########################

    counter = 0;
    for f in os.listdir(indir):
        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_files_len:
                counter_str = "0" + counter_str;
            print ("> " + hpcore.getDateTime() + " " + counter_str + " / " + num_files_str);
        counter += 1;
        # Track progress

        ##########################

        cur_json_file = os.path.join(indir, f);
        if not os.path.isfile(cur_json_file) or not cur_json_file.endswith(".json"):
            continue;
        if os.stat(cur_json_file).st_size == 0:
            num_unfinished +=1 ;
            continue;
        # Get the current json file.

        ##########################

        #print(f);
        if features:
            if "-" in f:
                fid = f.split("-")[0];
            elif "." in f:
                fid = f.split(".")[0];
            else:
                fid = f;
            cur_feature = features[fid];
            # Look up transcript/gene info for current gene to save in output, if metadata is provided.

        ##########################

        with open(cur_json_file) as json_data:
            cur_data = json.load(json_data);
        # Read the Hyphy json file.

        #print(cur_data)

        if features:
            summary_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                "num.sites" : "NA", "num.sites.alpha.ne.beta" : "NA", "num.sites.ps" : "NA" };
        else:
            summary_info = { "num.sites" : "NA", "num.sites.alpha.ne.beta" : "NA", "num.sites.ps" : "NA" };   
        # Initialize the output dictionary for the current branch.

        ##########################

        num_sites, num_sites_tested, num_sites_ps = 0, 0, 0;

        for site in cur_data["MLE"]["content"]["0"]:
            num_sites += 1;
            if site[3] != 0:
                num_sites_tested += 1;
                if float(site[4]) < adj_alpha:
                    num_sites_ps += 1;

        summary_info['num.sites'] = str(num_sites);
        summary_info['num.sites.alpha.ne.beta'] = str(num_sites_tested);
        summary_info['num.sites.ps'] = str(num_sites_ps);

        ##########################

        summary_outline = [f] + [ summary_info[h] for h in summary_headers if h not in ["file"] ];
        outfile.write(",".join(summary_outline) + "\n");
        # Ouput rate estimates for both the gene.

        ##########################

    hpcore.PWS("# ----------------", outfile);
    #hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);