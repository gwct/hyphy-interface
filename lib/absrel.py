############################################################
# Functions for the aBSREL model from Hyphy.
# 12.2020
############################################################

import sys, os, json, re, lib.hpcore as hpcore, lib.hpseq as hpseq, lib.hptree as tp
from collections import defaultdict

############################################################

def assignTargetClade(aln, tree_file, targets, target_summary_dict, target_gene_dict):
# Function that reads a tree into a dictionary and determines the target branch from a set of tip labels

    tree_str = open(tree_file, "r").read().strip();
    tinfo, t, root = tp.treeParse(tree_str);
    # Parse tree into dictionary
    
    target_labels = [];
    # A list of internal node labels that correspond to target clades

    targets_found = False;
    # A boolean to tell whether any targets were found in the current gene tree

    for target in targets:
        max_split1 = [];
        max_split2 = [];
        max_node = "";
        # Variables to store the information for the maximal subset clade

        for node in tinfo:
            split1 = tp.getClade(node, tinfo);
            split1.sort();
            split2 = [ n2 for n2 in tinfo if tinfo[n2][2] == 'tip' and n2 not in split1 ];
            split2.sort();

            splits = [split1, split2];
            # Get the splits from the current branch in the gene's gene tree

            for s in range(len(splits)):
            # For each split, we want to check whether it is a subset of the current branch

                if set(splits[s]) == set(targets[target]) or set(splits[s]) <= set(targets[target]):
                # Check if the current split is an equivalent set or subset of the targets

                    if len(splits[s]) > len(max_split1):
                    # Check if the current split that is a subset of the targets contains more species
                    # than the current maximal subset. If so, it becomes the maximal subset.

                        max_split1 = splits[s];
                        max_node = node;
                        if s == 0:
                            max_split2 = splits[1];
                        else:
                            max_split2 = splits[0];
                        # Assign the maximal subset to max_split1 and the other split to max_split2
        # Check every node in the tree for target species

        if max_split1 == [] or any(spec in max_split2 for spec in targets[target]):
            target_summary_dict[target] += 1;
        # If no species from the current branch are found then this clade doesn't exist in this gene OR
        # If species from the branch are present in the second split (not the maximal subset), then this is a case of
        # discordance and the clade truly doesn't exist in this gene as a monophyly. Do not add to the target labels
        # dictionary, but rather increment the summary dict that counts how many genes DON'T contain a clade.

        else:
            target_gene_dict[target][aln] = { 'node' : max_node, 'num-species' : str(len(max_split1)), 'species' : max_split1 };
            # Add this node for this target to the target_gene dict.

            target_labels.append(max_node);
            targets_found = True;
        # Otherwise save the internal node that label that corresponds to the clade and switch targets_found to True

    return t, target_labels, targets_found, target_summary_dict, target_gene_dict;

############################################################

def generate(indir, tree_input, gt_opt, gt_extension, aln_id_delim, target_clades, hyphy_path, outdir, target_tree_dir, logdir, outfile):
    if aln_id_delim:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : f.split(aln_id_delim)[0], "tree" : False, "targets" : [] } for f in os.listdir(indir) if f.endswith(".fa") };
    else:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : "NA", "tree" : False, "targets" : [] } for f in os.listdir(indir) if f.endswith(".fa") };
    # Read and sort the alignment file names

    if target_clades:
        target_summary = { t : 0 for t in target_clades };
        target_gene_key = { t : {} for t in target_clades };
        # Two dictionaries to keep track of target stats:
        # target_summary is just a count of how many genes DON'T have each target branch: <input target subtree filename> : <# of genes without subtree>
        # target_gene_dict keeps track of which node corresponds to each target branch in each clade: <input target subtree filename> : { <alignment id> : { <node that corresponds to target clade>, <# of species in target clade found> < list of species in target clade found> }}

    if target_clades and not gt_opt:
        target_tree, target_labels, targets_found, target_summary, target_gene_key = assignTargetClade("all", tree_input, target_clades, target_summary, target_gene_key);
        # Assign the target labels for the input tree
        if not targets_found:
            sys.exit(" * ERROR: None of the provided target clades exist in the input tree.");
        # If the target branch doesn't exist in the single input tree, exit
        else:
            tree_input = os.path.basename(tree_input);
            tree_input = os.path.join(target_tree_dir, tree_input);
        # Write the labeled input tree for the run
    # This block handles the target branch identification for a single input tree

    no_target_trees = 0;
    # A count of the number of trees that don't have the target branch, for gene tree input
    i = 0;
    tree_skipped, stop_skipped = 0, 0;
    for aln in aligns:
        # if i == 10:
        #     break;
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

            if target_clades:
                target_tree, target_labels, targets_found, target_summary, target_gene_key = assignTargetClade(aln, tree_file, target_clades, target_summary, target_gene_key);
                # Identify the target nodes in the current gene tree

                if not targets_found:
                    outfile.write(" # No target branches found in input tree. Skipping: " + tree_file + "\n");
                    no_target_trees += 1;
                    continue;
                # Count the tree if none of the target clades are found
                else:
                    tree_file = os.path.basename(tree_file).replace(".treefile", "-labeled.treefile");
                    tree_file = os.path.join(target_tree_dir, tree_file);
                    for target_label in target_labels:
                        target_tree = target_tree.replace(target_label, target_label + "{Foreground}");
                    # Add the "{Foreground}" label to the test branches
                    with open(tree_file, "w") as tout:
                        tout.write(target_tree);
                    aligns[aln]['tree'] = tree_file;
                    aligns[aln]['targets'] = target_labels;
                # Re-write the labeled tree file if any target clades are found
            # If target clades are provided, check if they are in the current tree
        # This handles individual gene trees for each alignment

        else:
            aligns[aln]['tree'] = tree_input;
            aligns[aln]['targets'] = target_labels;
        # If a single tree is input, we need to set the tree and targets to whatever was determined above in the
        # single tree block.

        # print(target_gene_key);
        # sys.exit();



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

        cur_jsonfile = os.path.join(outdir, aln + ".json");
        #cur_outfile = os.path.join(cur_outdir, align + "-out.txt");
        cur_logfile = os.path.join(logdir, aln + ".log");
        # Get the control and output file names

        hyphy_cmd = "hyphy absrel --alignment " + aligns[aln]['aln-file'] + " --tree " +  aligns[aln]['tree'];
        if aligns[aln]['targets']:
            hyphy_cmd += " --branches Foreground";
        # Add the target node if specified and found
        hyphy_cmd += " --output " + cur_jsonfile + " &> " + cur_logfile 
        outfile.write(hyphy_cmd + "\n");
        # Construct and write the hyphy command

    hpcore.PWS("# END CMDS", outfile);
    hpcore.PWS("# ----------------", outfile);
    hpcore.PWS("# Num skipped because tree file not found     : " + str(tree_skipped), outfile);
    hpcore.PWS("# Num skipped because of premature stop codons: " + str(stop_skipped), outfile);
    if target_clades:
        hpcore.PWS("# TARGET SUMMARY", outfile);
        hpcore.PWS("# Num trees with no target branches           : " + str(no_target_trees), outfile);
        summary_pad = 18;

        hpcore.PWS(hpcore.spacedOut("# Clade num", summary_pad) + hpcore.spacedOut("trees not found", summary_pad) + "clade list", outfile);
        for clade in target_clades:
            outline = hpcore.spacedOut("# " + str(clade), summary_pad) + hpcore.spacedOut(str(target_summary[clade]), summary_pad) + ",".join(target_clades[clade]);
            hpcore.PWS(outline, outfile);

        target_key_file = os.path.join(logdir, "target-gene-keys.tab");
        hpcore.PWS(hpcore.spacedOut("# Writing target gene key file: ", summary_pad) + target_key_file, outfile);
        with open(target_key_file, "w") as keyfile:
            for target in target_gene_key:
                for aln in target_gene_key[target]:
                    outline = "\t".join([target, aln, target_gene_key[target][aln]['node'], target_gene_key[target][aln]['num-species'], ",".join(target_gene_key[target][aln]['species'])]);
                    keyfile.write(outline + "\n");
        # This block writes the target_gene_dict info to a file.
    # Report some stats

############################################################

def parseOLD(indir, features, outfile, pad):

    if features:
        headers = ["file","id","chr","start","end","num ps branches", "ps pvals"];
    else:
        headers = ["file","branches","num ps branches","ps pvals"];
    outfile.write(",".join(headers) + "\n");
    # Write the output headers 

    hyphy_files = os.listdir(indir);
    num_files = len(hyphy_files);
    num_files_str = str(num_files);
    num_files_len = len(num_files_str);
    # Read align file names from input directory

    num_unfinished = 0;
    # A count of the number of unfinished hyphy runs as determined by empty output files

    counter = 0;
    for f in os.listdir(indir):
        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_files_len:
                counter_str = "0" + counter_str;
            print ("> " + hpcore.getDateTime() + " " + counter_str + " / " + num_files_str);
        counter += 1;
        # Track progress

        cur_json_file = os.path.join(indir, f);
        if not os.path.isfile(cur_json_file) or not cur_json_file.endswith(".json"):
            continue;
        if os.stat(cur_json_file).st_size == 0:
            num_unfinished +=1 ;
            continue;
        # Get the current output file.

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

        with open(cur_json_file) as json_data:
            cur_data = json.load(json_data);
        # Read the Hyphy json file.

        #print(cur_data)

        if features:
            gene_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                "num ps branches" : 0, "ps pvals" : [] };
        else:
            gene_info = { "branches" : [], "num ps branches" : 0, "ps pvals" : [] };   
        # Initialize the output dictionary for the current branch.

        #gene_info["dn/ds"] = str(cur_data["fits"]["Standard MG94"]["Rate Distributions"]["non-synonymous/synonymous rate ratio"]);
        for node in cur_data["branch attributes"]["0"]:
            if float(cur_data["branch attributes"]["0"][node]["Corrected P-value"]) < 0.01:
                gene_info["branches"].append(str(node)); 
                gene_info["num ps branches"] += 1;
                gene_info["ps pvals"].append(str(cur_data["branch attributes"]["0"][node]["Corrected P-value"]));
        # Retrieve the rate estimates from the json data.

        gene_info["branches"] = ";".join(gene_info["branches"]);
        gene_info["ps pvals"] = ";".join(gene_info["ps pvals"]);
        gene_info["num ps branches"] = str(gene_info["num ps branches"]);
        gene_outline = [f] + [ gene_info[h] for h in headers if h not in ["file"] ];
        outfile.write(",".join(gene_outline) + "\n");
        # Ouput rate estimates for both the gene.

    hpcore.PWS("# ----------------", outfile);
    #hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);


############################################################


def parse(indir, features, outfile, pad):

    alpha = 0.01;
    print("alpha: ", alpha);

    out_mode = "splits";
    print("out mode: " + out_mode);
    # clades or splits

    ##########################

    outdir = os.path.join(indir, "csv");
    if not os.path.isdir(outdir):
        os.system("mkdir " + outdir);
    # Make the csv directory if it doesn't exist
    ##########################

    if features:
        if out_mode == "splits":
            gene_headers = ["file","id","chr","start","end","baseline mean omega","baseline median omega","num branches","num branches pval less than alpha"];
            branch_headers = ["file","id","chr","start","end","branch","split1","split2","baseline omega","pval","lrt","rate classes"];
        elif out_mode == "clades":
            gene_headers = ["file","id","chr","start","end","baseline mean omega","baseline median omega","num branches","num branches pval less than alpha"];
            branch_headers = ["file","id","chr","start","end","branch","clade","baseline omega","pval","lrt","rate classes"];
    else:
        if out_mode == "splits":
            gene_headers = ["file","baseline mean omega","baseline median omega","num branches","num branches pval less than alpha"];
            branch_headers = ["file","branch","split1","split2","baseline omega","pval","lrt","rate classes"];
        elif out_mode == "clades":
            gene_headers = ["file","baseline mean omega","baseline median omega","num branches","num branches pval less than alpha"];
            branch_headers = ["file","branch","clade","baseline omega","pval","lrt","rate classes"];
    outfile.write(",".join(gene_headers) + "\n");

    # Write the output headers 
    ##########################


    hyphy_files = os.listdir(indir);
    num_files = len(hyphy_files);
    num_files_str = str(num_files);
    num_files_len = len(num_files_str);

    print("num files: ", num_files);

    adj_alpha = alpha / num_files;
    print("adjusted alpha: ", adj_alpha)

    # Read align file names from input directory
    ##########################

    num_unfinished = 0;
    # A count of the number of unfinished hyphy runs as determined by empty output files

    num_skipped = 0;
    # Number of files skipped because the tree is too small

    num_sig_branches = 0;
    num_sig_genes = 0;

    counter = 0;
    for f in os.listdir(indir):
    # Loop over every file in the input directory

        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_files_len:
                counter_str = "0" + counter_str;
            print ("> " + hpcore.getDateTime() + " " + counter_str + " / " + num_files_str);
        counter += 1;
        # Track progress

        ##########################

        #print(f);

        cur_json_file = os.path.join(indir, f);
        if not os.path.isfile(cur_json_file) or not cur_json_file.endswith(".json"):
            continue;
        if os.stat(cur_json_file).st_size == 0:
            num_unfinished +=1 ;
            continue;
        # Get the current json file.

        ##########################

        with open(cur_json_file) as json_data:
            cur_data = json.load(json_data);
        # Read the Hyphy json file.

        cur_tree = cur_data['input']['trees']['0'];
        tinfo, tree, root = tp.treeParse(cur_tree);
        # Read the tree from the json data.

        tips = [ n for n in tinfo if tinfo[n][2] == 'tip' ];
        if len(tips) < 3:
            num_skipped += 1;
            continue;
        # If the current tree has fewer than 3 tips, skip

        ##########################

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

        gene_info = { h : 0 for h in gene_headers };
        # The headers for the current gene in the main output file

        gene_info["file"] = f;
        gene_info["baseline mean omega"] = cur_data["fits"]["Baseline MG94xREV"]["Rate Distributions"]["Per-branch omega"]["Mean"];
        gene_info["baseline median omega"] = cur_data["fits"]["Baseline MG94xREV"]["Rate Distributions"]["Per-branch omega"]["Median"];
        # Lookup the gene info in the json

        ##########################

        sig = False;

        branch_outfile = os.path.join(outdir, f.replace(".json", ".csv"));
        # Name the csv output file for the current gene

        with open(branch_outfile, "w") as goutfile:
        # Open the current gene's output csv file

            goutfile.write(",".join(branch_headers) + "\n");
            # Write the headers for the current csv file

            splits = defaultdict(list);
            # A dictionary that holds the splits for each branch
            # <key> : <value> = <hyphy node label> : <[split1, split2]>
            for n in tinfo:
                if tinfo[n][2] == 'root':
                    continue;
                if tinfo[n][2] == 'tip':
                    node_label = n;
                else:
                    node_label = tinfo[n][3];
                # Skip the root, for tips use the tip label, for internal nodes use the hyphy label

                split1 = sorted(tp.getClade(n, tinfo), key=str.casefold);
                # Get the current clade and sort it

                split2 = [ n2 for n2 in tinfo if tinfo[n2][2] == 'tip' and n2 not in split1 ];
                split2 = sorted(split2, key=str.casefold);
                # Get all other tips as the second split and sort it

                splits[node_label].append(";".join(split1));
                splits[node_label].append(";".join(split2));
                # For unrooted trees there is no directionality, so get clades from both sides of each branch.
            # For each node/branch in the tree, save the tip species that make up the clade.

            ##########################

            branch_info = {};
            for node in splits:
                #print("node1: " + node);
                if features:
                    if out_mode == "splits":
                        branch_info[node] = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                            "branch" : node, "split1" : "NA", "split2" : "NA", 'baseline omega' : "NA", 'pval' : "NA", 'lrt' : "NA", 'rate classes' : "NA" };
                    elif out_mode == "clades":
                        branch_info[node] = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                            "branch" : node, "clade" : "NA", 'baseline omega' : "NA", 'pval' : "NA", 'lrt' : "NA", 'rate classes' : "NA" };                        
                else:
                    if out_mode == "splits":
                        branch_info[node] = { "branch" : node ,"split1" : "NA", "split2" : "NA", 'baseline omega' : "NA", 'pval' : "NA", 'lrt' : "NA", 'rate classes' : "NA" };   
                    elif out_mode == "clades":
                        branch_info[node] = { "branch" : node ,"clades" : "NA", 'baseline omega' : "NA", 'pval' : "NA", 'lrt' : "NA", 'rate classes' : "NA" };   
            # Initialize the output dictionary for the current branch.
            # TODO: Generate via comprehension

            ##########################

            for node in cur_data["branch attributes"]["0"]:
            # Loop over every node in the tree in the json data
                #print("node2: " + node);
                #print(cur_data["branch attributes"]["0"][node]);

                gene_info["num branches"] += 1;

                branch_info[node]['baseline omega'] = cur_data["branch attributes"]["0"][node]["Baseline MG94xREV omega ratio"];
                branch_info[node]['pval'] = cur_data["branch attributes"]["0"][node]["Corrected P-value"];
                if cur_data["branch attributes"]["0"][node]["Corrected P-value"] <= adj_alpha:
                    sig = True;
                    num_sig_branches += 1;
                    gene_info["num branches pval less than alpha"] += 1;
                    
                branch_info[node]['lrt'] = cur_data["branch attributes"]["0"][node]["LRT"];
                branch_info[node]['rate classes'] = cur_data["branch attributes"]["0"][node]["Rate classes"];

                # dist_num = 1;
                # for rate_dist in cur_data["branch attributes"]["0"][node]["Rate Distributions"]:
                #     print(rate_dist);
                
                if out_mode == "splits":
                    branch_info[node]["split1"] = splits[node][0];
                    branch_info[node]["split2"] = splits[node][1];

                    branch_outline = [f] + [ str(branch_info[node][h]) for h in branch_headers if h not in ["file"] ];
                    #outfile.write(",".join(branch_outline) + "\n");
                    goutfile.write(",".join(branch_outline) + "\n");

                if out_mode == "clades":
                    for clade in splits[node]:
                        branch_info[node]["clade"] = clade;
                        branch_outline = [f] + [ branch_info[node][h] for h in branch_headers if h not in ["file"] ];
                        #outfile.write(",".join(branch_outline) + "\n");
                        goutfile.write(",".join(branch_outline) + "\n");
                # Ouput rate estimates for both clades of the current branch.
            ## End json node loop
            ##########################
        ## Close gene csv file

        if sig:
            num_sig_genes += 1;
        gene_outline = [ str(gene_info[h]) for h in gene_headers ];
        outfile.write(",".join(gene_outline) + "\n");
        # Write the output for the current gene

    ## End file loop

    hpcore.PWS("# ----------------", outfile);
    hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);
    hpcore.PWS(hpcore.spacedOut("# Number skipped with only 2 species:", pad) + str(num_skipped), outfile);
    hpcore.PWS(hpcore.spacedOut("# Number branches below alpha:", pad) + str(num_sig_branches), outfile);
    hpcore.PWS(hpcore.spacedOut("# Number genes with at least 1 branch below alpha:", pad) + str(num_sig_genes), outfile);
    # Report some stats

############################################################












































'''
def generate(indir, tree_input, gt_opt, targets, aln_id_delim, hyphy_path, outdir, logdir, outfile):

    if aln_id_delim:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : f.split(aln_id_delim)[0], "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    else:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : "NA", "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    # Read and sort the alignment file names

    for aln in aligns:
        if gt_opt:
            tree_file = os.path.join(tree_input, aln, aln + "-rooted.treefile");
        else:
            tree_file = tree_input;

        if os.path.isfile(tree_file):
            cur_tree = open(tree_file, "r").read().strip();
            cur_tree = re.sub("\)[\de.-]+:", "):", cur_tree);
            aligns[aln]['tree'] = cur_tree;
        # Read the tree and remove any bootstrap node labels.
    # Read the appropriate tree depending on the -tree and -genetree options.

    tree_skipped, target_skipped, stop_skipped = 0, 0, 0;
    for aln in aligns:
        if not aligns[aln]['tree']:
            outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
            tree_skipped += 1;
            continue;
        # Check the tree file.          

        target_found, cur_tree = assignTargetClade(aligns[aln]['tree'], targets);
        if not target_found:
            outfile.write(" # Target clade not found. Skipping: " + aligns[aln]['aln-file'] + "\n");
            target_skipped += 1;
            continue;
        # Assign the target lineages to the current tree.

        seq_dict = pseq.fastaGetDict(aligns[aln]['aln-file']);
        prem_stop_flag = False
        for title in seq_dict:
            prem_stop, new_seq = pseq.premStopCheck(seq_dict[title], allowlastcodon=True, rmlast=True);
            if prem_stop:
                prem_stop_flag = True;
            seq_dict[title] = new_seq;
        # Read the sequence and check for premature stop codons.

        if prem_stop_flag:
            outfile.write(" # Premature stop found. Skipping: " + aln + "\n");
            stop_skipped += 1;
            continue;
        # Check the alignment for premature stop codons (which PAML hangs on)

        cur_outdir = os.path.join(outdir, aln);
        if not os.path.isdir(cur_outdir):
            os.system("mkdir " + cur_outdir);
        # Make the output directory for this alignment

        new_treefile = os.path.join(cur_outdir, "codeml.tre");
        with open(new_treefile, "w") as treefile:
            treefile.write(cur_tree);
        # Make the tree file for this alignment

        new_seqfile = os.path.join(cur_outdir, "codeml.fa");
        with open(new_seqfile, "w") as seqfile:
            for title in seq_dict:
                seqfile.write(title + "\n");
                seqfile.write(seq_dict[title] + "\n");
        # Write the sequences for this alignment

        cur_ctlfile = os.path.join(cur_outdir, "codeml.ctl");
        cur_outfile = os.path.join(cur_outdir, "codeml.out");
        #cur_logfile = os.path.join(cur_outdir, base_input + "-codeml.log");
        # Get the control and output file names

        with open(cur_ctlfile, "w") as ctlfile:
            ctlfile.write(ctlfile_template.format(infile=new_seqfile, treefile=new_treefile, outfile=cur_outfile));
        # Write the control file

        codeml_cmd = "cd " + cur_outdir + "; " + paml_path + " codeml.ctl";
        outfile.write(codeml_cmd + "\n");
        # Construct and write the codeml command

    pcore.PWS("# Num skipped because tree file not found     : " + str(tree_skipped), outfile);
    pcore.PWS("# Num skipped because target clade not found  : " + str(target_skipped), outfile);
    pcore.PWS("# Num skipped because of premature stop codons: " + str(stop_skipped), outfile);

############################################################

def parse(indir, features, outfile, pad):

    if features:
        headers = ["file","id","chr","start","end","lnl","k","dn sum","ds sum","dn avg","ds avg","clade","clade dn/ds","clade dn","clade ds"];
    else:
        headers = ["file","lnl","k","dn sum","ds sum","dn avg","ds avg","clade","clade dn/ds","clade dn","clade ds"];
    outfile.write(",".join(headers) + "\n");
    # Write the output headers 

    num_unfinished = 0;
    # A count of the number of unfinished codeml runs as determined by incomplete output files

    num_nonmono = 0;
    # A count of the number of trees where the target clade is not monophyletic.

    codeml_dirs = os.listdir(indir);
    num_dirs = len(codeml_dirs);
    num_dirs_str = str(num_dirs);
    num_dirs_len = len(num_dirs_str);
    # Read align file names from input directory

    counter = 0;
    for d in os.listdir(indir):
        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_dirs_len:
                counter_str = "0" + counter_str;
            print ("> " + pcore.getDateTime() + " " + counter_str + " / " + num_dirs_str);
        counter += 1;
        # Track progress

        cur_dir = os.path.join(indir, d);
        # Get the current gene ID and codeml directory.

        if features:
            if "-" in d:
                fid = d.split("-")[0];
            elif "." in d:
                fid = d.split(".")[0];
            else:
                fid = d;
            cur_feature = features[fid];
            # Look up transcript/gene info for current gene to save in output

            gene_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                'lnl' : "NA", 'k' : "NA", 'dnds' : "NA", 'dn sum' : "NA", 'ds sum' : "NA", 'dn avg' : "NA", 'ds avg' : "NA" };
        else:
            gene_info = { 'lnl' : "NA", 'k' : "NA", 'dnds' : "NA", 'dn sum' : "NA", 'ds sum' : "NA", 'dn avg' : "NA", 'ds avg' : "NA" }
        # Add meta headers to the output dict if metadata is provided

        branch_info = {};
        # Initialize the output dictionary for the current gene

        read_paml_ancs = False;
        read_paml_tree = False;
        read_orig_tree = False;
        # Flags to indicate when we've reached certain points in the PAML output

        num_branches = 0.0;
        # Count the number of branches in the tree to calculate per branch averages

        paml_ancs = {};
        # This will be the main tree info dictionary, akin to the treeParse dictionary
        # <PAML node label> : [<branch length>, <ancestral node>, <node type>, <original tip node label>]

        mono_skipped = False;

        cur_codeml_file = os.path.join(cur_dir, "codeml.out");
        for line in open(cur_codeml_file):
            if line.startswith("lnL"):
                line = list(filter(None, line.strip().split(" ")));
                gene_info['lnl'] = line[4];

                read_paml_ancs = True;
                continue;
            # This line indicates that the next line will contain the PAML ancestral relationship labels (ANC..DESC)

            if read_paml_ancs:
            # PAML ancestral relationship labels line.

                nodes = list(filter(None, line.strip().split(" ")));
                # Get the relationships into a list

                for node in nodes:
                    anc, n = node.split("..");
                    paml_ancs[n] = ['', anc, '', ''];
                # Read through every PAML ANC..DESC relationship and add it to the main dictionary.

                for node in nodes:
                    anc, n = node.split("..");
                    if anc not in paml_ancs:
                        paml_ancs[anc] = ['NA', "NA", 'root', ''];
                        root = anc;
                # The trifurcating 'root' node will not be included in the loop above (since it has no ancestor). Add it to the
                # dictionary as 'root' here.

                read_paml_ancs = False;
                # Need to indicate that we don't want to run this block on the subsequent lines


            if line.startswith("tree length ="):
                read_paml_tree = True;
                continue;
            # This line indicates that the next line will contain the tree with PAML tip labels

            if read_paml_tree and line != "\n":
            # Line with tree with PAML tip labels

                for node in paml_ancs:
                # Parse every node in the main dictionary

                    if node + ": " in line:
                    # If the node has a branch length in the tree string, do the following

                        paml_ancs[node][2] = 'tip'
                        # Add this node as a tip in the main dictionary

                        #node_ind = line.index(node + ": ");

                        bl_match = re.findall(node + ': [\d.]+\)', line);
                        if bl_match == []:
                            bl_match = re.findall(node + ': [\d.]+,', line);
                        # Retrieve the branch length for this node based on pattern "<NODE> :<BL>)" OR "<NODE> :<BL>,"
                        
                        cur_bl = bl_match[0];
                        cur_bl = cur_bl[cur_bl.index(" ")+1:-1];
                        paml_ancs[node][0] = cur_bl;
                        # Parse the branch length retrieved and add it to the dictionary

                    elif paml_ancs[node][2] != 'root':
                        paml_ancs[node][2] = 'internal'
                    # If the node isn't a tip or the trifurcating 'root', add it as internal with no branch length saved

                paml_label_order = [ l[:l.index(":")] for l in re.findall('[\d]+: ', line) ];
                # Get all the tip node labels in the order they appear in the tree string

                read_paml_tree = False;
                read_orig_tree = True;
                continue;
                # Indicate that we don't want to run this block again, and that the next line contains the tree with 
                # the original tip labels
                
            if read_orig_tree and line != "\n":
            # Line with the tree with the original tip labels
                #print(line);
                orig_label_order = [ l[:l.index(":")] for l in re.findall('[\w]+: ', line) ];
                # Get the tip node labels in the order they appear in the tree string

                assert len(paml_label_order) == len(orig_label_order), "\nUnequal number of labels: " + d + "\nPAML: " + ",".join(paml_label_order) + "\nOrig:" + ",".join(orig_label_order)
                # Make sure we have the same number of labels in the PAML and original trees

                for n in range(len(paml_label_order)):
                    node = paml_label_order[n];
                    paml_ancs[node][3] = orig_label_order[n];
                # Add the original label into the label field of the main dictionary

                max_clade_node, max_clade_count = "", 0;
                for n in paml_ancs:
                    if paml_ancs[n][2] != 'tip':
                        cur_clade_paml = ptree.getClade(n, paml_ancs);
                        cur_clade_orig = sorted([ paml_ancs[tip][3] for tip in cur_clade_paml ]);
                        paml_ancs[n][3] = cur_clade_orig;

                        if paml_ancs[n][2] != 'root' and len(cur_clade_paml) > max_clade_count:
                            max_clade_node, max_clade_count = n, len(cur_clade_paml);
                        # Get the longest clade other than the "root" clade
                # Get the clades for each node

                for n in paml_ancs:
                    if paml_ancs[n][2] == 'tip':
                        paml_ancs[n][3] = [paml_ancs[n][3]];
                # Convert the tip labels into lists    

                # clade_len_counts = [ len(paml_ancs[n][3]) for n in paml_ancs ];
                # if clade_len_counts.count(max_clade_count) != 1:
                #     print("More than one max clade length found: ");
                #     print(tid);
                #     print(paml_ancs);
                #     print(max_clade_node);
                #     sys.exit(1);
                #assert clade_len_counts.count(max_clade_count) == 1, "\nMore than one max clade length found: " + gid + "\n" + paml_ancs;

                root_clade = sorted(list(set(paml_ancs[root][3]) - set(paml_ancs[max_clade_node][3])));
                paml_ancs[root][3] = root_clade;
                # Subtract out the longest clade from the 'root' clade to get the remaining clade.

                read_orig_tree = False;
                # Indicate that we don't need to run this block again
                continue;

            if line.startswith("kappa (ts/tv)"):
                line = line.strip().split(" ");
                #print(line);
                gene_info['k'] = line[-1];
                continue;
            # Get the kappa value for the gene

            if line.count("..") == 1 and "check convergence" not in line:
            # These lines report the brance specific results

                line = list(filter(None, line.strip().split(" ")));
                # Convert the line into a list

                nodes = line[0].split("..")

                node = nodes[1];
                anc = nodes[0];

                dnds, dn, ds = line[4], line[5], line[6];
                # Parse relevant info from line

                cur_clade = ";".join(paml_ancs[node][3]);
                anti_clade = ";".join(sorted([ s for s in orig_label_order if s not in paml_ancs[node][3] ]));
                branch_info[cur_clade] = { 'dnds' : dnds, 'dn' : dn, 'ds' : ds };
                branch_info[anti_clade] = { 'dnds' : dnds, 'dn' : dn, 'ds' : ds };
                # Since unrooted trees lack directionality, get both clades on both sides of the branch

                # if node == max_clade_node and anc == root:
                #     anc = max_clade_node;
                #     node = root;
                #     cur_clade = ";".join(paml_ancs[node][3]);
                #     branch_info[cur_clade] = { 'dnds' : dnds, 'dn' : dn, 'ds' : ds };
                # For the branch between the 'root' and its opposite, max clade node, reverse the directonality and output that clade as well

                num_branches += 1;
                continue;

            if line.startswith("tree length for dN:"):
                line = line.strip().split(" ");
                #print(line);
                gene_info['dn sum'] = line[-1];
                continue;
            # Sum of dN for all branches

            if line.startswith("tree length for dS:"):
                line = line.strip().split(" ");
                #print(line);
                gene_info['ds sum'] = line[-1];
                continue;
            # Sum of dS for all branches

        try:
            gene_info['dn avg'] = str(float(gene_info['dn sum']) / num_branches);
            gene_info['ds avg'] = str(float(gene_info['ds sum']) / num_branches);   
        except:
            #print(gid);
            if not mono_skipped:
                num_unfinished += 1;
        # Try to calculate average dN and dS for all branches. If this fails, codeml likely didn't finish

        gene_outline = [d] + [ gene_info[h] for h in headers if h not in ["file","clade","clade dn/ds","clade dn","clade ds"] ];

        # gene_outline = [ gene_info['gid'], gene_info['mtid'], gene_info['mchr'], gene_info['mstart'], gene_info['mend'],
        #             gene_info['ptid'], gene_info['pchr'], gene_info['pscaff'], gene_info['pstart'], gene_info['pend'],
        #             gene_info['lnl'], gene_info['k'], 
        #             gene_info['dn-sum'], gene_info['ds-sum'], gene_info['dn-avg'], gene_info['ds-avg'] ];

        for branch in branch_info:
            branch_outline = [ branch, branch_info[branch]['dnds'], branch_info[branch]['dn'], branch_info[branch]['ds'] ];
            outline = gene_outline + branch_outline;
            outfile.write(",".join(outline) + "\n");
        # Output all info for each branch for current gene

    pcore.PWS("# ----------------", outfile);
    pcore.PWS(pcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);
'''
