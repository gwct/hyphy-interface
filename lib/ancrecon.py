############################################################
# Functions for the free-ratio (MG94) model from Hyphy.
# 11.2020
############################################################

import sys, os, json, re, lib.hpcore as hpcore, lib.hpseq as hpseq, lib.hptree as tp
from collections import defaultdict

############################################################

def generate(indir, model_file, hyphy_path, outdir, logdir, outfile):
    for f in os.listdir(indir):
        if f.endswith(".fit"):
            infilename = os.path.join(indir, f);
            outfilename = os.path.join(outdir, f.replace(".fit", ".anc"));
            logfilename = os.path.join(logdir, f.replace(".fit", ".log"));

            hyphy_cmd = "hyphy " + model_file + " --fit " + infilename + " --output " + outfilename + " &> " + logfilename;
            outfile.write(hyphy_cmd + "\n");
            # Construct and write the hyphy command

############################################################

def parse(indir, features, outfile, pad):

    out_mode = "splits";
    print("out mode: " + out_mode);
    # clades or splits
    # NOTE: Only splits supported now
    # NOTE: Since these are done on rooted trees, we are essentially unrooting the tree for each branch, but we will
    #       still assume the directionality of the branch in the branch rates script

    outdir = os.path.join(indir, "csv");
    if not os.path.isdir(outdir):
        os.system("mkdir " + outdir);
    # Create the locus output directory within the input directory

    ####################

    #codon_sub_types = ["syn", "non-syn", "one-each", "amb"];
    #codon_sub_types = ["s", "n", "a"];
    # codon_sub_types = ['s-s', 's-n', 's-a', 'n-n', 'n-a', 'a-a'];
    # codon_sub_classes = ["con", "par", "div"];
    # aa_sub_classes = ["con", "par", "div", "none", "na"];
    # Various substitution types and classes to count among pairs of branches

    headers = ["file","branch","node type","clade","split2","total subs","total mns","total s","total n","total a","s mns","n mns","a mns"];
    # Headers for the output file
    
    # for codon_sub_type1 in codon_sub_types:
    #     for codon_sub_class in codon_sub_classes:
    #             for aa_sub_class in aa_sub_classes:
    #                 sub_label = "-".join([codon_sub_type1, codon_sub_class, aa_sub_class]);
    #                 headers.append(sub_label);
    # Generate the headers for each substitution label   

    # if features:
    #     if out_mode == "splits":
    #         #headers = ["file","id","chr","start","end","branch","split1","split2","num subs","num ins","num dels","num mns","sites","subs"];
    #         headers = ["file","id","chr","start","end","branch","split1","split2","num subs","syn subs","nonsyn subs","ambiguous subs","num convergent subs","syn convergent subs","nonsyn convergent subs","ambiguous convergent subs","num mns","syn mns","nonsyn mns","ambiguous mns","sites","subs"];
    #     # elif out_mode == "clades":
    #     #     headers = ["file","id","chr","start","end","branch","clade","dn","ds","dn/ds","nonsynonymous bl","synonymous bl"];
    # else:
    #     if out_mode == "splits":
    #         #headers = ["file","branch","split1","split2","num subs","num ins","num dels","num mns","sites","subs"];
    #         headers = ["file","branch","split1","split2","num subs","syn subs","nonsyn subs","ambiguous subs","num convergent subs","syn convergent subs","nonsyn convergent subs","ambiguous convergent subs","num mns","syn mns","nonsyn mns","ambiguous mns","sites","subs"];
    #     # elif out_mode == "clades":
    #     #     headers = ["file","branch","clade","dn","ds","dn/ds","nonsynonymous bl","synonymous bl"];
    outfile.write(",".join(headers) + "\n");
    # Compile the headers and write them to the log file

    ####################

    hyphy_files = [f for f in os.listdir(indir) if f.endswith(".anc") ];
    num_files = len(hyphy_files);
    num_files_str = str(num_files);
    num_files_len = len(num_files_str);
    # Read align file names from input directory

    num_unfinished = 0;
    # A count of the number of unfinished hyphy runs as determined by empty output files

    debug = True;

    counter = 0;
    for f in os.listdir(indir):
        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_files_len:
                counter_str = "0" + counter_str;
            print ("> " + hpcore.getDateTime() + " " + counter_str + " / " + num_files_str);
        counter += 1;
        # Track progress

        ####################

        cur_json_file = os.path.join(indir, f);
        if not os.path.isfile(cur_json_file) or not cur_json_file.endswith(".anc"):
            continue;
        if os.stat(cur_json_file).st_size == 0:
            num_unfinished +=1 ;
            continue;
        # Get the current input file

        gene_outfile = os.path.join(outdir, f.replace(".anc", ".csv"));
        # The name of the output file for the current gene

        ####################

        with open(gene_outfile, "w") as goutfile:
            goutfile.write(",".join(headers) + "\n");
            # Open the output file for the current gene and write the headers

            #print(f);
            if features:
                if "-" in f:
                    fid = f.split("-")[0];
                elif "." in f:
                    fid = f.split(".")[0];
                else:
                    fid = f;
                cur_feature = features[fid];
            # Look up transcript/gene info for current gene to save in output, if metadata is provided

            ####################

            with open(cur_json_file) as json_data:
                cur_data = json.load(json_data);
            # Read the HyPhy json file with ancestral reconstructions

            ####################

            cur_tree = cur_data['tree'];
            tinfo, tree, root = tp.treeParse(cur_tree);
            # Read the tree from the json data

            ####################

            #splits = defaultdict(list);
            # The dictionary that stores the splits for each branch
            # key:value = :[ [species from one split], [species from other split] ]
            # NOTE: The first split is the descending clade and the second one is the ancestral clade

            branch_info = {};
            # The dictionary for results for each branch in the current gene

            ancs = {};
            # Dictionary to keep track of ancestors to skip comparisons

            for n in tinfo:
            # For each node/branch in the tree, save the tip species that make up the clade.
                if tinfo[n][2] == 'root':
                    continue;
                if tinfo[n][2] == 'tip':
                    node_label = n;
                else:
                    node_label = tinfo[n][3];
                # For the splits dict, replace the node labels with thos provided by HyPhy

                split1 = tp.getClade(n, tinfo);
                split1 = ";".join(sorted(split1, key=str.casefold))

                #print(cur_clade);
                split2 = [ n2 for n2 in tinfo if tinfo[n2][2] == 'tip' and n2 not in split1 ];
                split2 = ";".join(sorted(split2, key=str.casefold));

                #splits[node_label].append(";".join(split1));
                #splits[node_label].append(";".join(split2));
                # For unrooted trees there is no directionality, so get clades from both sides of each branch.

                anc = tinfo[n][1];
                anc_label = tinfo[anc][3];
                ancs[node_label] = anc_label;

                branch_info[node_label] = { "branch" : node_label, "node type" : tinfo[n][2], "clade" : split1, "split2" : split2, "total subs" : 0 ,"total s" : 0, "total n" : 0, "total a" : 0, "total mns" : 0, "s mns" : 0, "n mns" : 0, "a mns" :0 };
                # Initialize output dict for this node
                # for codon_sub_type1 in codon_sub_types:
                #     for codon_sub_class in codon_sub_classes:
                #         for aa_sub_class in aa_sub_classes:
                #             sub_label = "-".join([codon_sub_type1, codon_sub_class, aa_sub_class]);
                #             branch_info[node_label][sub_label] = 0;
                # Loop over all the possible combinations of sub types and classes to get labels for each count
            ## End split loop

            ####################

            #for node in splits:
            # Loop over every node in the splits dict

                # if features:
                #     if out_mode == "splits":
                #         branch_info[node] = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                #             "branch" : node, "split1" : splits[node][0], "split2" : splits[node][1], "num subs" : 0 ,"syn subs" : 0,"nonsyn subs" : 0, 'ambiguous subs' : 0,
                #             "num convergent subs" : 0, "syn convergent subs" : 0, "nonsyn convergent subs" : 0, "ambiguous convergent subs": 0,
                #             "num mns" : 0, "syn mns" : 0, "nonsyn mns" : 0, 'ambiguous mns' : 0,
                #             "sites" : [], "subs" : [] };
                #     # elif out_mode == "clades":
                #     #     branch_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                #     #         "branch" : node, "clade" : "NA", "dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA", "nonsynonymous bl" : "NA", "synonymous bl" : "NA" };                        
                # else:
                #     if out_mode == "splits": 
                #         branch_info[node] = { "branch" : node, "split1" : splits[node][0], "split2" : splits[node][1], "num subs" : 0 ,"syn subs" : 0, "nonsyn subs" : 0, 'ambiguous subs' : 0,
                #             "num convergent subs" : 0, "syn convergent subs" : 0, "nonsyn convergent subs" : 0, "ambiguous convergent subs": 0,
                #             "num mns" : 0, "syn mns" : 0, "nonsyn mns" : 0, 'ambiguous mns' : 0,
                #             "sites" : [], "subs" : [] };
                #     # elif out_mode == "clades":
                #     #     branch_info = { "branch" : node ,"clades" : "NA", "dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA", "nonsynonymous bl" : "NA", "synonymous bl" : "NA" };   
                # # Initialize the output dictionary for the current branch.
            ## End split header loop

            ####################

            for site in cur_data['substitution_map']:
            # Loop over every site (codon) with substitutions in at least one branch

                #cur_subs = {};
                # Keep track of all substitutions at the current site to pair them up later

                for anc_codon in cur_data['substitution_map'][site]:
                # Loop over every ancestral state at the current site
                    for der_codon in cur_data['substitution_map'][site][anc_codon]:
                    # Loop over every derived state at the current site

                        sub_type = "";
                        # Re-initialize sub_type for each derived codon so the previous "ambiguous" doesn't stick

                        if not all(base in "ATCGatcg" for base in anc_codon) or not all(base in "ATCGatcg" for base in der_codon):
                            sub_type = "a";
                            anc_aa = "a";
                            der_aa = "a";
                        else:
                            anc_aa = hpseq.codonToAA(anc_codon);
                            der_aa = hpseq.codonToAA(der_codon);
                        # Check if either ancestral or derived codon is ambiguous and if not get the corresponding AA

                        #print(anc_codon, anc_aa);
                        #print(der_codon, der_aa);

                        if sub_type == "a":
                            pass;
                        else:
                            if anc_aa == der_aa:
                                sub_type = "s";
                            else:
                                sub_type = "n";
                        # Classify the substitution type based on the AAs

                        #print(sub_type);

                        for key,branch in cur_data['substitution_map'][site][anc_codon][der_codon].items():
                        # Loop over every branch at the site that has this particular substitution

                            branch_info[branch]["total subs"] += 1;
                            branch_info[branch]["total " + sub_type] += 1;

                            num_diffs = sum(1 for a, b in zip(anc_codon, der_codon) if a != b);
                            if num_diffs > 1:
                                branch_info[branch]['total mns'] += 1;
                                branch_info[branch][sub_type + " mns"] += 1;

                            #print(branch);
                            #print(branches[f]);
                            #print(node_label);
                            #print(cur_clades[branch]);

                            #cur_subs[branch] = { 'anc-codon' : anc_codon, 'anc-aa' : anc_aa, 'der-codon' : der_codon, 'der-aa' : der_aa, 'sub-type' : sub_type };
                            # Only get subs from branches that map to the species tree as read in from branch-counts.csv
                        ## End branch loop
                    ## End derived state loop
                ## End ancestral state loop

                ####################

                # if not cur_subs:
                #     continue;

                # done = [];
                # # We only want to count each pair of branches once, so keep track of the ones we've already done

                # for sub1 in cur_subs:
                # # For every branch with a substitution at the current site        
                                    
                #     branch_info[sub1]["total subs"] += 1;
                #     branch_info[sub1]["total " + cur_subs[sub1]['sub-type']] += 1;

                #     num_diffs = sum(1 for a, b in zip(cur_subs[sub1]['anc-codon'], cur_subs[sub1]['der-codon']) if a != b);
                #     if num_diffs > 1:
                #         mns_label = cur_subs[sub1]['sub-type'] + " mns";
                #         branch_info[sub1]['total mns'] += 1;
                #         branch_info[sub1][mns_label] += 1;

                #     for sub2 in cur_subs:
                #     # Compare to all other branches with a substitution at the current site

                #         if sub1 == sub2 or [sub1, sub2] in done or [sub2, sub1] in done or sub1 == ancs[sub2] or sub2 == ancs[sub1]:
                #             continue;
                #         # If the branches are the same, or we've already made this comparison, skip

                #         done.append([sub1, sub2]);
                #         done.append([sub2, sub1]);
                #         # Add both branch orderings to the list of done comparisons

                #         if debug:
                #             print("subs", sub1, sub2);
                #             print(cur_subs[sub1]);
                #             print(cur_subs[sub2]);

                #         ####################

                #         codon_sub_type = cur_subs[sub1]['sub-type'] + "-" + cur_subs[sub2]['sub-type'];
                #         if codon_sub_type not in codon_sub_types:
                #             codon_sub_type = cur_subs[sub2]['sub-type'] + "-" + cur_subs[sub1]['sub-type'];

                #         # if cur_subs[sub1]['sub-type'] == "amb" or cur_subs[sub2]['sub-type'] == "amb":
                #         #     codon_sub_type = "amb";
                #         # elif cur_subs[sub1]['sub-type'] == "syn" or cur_subs[sub2]['sub-type'] == "non-syn":
                #         #     codon_sub_type = "one-each";
                #         # elif cur_subs[sub1]['sub-type'] == "non-syn" or cur_subs[sub2]['sub-type'] == "syn":
                #         #     codon_sub_type = "one-each";
                #         # elif cur_subs[sub1]['sub-type'] == "syn" or cur_subs[sub2]['sub-type'] == "syn":
                #         #     codon_sub_type = "syn";
                #         # elif cur_subs[sub1]['sub-type'] == "non-syn" or cur_subs[sub2]['sub-type'] == "non-syn":
                #         #     codon_sub_type = "non-syn"
                #         # Classify the substitution type based on whether one or both are ambiguous, synonymous, or non-synonymous

                #         ####################

                #         if cur_subs[sub1]['anc-codon'] == cur_subs[sub2]['anc-codon']:
                #             if cur_subs[sub1]['der-codon'] == cur_subs[sub2]['der-codon']:
                #                 codon_sub_class = "par";
                #             else:
                #                 codon_sub_class = "div";
                #         # If the ancestral codons are the same, the subsitutions can be parallel or divergent depending on the derived codons

                #         else:
                #             if cur_subs[sub1]['der-codon'] == cur_subs[sub2]['der-codon']:
                #                 codon_sub_class = "con";
                #             else:
                #                 codon_sub_class = "div";
                #         # If the ancestral codons are different, the substitutions can be convergent or divergent depending on the derived codons

                #         ####################

                #         #print(cur_subs[sub1]['anc-aa'], cur_subs[sub2]['anc-aa'], cur_subs[sub1]['anc-aa'] == cur_subs[sub2]['anc-aa']);
                #         #print(cur_subs[sub1]['der-aa'], cur_subs[sub2]['der-aa'], cur_subs[sub1]['der-aa'] == cur_subs[sub2]['der-aa']);
                #         #print(cur_subs[sub1]['anc-aa'], cur_subs[sub1]['der-aa'], cur_subs[sub1]['anc-aa'] == cur_subs[sub1]['der-aa']);
                #         if codon_sub_type == "s-s":
                #             aa_sub_class = "none";

                #         elif "a" in codon_sub_type:
                #             aa_sub_class = "na";

                #         elif cur_subs[sub1]['anc-aa'] == cur_subs[sub2]['anc-aa']:
                #             if cur_subs[sub1]['der-aa'] == cur_subs[sub2]['der-aa']:
                #                 if cur_subs[sub1]['anc-aa'] == cur_subs[sub1]['der-aa']:
                #                     aa_sub_class = "none";
                #                 else:
                #                     aa_sub_class = "par";
                #             else:
                #                 aa_sub_class = "div";
                #         # If the ancestral AAs are the same, the subsitutions can be parallel or divergent depending on the derived AAs

                #         else:
                #             if cur_subs[sub1]['der-aa'] == cur_subs[sub2]['der-aa']:
                #                 aa_sub_class = "con";
                #             else:
                #                 aa_sub_class = "div";
                #         # If the ancestral AAs are different, the substitutions can be convergent or divergent depending on the derived AAs

                #         ####################

                #         final_sub_label = "-".join([codon_sub_type, codon_sub_class, aa_sub_class]);
                #         # The final substitution label combining whether it is convergent, parallel, or divergent and ambiguous, synonymous, non-synonymous, or one of each

                #         if debug and "s-s-par-none" in final_sub_label:
                #             print("final sub type", final_sub_label);

                #         branch_info[sub1][final_sub_label] += 1;
                #         # Add to the count of substitutions for the current pair of branches for the current substitution type

                #         if debug and "s-s-par-none" in final_sub_label:
                #             print("updated cur_pairs_dict", branch_info[sub1]);
                #             print("------");

                        # if debug and "s-s-par-none" in final_sub_label:
                        #     sys.exit();

                        ####################
                    ## End substitution 2 loop
                ## End substitution 1 loop
            ## End site loop

            ####################

            for node in branch_info:
            # For every node in the current tree, output the counts
                branch_outline = [f] + [ str(branch_info[node][h]) for h in headers if h not in ["file"] ];
                goutfile.write(",".join(branch_outline) + "\n");
                # Format and write the current output line
            ## End node output loop  
        ## Close gene output file (goutfile)
    hpcore.PWS("# ----------------", outfile);
    hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);

############################################################