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

    if features:
        if out_mode == "splits":
            #headers = ["file","id","chr","start","end","branch","split1","split2","num subs","num ins","num dels","num mns","sites","subs"];
            headers = ["file","id","chr","start","end","branch","split1","split2","num subs","syn subs","nonsyn subs","ambiguous subs","num convergent subs","syn convergent subs","nonsyn convergent subs","ambiguous convergent subs","num mns","syn mns","nonsyn mns","ambiguous mns","sites","subs"];
        # elif out_mode == "clades":
        #     headers = ["file","id","chr","start","end","branch","clade","dn","ds","dn/ds","nonsynonymous bl","synonymous bl"];
    else:
        if out_mode == "splits":
            #headers = ["file","branch","split1","split2","num subs","num ins","num dels","num mns","sites","subs"];
            headers = ["file","branch","split1","split2","num subs","syn subs","nonsyn subs","ambiguous subs","num convergent subs","syn convergent subs","nonsyn convergent subs","ambiguous convergent subs","num mns","syn mns","nonsyn mns","ambiguous mns","sites","subs"];
        # elif out_mode == "clades":
        #     headers = ["file","branch","clade","dn","ds","dn/ds","nonsynonymous bl","synonymous bl"];
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

            splits = defaultdict(list);
            # The dictionary that stores the splits for each branch
            # key:value = :[ [species from one split], [species from other split] ]
            # NOTE: The first split is the descending clade and the second one is the ancestral clade

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
                split1.sort();
                #print(cur_clade);
                split2 = [ n2 for n2 in tinfo if tinfo[n2][2] == 'tip' and n2 not in split1 ];
                split2.sort();

                splits[node_label].append(";".join(split1));
                splits[node_label].append(";".join(split2));
                # For unrooted trees there is no directionality, so get clades from both sides of each branch.
            ## End split loop

            ####################

            branch_info = {};
            # The dictionary for results for each branch in the current gene

            for node in splits:
            # Loop over every node in the splits dict
                if features:
                    if out_mode == "splits":
                        branch_info[node] = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                            "branch" : node, "split1" : splits[node][0], "split2" : splits[node][1], "num subs" : 0 ,"syn subs" : 0,"nonsyn subs" : 0, 'ambiguous subs' : 0,
                            "num convergent subs" : 0, "syn convergent subs" : 0, "nonsyn convergent subs" : 0, "ambiguous convergent subs": 0,
                            "num mns" : 0, "syn mns" : 0, "nonsyn mns" : 0, 'ambiguous mns' : 0,
                            "sites" : [], "subs" : [] };
                    # elif out_mode == "clades":
                    #     branch_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                    #         "branch" : node, "clade" : "NA", "dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA", "nonsynonymous bl" : "NA", "synonymous bl" : "NA" };                        
                else:
                    if out_mode == "splits": 
                        branch_info[node] = { "branch" : node, "split1" : splits[node][0], "split2" : splits[node][1], "num subs" : 0 ,"syn subs" : 0, "nonsyn subs" : 0, 'ambiguous subs' : 0,
                            "num convergent subs" : 0, "syn convergent subs" : 0, "nonsyn convergent subs" : 0, "ambiguous convergent subs": 0,
                            "num mns" : 0, "syn mns" : 0, "nonsyn mns" : 0, 'ambiguous mns' : 0,
                            "sites" : [], "subs" : [] };
                    # elif out_mode == "clades":
                    #     branch_info = { "branch" : node ,"clades" : "NA", "dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA", "nonsynonymous bl" : "NA", "synonymous bl" : "NA" };   
                # Initialize the output dictionary for the current branch.
            ## End split header loop

            ####################

            for site in cur_data['substitution_map']:
            # Loop over every site (codon) with substitutions in at least one branch

                for anc_codon in cur_data['substitution_map'][site]:
                # Loop over every ancestral state at the current site

                    for der_codon in cur_data['substitution_map'][site][anc_codon]:
                    # Loop over every derived state at the current site

                        sub_str = anc_codon + ":" + der_codon;
                        # Get the string to save in the subs list by combining the two codons

                        if not all(base in "ATCGatcg" for base in anc_codon) or not all(base in "ATCGatcg" for base in der_codon):
                            sub_type = "ambiguous";
                        # If either codon contains ambiguous bases, the current substitution type is also ambiguous
                        # TODO: Could actually do a better check: if all possible codons result in synonymous subs then sub type is synonymous. Same for non-synonymous

                        elif hpseq.codonToAA(anc_codon) == hpseq.codonToAA(der_codon):
                            sub_type = "syn";
                        # If the change results in the same amino acid, type is synonymous
                        else:
                            sub_type = "nonsyn";
                        # If the change results in a different amino acid, type is nonsynonymous
                        ## End sub type block

                        if len(cur_data['substitution_map'][site][anc_codon][der_codon]) > 1:
                            convergent = True;
                        else:
                            convergent = False;

                        for key,branch in cur_data['substitution_map'][site][anc_codon][der_codon].items():
                        # Loop over every branch at the site that has this particular substitution

                            branch_info[branch]['num subs'] += 1;
                            branch_info[branch][sub_type + ' subs'] += 1;
                            # Increment counts for substitutions

                            if (anc_codon[0] != der_codon[0] and anc_codon[1] != der_codon[1]) or (anc_codon[1] != der_codon[1] and anc_codon[2] != der_codon[2]):
                                branch_info[branch]['num mns'] += 1;
                                branch_info[branch][sub_type + ' mns'] += 1;
                            # If the substitution contains 2 adjacent changes, count it as a multi-nucleotide substitution

                            if convergent:
                                branch_info[branch]['num convergent subs'] += 1;
                                branch_info[branch][sub_type + ' convergent subs'] += 1;

                            branch_info[branch]['sites'].append(site);
                            branch_info[branch]['subs'].append(sub_str);
                            # Append the site and the substitutions to the lists for this branch
                        ## End branch loop
                    ## End derived state loop
                ## End ancestral state loop
            ## End site loop

            ####################

            for node in branch_info:
            # For every node in the current tree, output the counts
                if out_mode == "splits":
                    branch_outline = [f] + [ str(branch_info[node][h]) for h in headers if h not in ["file"] ];
                    #outfile.write(",".join(branch_outline) + "\n");
                    goutfile.write(",".join(branch_outline) + "\n");
                # Format and write the current output line
                      
                # if out_mode == "clades":
                #     for clade in splits[node]:
                #         branch_info["clade"] = clade;
                #         branch_outline = [f] + [ branch_info[h] for h in headers if h not in ["file"] ];
                #         outfile.write(",".join(branch_outline) + "\n");
                #         goutfile.write(",".join(branch_outline) + "\n");
                # # Ouput rate estimates for both clades of the current branch.

            ## End node output loop  
        ## Close gene output file (goutfile)
    hpcore.PWS("# ----------------", outfile);
    hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);