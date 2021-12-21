#!/usr/bin/env python3
#############################################################################
# Counts pairwise convergent substitutions between all branches in a species
# tree from HyPhy's json output from the standalone ancestral reconstruction
# analysis
#
# Gregg Thomas, Fall 2021
#############################################################################

# pairwise_convergence.py -i full-coding-mg94-local-rooted-anc/ -b TEST/branch-counts.csv -t ../docs/data/trees/full_coding_iqtree_astral.cf.rooted.tree -o test.csv --overwrite

import sys, os, json, argparse, timeit, lib.hpcore as hpcore, lib.hptree as tp, lib.hpseq as hpseq
from collections import defaultdict

############################################################
# Options

prog_start_time = timeit.default_timer();

parser = argparse.ArgumentParser(description="Parse Hyphy json output for pairwise convergence");
parser.add_argument("-i", dest="input", help="Directory containing hyphy json output files from AncestralReconstruction.bf.", default=False);
parser.add_argument("-b", dest="branches", help="The branch-counts.csv file from a branch_rates.py --rooted run.", default=False);
parser.add_argument("-t", dest="tree", help="The species tree file used to run HyPhy", default=False);
parser.add_argument("-o", dest="output", help="An output .csv file.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output file already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# IO options
args = parser.parse_args();

if not args.input or not os.path.isdir(args.input):
    sys.exit(" * Error 1: Please provide a valid input directory (-i).");

if not args.branches or not os.path.isfile(args.branches):
    sys.exit(" * Error 2: Please provide a branch-counts.csv file (-b).");

if not args.tree or not os.path.isfile(args.tree):
    sys.exit(" * Error 3: Please provide a rooted species tree in Newick format (-t).");

if not args.output:
    sys.exit(" * Error 4: Please provide the name of an output file (-o).")

if os.path.isfile(args.output) and not args.overwrite:
    sys.exit( " * Error 5: Output file (-o) already exists! Explicity specify --overwrite to overwrite it.");

pad = 25;
debug = False;

with open(args.output, "w") as outfile:
    hpcore.runTime("# Pairwise convergent substitution counting", outfile);
    hpcore.PWS("# IO OPTIONS", outfile);
    hpcore.PWS(hpcore.spacedOut("# Input directory:", pad) + args.input, outfile);
    hpcore.PWS(hpcore.spacedOut("# Branch counts file:", pad) + args.branches, outfile);
    hpcore.PWS(hpcore.spacedOut("# Tree file:", pad) + args.tree, outfile);
    hpcore.PWS(hpcore.spacedOut("# Output file:", pad) + args.output, outfile);
    if args.overwrite:
        hpcore.PWS(hpcore.spacedOut("# --overwrite set:", pad) + "Overwriting previous output file.", outfile);
    hpcore.PWS("# ----------------", outfile);

    step_start_time = hpcore.report_step("", "", "", prog_start_time, start=True);
    # Initialize the step headers

    ####################

    step = "Reading species tree";
    step_start_time = hpcore.report_step(step, False, "In progress...", prog_start_time);
    # Status update

    try:
        tree_string = open(args.tree, "r").read().strip();
        tinfo, labeled_tree, root = tp.treeParse(tree_string)
    except:
        hpcore.errorOut(" * Error 6: reading species tree!");

    step_start_time = hpcore.report_step(step, step_start_time, "SUCCESS", prog_start_time);
    #print(labeled_tree);
    # Status update

    ####################

    headers = ["node1", "node2"];
    # Headers for the output file

    sum_headers = ["s-s-con-none", "s-s-par-none", "s-s-div-none"]

    codon_sub_types = ["s-n", "s-a", "n-n", "n-a", "a-a"];
    codon_sub_classes = ["con", "par", "div"];
    aa_sub_classes = ["con", "par", "div", "none", "na"];
    # Various substitution types and classes to count among pairs of branches

    for codon_sub_type in codon_sub_types:
        for codon_sub_class in codon_sub_classes:
            for aa_sub_class in aa_sub_classes:
                sub_label = "-".join([codon_sub_type, codon_sub_class, aa_sub_class]);
                sum_headers.append(sub_label);
    # Generate the headers for each substitution label   

    headers += sum_headers;
    headers += ["clade1", "clade2"];

    ####################

    step = "Pairing branches";
    step_start_time = hpcore.report_step(step, False, "In progress...", prog_start_time);
    # Status update

    pairs_dict = {};

    ancs = {};
    # Dictionary to keep track of ancestors to skip comparisons

    for n1 in tinfo:
        if n1 == root:
            continue;

        for n2 in tinfo:
            if n2 == n1 or n2 == root:
                continue;

            if tinfo[n1][1] == n2 or tinfo[n2][1] == n1:
                continue;
            # Skip if one of the branches is a direct descendant of the other

            clade1 = sorted(tp.getClade(n1, tinfo), key=str.casefold);
            clade2 = sorted(tp.getClade(n2, tinfo), key=str.casefold);

            # if set(clade1) == set(['Arvicanthus_niloticus_MNHN1999201','Arvicanthus_neumanni_FMNH158037','Lemniscomys_striatus_TCD4711']):
            #     if set(clade2) == set(['Papagomys_armandvillei_WAMM32592']):
            #         print(clade1);
            #         print(clade2);
            #         print(tuple(sorted([";".join(clade1), ";".join(clade2)], key=str.casefold)));
            #         sys.exit();
            # if set(clade2) == set(['Arvicanthus_niloticus_MNHN1999201','Arvicanthus_neumanni_FMNH158037','Lemniscomys_striatus_TCD4711']):
            #     if set(clade1) == set(['Papagomys_armandvillei_WAMM32592']):
            #         print(clade2);
            #         print(clade1);
            #         print(tuple(sorted([";".join(clade1), ";".join(clade2)], key=str.casefold)));
            #         sys.exit();
            pair = sorted([";".join(clade1), ";".join(clade2)], key=str.casefold);

            pairs_dict[tuple(pair)] = { "node1" : n1, "node2" : n2 };

            for h in sum_headers:
                pairs_dict[tuple(pair)][h] = 0.0;
            # Generate the labels for each substitution type for each pair of branches   

    step_start_time = hpcore.report_step(step, step_start_time, "SUCCESS: " + str(len(pairs_dict)) + " pairs", prog_start_time);
    # Status update

    # print(pairs_dict[('Arvicanthus_niloticus_MNHN1999201;Arvicanthus_neumanni_FMNH158037;Lemniscomys_striatus_TCD4711', 'Papagomys_armandvillei_WAMM32592')]);
    # sys.exit();

    ####################            

    step = "Reading branch presence file";
    step_start_time = hpcore.report_step(step, False, "In progress...", prog_start_time);
    # Status update

    branches = defaultdict(dict);
    for line in open(args.branches):
        line = line.strip().split(",");
        f = line[0].replace(".csv", ".anc");
        branches[f][line[5]] = line[2];#     .append({ line[5], line[2] ]);
    # Gene : { gene tree clade : species tree clade } 

    # for gene in branches:
    #     print(gene);
    #     print(branches[gene]);
    #     sys.exit();

    # print(branches);
    # sys.exit();

    step_start_time = hpcore.report_step(step, step_start_time, "SUCCESS: " + str(len(branches)) + " genes read", prog_start_time);
    # Status update

    ####################

    step = "Counting HyPhy files";
    step_start_time = hpcore.report_step(step, False, "In progress...", prog_start_time);
    # Status update

    hyphy_files = [f for f in os.listdir(args.input) if f.endswith(".anc") ];
    num_files = len(hyphy_files);
    num_files_str = str(num_files);
    num_files_len = len(num_files_str);
    # Read align file names from input directory

    step_start_time = hpcore.report_step(step, step_start_time, "SUCCESS: " + num_files_str + " files counted", prog_start_time);
    # Status update

    ####################

    step = "Counting substitutions";
    #step_start_time = hpcore.report_step(step, False, "Processed 0 / " + num_files_str + " loci...", prog_start_time, full_update=True);
    # Status update

    num_unfinished = 0;
    # A count of the number of unfinished hyphy runs as determined by empty output files

    counter = 0;
    for f in hyphy_files:
    # Count substitutions for every file

        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_files_len:
                counter_str = "0" + counter_str;
            #print ("> " + hpcore.getDateTime() + " " + counter_str + " / " + num_files_str);
            step_start_time = hpcore.report_step(step, step_start_time, "Processed " + counter_str + " / " + num_files_str + " loci...", prog_start_time, full_update=True);
            # Status update
        counter += 1;
        # Track progress

        ####################

        cur_json_file = os.path.join(args.input, f);
        if not os.path.isfile(cur_json_file) or not cur_json_file.endswith(".anc"):
            continue;
        if os.stat(cur_json_file).st_size == 0:
            num_unfinished +=1 ;
            continue;
        # Get the current input file

        ####################

        with open(cur_json_file) as json_data:
            cur_data = json.load(json_data);
        # Read the HyPhy json file with ancestral reconstructions

        ####################

        cur_tree = cur_data['tree'];
        tinfo, tree, root = tp.treeParse(cur_tree);
        # Read the tree from the json data

        ####################

        cur_clades = {};
        # A dict of clades present in the current gene tree
        # node label : clade

        for n in tinfo:
        # For each node/branch in the tree, save the tip species that make up the clade.
            if tinfo[n][2] == 'root':
                continue;
            if tinfo[n][2] == 'tip':
                node_label = n;
            else:
                node_label = tinfo[n][3];
            # For the splits dict, replace the node labels with thos provided by HyPhy

            clade = sorted(tp.getClade(n, tinfo), key=str.casefold);

            cur_clades[node_label] = ";".join(clade);
        ## End clade loop

        ####################

        for site in cur_data['substitution_map']:
        # Loop over every site (codon) with substitutions in at least one branch

            cur_subs = {};
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

                        #print(branch);
                        #print(branches[f]);
                        #print(node_label);
                        #print(cur_clades[branch]);
                        if cur_clades[branch] in branches[f]:
                            cur_subs[branch] = { 'anc-codon' : anc_codon, 'anc-aa' : anc_aa, 'der-codon' : der_codon, 'der-aa' : der_aa, 'sub-type' : sub_type, 'st-clade' : branches[f][cur_clades[branch]] };
                        # Only get subs from branches that map to the species tree as read in from branch-counts.csv
                    ## End branch loop
                ## End derived state loop
            ## End ancestral state loop
        
            ####################

            if len(cur_subs) == 1:
                continue;
            # If there is only one substitution at this site, there's no possibility for convergence, so skip

            else:
            # Otherwise, count pairs of substitutions

                #print(f);
                #print("cursubs", cur_subs, len(cur_subs));
                
                done = [];
                # We only want to count each pair of branches once, so keep track of the ones we've already done

                for sub1 in cur_subs:
                # For every branch with a substitution at the current site                    

                    for sub2 in cur_subs:
                    # Compare to all other branches with a substitution at the current site

                        #print("subs", sub1, sub2);
                        if sub1 == sub2 or [sub1, sub2] in done or [sub2, sub1] in done:
                            continue;
                        # If the branches are the same, or we've already made this comparison, skip

                        done.append([sub1, sub2]);
                        done.append([sub2, sub1]);
                        # Add both branch orderings to the list of done comparisons

                        if debug:
                            print("subs", sub1, sub2);
                            print(cur_subs[sub1]);
                            print(cur_subs[sub2]);

                        ####################

                        st_clade1 = cur_subs[sub1]['st-clade'];
                        st_clade2 = cur_subs[sub2]['st-clade'];
                        # Retrieve the species tree clades that the current branches map to

                        cur_pair = tuple(sorted([st_clade1, st_clade2], key=str.casefold));
                        # Make the pairing of the mapped branches to lookup in the pairs_dict

                        if debug:
                            print("cur_pair", cur_pair);
                            print("cur_pairs_dict", pairs_dict[cur_pair]);

                        if cur_pair not in pairs_dict:
                            if debug:
                                print("SKIPPING PAIR");
                                print("------");
                            continue;

                        ####################

                        codon_sub_type = cur_subs[sub1]['sub-type'] + "-" + cur_subs[sub2]['sub-type'];
                        if codon_sub_type not in codon_sub_types:
                            codon_sub_type = cur_subs[sub2]['sub-type'] + "-" + cur_subs[sub1]['sub-type'];

                        # if cur_subs[sub1]['sub-type'] == "amb" or cur_subs[sub2]['sub-type'] == "amb":
                        #     codon_sub_type = "amb";
                        # elif cur_subs[sub1]['sub-type'] == "syn" or cur_subs[sub2]['sub-type'] == "non-syn":
                        #     codon_sub_type = "one-each";
                        # elif cur_subs[sub1]['sub-type'] == "non-syn" or cur_subs[sub2]['sub-type'] == "syn":
                        #     codon_sub_type = "one-each";
                        # elif cur_subs[sub1]['sub-type'] == "syn" or cur_subs[sub2]['sub-type'] == "syn":
                        #     codon_sub_type = "syn";
                        # elif cur_subs[sub1]['sub-type'] == "non-syn" or cur_subs[sub2]['sub-type'] == "non-syn":
                        #     codon_sub_type = "non-syn"
                        # Classify the substitution type based on whether one or both are ambiguous, synonymous, or non-synonymous

                        ####################

                        if cur_subs[sub1]['anc-codon'] == cur_subs[sub2]['anc-codon']:
                            if cur_subs[sub1]['der-codon'] == cur_subs[sub2]['der-codon']:
                                codon_sub_class = "par";
                            else:
                                codon_sub_class = "div";
                        # If the ancestral codons are the same, the subsitutions can be parallel or divergent depending on the derived codons

                        else:
                            if cur_subs[sub1]['der-codon'] == cur_subs[sub2]['der-codon']:
                                codon_sub_class = "con";
                            else:
                                codon_sub_class = "div";
                        # If the ancestral codons are different, the substitutions can be convergent or divergent depending on the derived codons

                        ####################

                        if codon_sub_type == "s-s":
                            aa_sub_class = "none";

                        elif "a" in codon_sub_type:
                            aa_sub_class = "na";

                        elif cur_subs[sub1]['anc-aa'] == cur_subs[sub2]['anc-aa']:
                            if cur_subs[sub1]['der-aa'] == cur_subs[sub2]['der-aa']:
                                if cur_subs[sub1]['anc-aa'] == cur_subs[sub1]['der-aa']:
                                    aa_sub_class = "none";
                                else:
                                    aa_sub_class = "par";
                            else:
                                aa_sub_class = "div";
                        # If the ancestral AAs are the same, the subsitutions can be parallel or divergent depending on the derived AAs

                        else:
                            if cur_subs[sub1]['der-aa'] == cur_subs[sub2]['der-aa']:
                                aa_sub_class = "con";
                            else:
                                aa_sub_class = "div";
                        # If the ancestral AAs are different, the substitutions can be convergent or divergent depending on the derived AAs

                        ####################

                        final_sub_label = "-".join([codon_sub_type, codon_sub_class, aa_sub_class]);
                        # The final substitution label combining whether it is convergent, parallel, or divergent and ambiguous, synonymous, non-synonymous, or one of each
                        
                        if debug and "s-s-par-none" in final_sub_label:
                            print("final sub type", final_sub_label);

                        pairs_dict[cur_pair][final_sub_label] += 1;
                        # Add to the count of substitutions for the current pair of branches for the current substitution type

                        if debug and "s-s-par-none" in final_sub_label:
                            print("updated cur_pairs_dict", pairs_dict[cur_pair]);
                            print("------");

                        ####################
                    ## End substitution 2 loop
                ## End substitution 1 loop
                #sys.exit();
            ## End sub counting block
            #print("-----------------");
        ## End site loop
    ## End file loop        

    step_start_time = hpcore.report_step(step, step_start_time, "SUCCESS", prog_start_time, full_update=True);
    # Status update

    ####################

    step = "Counting HyPhy files";
    step_start_time = hpcore.report_step(step, False, "In progress...", prog_start_time);
    # Status update

    outfile.write(",".join(headers) + "\n");

    for pair in pairs_dict:
        outline = [];
        for h in headers[:-2]:
            outline.append(str(pairs_dict[pair][h]));

        outline.append(pair[0]);
        outline.append(pair[1]);

        outfile.write(",".join(outline) + "\n");
    ## End node output loop  

    step_start_time = hpcore.report_step(step, step_start_time, "SUCCESS", prog_start_time);
    # Status update

    ####################