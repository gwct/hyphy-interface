#############################################################################
# Function for determining whether a branch exists in a gene tree
#############################################################################

import sys
import os
import lib.tree as TREE

#############################################################################

def branchSum(branch_sum_item):
# Retrieve rates/counts for each branch in the species tree for each gene, if that branch exists
# in the gene

    cur_branch_dict, cur_branch, rates_dir, csv_files, filter_files, subset_files, rooted = branch_sum_item;

    # print(cur_branch);
    # print(cur_branch_dict);
    #core.PWS("# " + core.getDateTime() + " Starting branch " + cur_branch);

    if cur_branch_dict["node.type"] == "ROOT":
            return [cur_branch, cur_branch_dict];
    # Skip the root node because it doesn't have an associated branch

    branch_list = cur_branch.split(";");
    # Make a list of the species descendant from the current branch

    num_genes = len(csv_files);
    num_genes_str = str(num_genes);
    # The number of branches and genes, for progress updates, str here to minimize calls to str

    gene_counter = 1;
    # The gene counter for progress updates, the branch counter is also converted to string here to 
    # minimize calls to str

    for f in csv_files:
    # For the current branch, go through every gene and get the maximal subset clade

        if f in filter_files:
            continue;
        # Skip the genes filtered based on dS

        if subset_files and f not in subset_files:
            continue;
        # Skip the genes not in the subset of genes we're interested in

        # if gene_counter == 1 or gene_counter % 100 == 0:
        #     print("# Branch : " + str(cur_branch_num) + " / " + num_branches + " -> gene: " + str(gene_counter) + " / " + num_genes_str);
        gene_counter += 1;
        # Progress updates

        # if gene_counter > 10:
        #     continue;

        infile = os.path.join(rates_dir, f);
        # Compile path to current file

        max_split1 = [];
        max_split2 = [];
        max_line = [];
        # Variables to store the information for the maximal subset clade

        first = True;
        for line in open(infile):
            if first:
                first = False;
                continue;
            # Skip the header line in the file

            line = line.strip().split(",");
            split1 = line[2].split(";");
            split2 = line[3].split(";");
            splits = [split1, split2];
            # Get the splits from the current branch in the gene's gene tree

            for s in range(len(splits)):
            # For each split, we want to check whether it is a subset of the current branch

                if set(splits[s]) == set(branch_list) or set(splits[s]) <= set(branch_list):
                # Check if the current split is an equivalent set or subset of the current branch

                    if len(splits[s]) > len(max_split1):
                    # Check if the current split that is a subset of the branch contains more species
                    # than the current maximal subset. If so, it becomes the maximal subset.

                        max_split1 = splits[s];
                        max_line = line;
                        if s == 0:
                            max_split2 = splits[1];
                        else:
                            max_split2 = splits[0];
                        # Assign the maximal subset to max_split1 and the other split to max_split2

        if max_split1 == [] or any(spec in max_split2 for spec in branch_list):
            cur_branch_dict['num.genes.no.clade'] += 1;
        # If no species from the current branch are found then this clade doesn't exist in this gene OR
        # If species from the branch are present in the second split (not the maximal subset), then this is a case of
        # discordance and the clade truly doesn't exist in this gene as a monophyly. Do not increment sums.

        elif len(branch_list) == len(max_split1):
            if not rooted:
                cur_ES = float(max_line[4]);
                cur_EN = float(max_line[5]);
                cur_S = float(max_line[6]);
                cur_N = float(max_line[7]);
                # Parse the rates from the line corresponding to the maximal split

                #if cur_ds_bl == 1e-10:
                #    cur_branch_dict['num.genes.no.ds'] += 1;
                #    cur_ds, cur_ds_bl, cur_dnds = 0, 0, 0;
                #elif cur_dn_bl == 1e-10:
                #    cur_branch_dict['num.genes.no.dn'] += 1;
                #    cur_dn, cur_dn_bl, cur_dnds = 0, 0, 0;

                cur_branch_dict['num.genes.full'] += 1;
                cur_branch_dict['ES.sum'] += cur_ES;
                cur_branch_dict['EN.sum'] += cur_EN;
                cur_branch_dict['S.sum'] += cur_S;
                cur_branch_dict['N.sum'] += cur_N;
            
            else:
                cur_cS = float(max_line[8]);
                cur_cN = float(max_line[9]);
                cur_cA = float(max_line[10]);
                cur_mS = float(max_line[11]);
                cur_mN = float(max_line[12]);
                cur_mA = float(max_line[13]);

                # print(cur_mS);
                # print(cur_branch_dict['mS']);
                # print("------")

                cur_branch_dict['num.genes.full'] += 1;
                cur_branch_dict['cS'] += cur_cS;
                cur_branch_dict['cN'] += cur_cN;
                cur_branch_dict['cA'] += cur_cA;
                cur_branch_dict['mS'] += cur_mS;
                cur_branch_dict['mN'] += cur_mN;
                cur_branch_dict['mA'] += cur_mA;
        # If the maximal subset of the current branch is equal to the current branch, increment sums for each new column.

        else:
            if not rooted:
                cur_ES = float(max_line[4]);
                cur_EN = float(max_line[5]);
                cur_S = float(max_line[6]);
                cur_N = float(max_line[7]);
                # Parse the rates from the line corresponding to the maximal split

                if len(branch_list) == 1:
                    print(cur_branch);
                    print(max_split1);
                    print(len(max_split1));
                    print(max_split2);
                    print(len(max_split2));
                    sys.exit();
                # This shouldn't happen ... it would mean a tip was found in the opposite split of the maximal clade (which should just be the tip)

                #if cur_ds_bl == 1e-10:
                #    cur_branch_dict['num.genes.no.ds'] += 1;
                #    cur_ds, cur_ds_bl, cur_dnds = 0, 0, 0;
                #elif cur_dn_bl == 1e-10:
                #    cur_branch_dict['num.genes.no.dn'] += 1;
                #    cur_dn, cur_dn_bl, cur_dnds = 0, 0, 0;
                # If there are no substitutions of a given type, the branch length will be 0, which I think Hyphy represents as 1e-10
                # Count these here and set appropriate rates to 0

                cur_branch_dict['num.genes.partial'] += 1;
                cur_branch_dict['ES.sum'] += cur_ES;
                cur_branch_dict['EN.sum'] += cur_EN;
                cur_branch_dict['S.sum'] += cur_S;
                cur_branch_dict['N.sum'] += cur_N;
            else:
                cur_cS = float(max_line[8]);
                cur_cN = float(max_line[9]);
                cur_cA = float(max_line[10]);
                cur_mS = float(max_line[11]);
                cur_mN = float(max_line[12]);
                cur_mA = float(max_line[13]);

                cur_branch_dict['num.genes.partial'] += 1;
                cur_branch_dict['cS'] += cur_cS;
                cur_branch_dict['cN'] += cur_cN;
                cur_branch_dict['cA'] += cur_cA;
                cur_branch_dict['mS'] += cur_mS;
                cur_branch_dict['mN'] += cur_mN;
                cur_branch_dict['mA'] += cur_mA;

        # If there are no species from the current branch in the second split (not the maximal subset), then this 
        # is a case of missing data and we can still count the branch as existing in this gene tree. Increment
        # the sums for each new column.

    #core.PWS("# " + core.getDateTime() + " Finishing branch " + cur_branch);
    #print(cur_branch_dict)

    return [cur_branch, cur_branch_dict];

#############################################################################

def branchSumNEW(branch_sum_item):
# Retrieve rates/counts for each branch in the species tree for each gene, if that branch exists
# in the gene

    cur_files, rates_dir, tree_dict, model, corrected_alpha = branch_sum_item;
    # Unpack the passed items

    branch_outlines = [];

    ##########################

    headers = ["node.label", "clade", "node.type", "orig.node.label"];
    if model == 'slac':
        count_headers = ["num.genes.full", "num.genes.partial", "num.genes.descendant.counted", "num.genes.discordant", "num.genes.missing"];
        sum_headers = ["ES", "EN", "S", "N"];   

    elif model == 'absrel':
        count_headers = ["num.genes.full", "num.genes.partial", "num.genes.descendant.counted", "num.genes.discordant", "num.genes.missing"];
        sum_headers = ["num.ps.genes", "ps.gene.list"];

    elif model == 'ancrecon':
        count_headers = ["num.genes.full", "num.genes.partial", "num.genes.descendant.counted", "num.genes.discordant", "num.genes.missing"];
        sum_headers = ["total subs", "total mns", "S", "N", "A", "mS", "mN", "mA"];

    headers += count_headers + sum_headers;
    # The headers/keys for each node in the species tree, depending on whether this is from a rooted ancestral sequences run or SLAC

    ##########################

    cur_dict = {};
    # The dictionary to keep track of counts for the current set of genes

    tips = [];
    internals = [];
    # Have to parse out the species tree nodes to sort them for a post-order traversal

    for node in tree_dict:
    # Setup the dictionary for each node in the species tree

        cur_dict[node] = {};
        # Add the current node to the tracker dict

        cur_clade = sorted(TREE.getClade(node, tree_dict), key=str.casefold);
        
        cur_dict[node]["clade"] = set(cur_clade);
        # Get the clade descending from the current node

        cur_dict[node]["node.label"] = node;
        cur_dict[node]["node.type"] = tree_dict[node][2];
        cur_dict[node]["orig.node.label"] = tree_dict[node][3];
        # Get the other info for the current node

        for h in count_headers + sum_headers:
            cur_dict[node][h] = 0.0;
        if model == 'absrel':
            cur_dict[node]['ps.gene.list'] = [];
        # Initialize all the counts for the current node

        if tree_dict[node][2] == "tip":
            tips.append(node);
        else:
            internals.append(int(node.replace("<", "").replace(">", "")));
        # Add the node to the appropriate list, and remove the extra chars from the internals so they can be sorted
    ## End species tree node loop

    internals.sort();
    # Sort the internal nodes for a post order traversal

    sorted_nodes = tips + internals;
    # Combine the tips and sorted internals

    ##########################

    for f in cur_files:
    # Loop over every file in the current chunk

        # f = "ENSMUSP00000051922-mafft-cds.filter.csv";
        # print(f);

        infile = os.path.join(rates_dir, f);
        # Compile path to current file        

        gt_splits = {};
        # Dictionary to store info for every node in the current gene tree

        first = True;
        for line in open(infile):
        # Read the current file line by line
            if first:
                #print(line);
                first = False;
                continue;
            # Skip the header line in the file

            line = line.strip().split(",");
            split1 = line[2].split(";");
            split2 = line[3].split(";");
            #splits = [split1, split2];
            # Get the splits from the current branch in the gene's gene tree

            gt_splits[line[1]] = { 'split1' : set(split1), 'split2' : set(split2) };
            # Store the splits in the gene tree dict

            if model == 'slac':
                gt_splits[line[1]]['ES'] = float(line[4]);
                gt_splits[line[1]]['EN'] = float(line[5]);
                gt_splits[line[1]]['S'] = float(line[6]);
                gt_splits[line[1]]['N'] = float(line[7]);

            elif model == 'absrel':
                gt_splits[line[1]]['pval'] = float(line[5]);

            elif model == 'ancrecon':
                sub_index = 5
                for sub_label in sum_headers:
                    gt_splits[line[1]][sub_label] = float(line[sub_index]);
                    sub_index += 1;
            # Get the rest of the information for the current branch depending on the type of run
        ## End current file loop

        ##########################

        gt_done = [];
        # A list of gene tree branches that have already been counted (so none are counted twice)

        for node in sorted_nodes:
        # Now go through every node in the species tree in a post-order traversal

            if node in internals:
                node_str = "<" + str(node) + ">";
            else:
                node_str = node;

            max_split1 = set();
            max_split2 = set();
            max_branch = "";
            # Variables to store the information for the maximal subset clade

            for branch in gt_splits:
            # For each branch in the gene tree, we want to check whether either split is a subset of the current branch

                splits = [gt_splits[branch]['split1'], gt_splits[branch]['split2']];
                # Get the splits from the current branch in the gene's gene tree

                for s in range(len(splits)):
                # For each split, we want to check whether it is a subset of the current branch

                    if splits[s] == cur_dict[node_str]["clade"] or splits[s] <= cur_dict[node_str]["clade"]:
                    # Check if the current split is an equivalent set or subset of the current branch

                        if len(splits[s]) > len(max_split1):
                        # Check if the current split that is a subset of the branch contains more species
                        # than the current maximal subset. If so, it becomes the maximal subset.

                            max_split1 = splits[s];
                            max_branch = branch;
                            if s == 0:
                                max_split2 = splits[1];
                            else:
                                max_split2 = splits[0];
                            # Assign the maximal subset to max_split1 and the other split to max_split2
                        ## End max split check block
                    ## End split subset check block
                ## End split loop
            ## End gene tree branch loop

            ##########################

            #if node_str == "Lophiomys_imhausi_UM5152":# "Phloeomys_cumingi_FMNH198472":
            # print(node_str);
            # print(cur_dict[node_str]["clade"]);
            # print(max_split1);
            # print(max_split2);
            # print(max_branch);
            
            clade_found = False;

            if max_branch in gt_done:
                #print("descendant counted");
                cur_dict[node_str]['num.genes.descendant.counted'] += 1;
            # Check that the chosen gene tree branch hasn't been counted before as one of its descendants. If so,
            # skip it.

            elif not max_split1:
                #print("missing");
                cur_dict[node_str]['num.genes.missing'] += 1;
            # If no species from the current branch are found then this clade doesn't exist in this gene

            elif cur_dict[node_str]["clade"] & max_split2:
                #print("discordant");
                cur_dict[node_str]['num.genes.discordant'] += 1;
            # If there exists a non-empty intersect between the current clade and the species not in the maximal sub-clade then
            # species from the branch are present in the second split, and this is a case of
            # discordance and the clade truly doesn't exist in this gene as a monophyly. Do not increment sums.

            elif cur_dict[node_str]["clade"] == max_split1:
                #print("full");
                cur_dict[node_str]['num.genes.full'] += 1;
                if model in ['slac', 'ancrecon']:
                    for h in sum_headers:
                        cur_dict[node_str][h] += gt_splits[max_branch][h];
                elif model == 'absrel':
                    if gt_splits[line[1]]['pval'] < corrected_alpha:
                        cur_dict[node_str]['num.ps.genes'] += 1;
                        cur_dict[node_str]['ps.gene.list'].append(f);
                gt_done.append(max_branch);
                clade_found = "full";
            # If the maximal subset of the current branch is exactly equal to the current branch, increment sums for each new column.

            else:
                #print("partial");
                cur_dict[node_str]['num.genes.partial'] += 1;
                if model in ['slac', 'ancrecon']:
                    for h in sum_headers:
                        cur_dict[node_str][h] += gt_splits[max_branch][h];
                elif model == 'absrel':
                    if gt_splits[line[1]]['pval'] < corrected_alpha:
                        cur_dict[node_str]['num.ps.genes'] += 1;
                        cur_dict[node_str]['ps.gene.list'].append(f);
                gt_done.append(max_branch);
                clade_found = "partial";
            # If there are no species from the current branch in the second split (not the maximal subset), then this 
            # is a case of missing data and we can still count the branch as existing in this gene tree. Increment
            # the sums for each new column.

            if clade_found:
                cur_clade = ";".join(sorted(list(cur_dict[node_str]["clade"]), key=str.casefold));
                num_clade_spec = str(len(cur_dict[node_str]["clade"]));

                found_clade = ";".join(sorted(list(max_split1), key=str.casefold));
                num_found_spec = str(len(max_split1));

                branch_outlines.append([f, node_str, cur_clade, num_clade_spec, clade_found, found_clade, num_found_spec]);
        ## End species tree node loop
    ## End files loop

    #core.PWS("# " + core.getDateTime() + " Finishing branch " + cur_branch);
    #print(cur_branch_dict)

    return cur_dict, branch_outlines;

#############################################################################