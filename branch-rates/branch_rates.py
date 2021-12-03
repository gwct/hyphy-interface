#!/usr/bin/env python3
#############################################################################
# Computes branch averages for rates from a directory of
# hyphy results as csv files.
# Modified from Gregg Thomas's script (31 August 2021) to parse ES, EN, S, and N instead of dN/dS (EK)
# 10.29.2021: Reformatted into a more project like script (GT)
#############################################################################

# time -p branch_avgs.py -i docs/data/trees/full-coding-astral-cf-rooted.csv -r 05-MolEvol/full-coding-slac/csv/ -s /mnt/beegfs/ek112884/murinae/DNDS_BY_SPERM_STAGE/HYPHY/SLAC/prot_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.txt -f 05-MolEvol/full-coding-mg94-local-ds-filter-0.95quant.csv -o 05-MolEvol/full-coding-slac/branch-avgs-astral-lz-induced/ -n 63 --overwrite

#############################################################################

import sys
import os
import multiprocessing as mp
import lib.core as CORE
import lib.params as params
import lib.opt_parse as OP
import lib.tree as TREE
import lib.branches as BRANCHES

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.
    
    print("\n" + " ".join(sys.argv) + "\n");

    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# branch_avgs version " + globs['version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.

    print("#");
    print("# " + "=" * 125);
    #print(CORE.welcome());
    #if "-h" not in sys.argv:
    #    print("            Degeneracy annotation of transcripts\n");
    # A welcome banner.

    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    ##########################

    step = "Reading files in input directory";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    globs['csv-files'] = [ f for f in os.listdir(globs['csv-rate-dir']) if f.endswith(".csv") ];
    # Read the filtered genes to skip them for averaging

    step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + str(len(globs['csv-files'])) + " files read");
    # Status update

    ##########################

    if globs['filter-file']:
        step = "Reading filter file";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        # Status update

        globs['filter-files'] = [ line.replace("\"", "").replace(".json", ".csv") for line in open(globs['filter-file']) ];
        # Read the filtered genes to skip them for averaging

        globs['csv-files'] = list( set(globs['csv-files']) - set(globs['filter-files']) );
        # Remove the files to filter from the list of input files

        step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + str(len(globs['filter-files'])) + " loci will be excluded");
        # Status update
    
    ##########################

    if globs['subset-file']:
        step = "Reading subset file";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        # Status update

        globs['subset-files'] = [ line.replace("\n","-mafft-cds.filter.csv") for line in open(globs['subset-file']) ];
        # Read in a subset of genes to include in analysis
        ## TODO: generalize this

        globs['csv-files'] = list( set(globs['csv-files']) & set(globs['subset-files']) );
        # Only take files in both lists

        step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + str(len(globs['subset-files'])) + " loci read");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + str(len(globs['csv-files'])) + " loci included in final counts (intersect of input files and subset list).")
        # Status update

    ##########################

    step = "Writing list of included loci";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    loci_file = os.path.join(globs['outdir'], "counted-loci.txt");
    with open(loci_file, "w") as locifile:
        for f in globs['csv-files']:
            locifile.write(f + "\n");
    step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: counted-loci.txt written");
    # Status update

    ########################## 

    step = "Reading species tree";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    globs['orig-tree-string'] = open(globs['tree-file']).read();

    try:
        globs['tree-dict'], globs['tree-string'], globs['root'] = TREE.treeParse(globs['orig-tree-string'])
    except:
        CORE.errorOut("M1", "Error reading species tree!", globs);

    step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS");
    # Status update

    # Read the tree
    ##########################

    step = "Initializing main output dictionary";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update   

    branches = {};

    headers = ["node.label", "clade", "node.type", "orig.node.label"];
    if globs['rooted']:
        count_headers = ["num.genes.full", "num.genes.partial", "num.genes.descendant.counted", "num.genes.discordant", "num.genes.missing", "cS", "cN", "cA", "mS", "mN", "mA"];
    else:
        count_headers = ["num.genes.full", "num.genes.partial", "num.genes.descendant.counted", "num.genes.discordant", "num.genes.missing", "ES", "EN", "S", "N"];
        avg_headers = ["avg.ES", "avg.EN", "avg.S", "avg.N", "dS", "dN", "dNdS"];
    headers += count_headers;
    # The headers/keys for each node in the species tree, depending on whether this is from a rooted ancestral sequences run or SLAC

    for node in globs['tree-dict']:
    # Setup the dictionary for each node in the species tree

        branches[node] = {};
        # Add the current node to the tracker dict

        cur_clade = sorted(TREE.getClade(node, globs['tree-dict']));
        branches[node]["clade"] = set(cur_clade);
        # Get the clade descending from the current node

        branches[node]["node.label"] = node;
        branches[node]["node.type"] = globs['tree-dict'][node][2];
        branches[node]["orig.node.label"] = globs['tree-dict'][node][3];
        # Get the other info for the current node

        for h in count_headers:
            branches[node][h] = 0.0;
        # Initialize all the counts for the current node

    step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS");
    # Status update

    # Initialize main dict
    ##########################

    step = "Splitting input files into chunks";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    chunks = CORE.chunks(globs['csv-files'], globs['chunk-size']);
    num_chunks = str(len(chunks));

    step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + num_chunks + "  chunks of size " + str(globs['chunk-size']));
    # Status update

    # Split the files into chunks
    ##########################

    step = "Summing values over branches";
    step_start_time = CORE.report_step(globs, step, False, "Processed 0 / " + num_chunks + " chunks...", full_update=True);
    # Status update

    # branch_num = 1;
    # for branch in branches:
    #     if branch != "Kadarsanomys_sodyi_MZB-Sample":
    #         continue;

    #     result = BRANCHES.branchSum((globs['branches'][branch], branch, globs['csv-rate-dir'], globs['filter-files'], globs['subset-files']));
    #     new_branches[result[0]] = result[1];
    #     #branches.update(result);
    #     cur_branch_time = CORE.report_step(globs, step, step_start_time, "Processed " + str(branch_num) + " / " + str(globs['num-branches']) + " branches...", full_update=True);
    #     branch_num += 1;
    ## Serial version for debugging

    with mp.Pool(processes=globs['num-procs']) as pool:
        chunk_num = 1;
        new_branches = {};
        for result in pool.imap_unordered(BRANCHES.branchSumNEW, ((chunk, globs['csv-rate-dir'], globs['tree-dict'], globs['rooted']) for chunk in chunks)):
            
            if result == "stop":
                sys.exit();

            for node in globs['tree-dict']:
                for h in count_headers:
                    branches[node][h] += result[node][h];
            # Sum the result

            cur_chunk_time = CORE.report_step(globs, step, step_start_time, "Processed " + str(chunk_num) + " / " + num_chunks + " chunks...", full_update=True);
            chunk_num += 1;
    ## Parallel version

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success", full_update=True);
    # Status update
    # sys.exit();
    # Get number of subs for each branch from each gene
    ##########################

    step = "Writing out new tree table";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    if not globs['rooted']:
        
        headers += avg_headers;

        for node in branches:
        # For every clade in branches, add the new headers into the branches dictionary

            for ah in avg_headers:
                if branches[node]["node.type"] == "ROOT":
                    branches[node][ah] = "NA";
                else:
                    branches[node][ah] = 0.0;
            # Add in the new headers and initialize at 0, except for the root node

    with open(globs['output-file'], "w") as outfile:
        outfile.write(",".join(headers) + "\n")
        # Open the output file and write the headers, which now contain the new columns

        for node in branches:
            if not globs['rooted']:
                if branches[node]["node.type"] != "ROOT" and (branches[node]['num.genes.full'] + branches[node]['num.genes.partial']) != 0:
                    branches[node]['avg.ES'] = branches[node]['ES']  / (branches[node]['num.genes.full'] + branches[node]['num.genes.partial']);
                    branches[node]['avg.EN'] = branches[node]['EN']  / (branches[node]['num.genes.full'] + branches[node]['num.genes.partial']);
                    branches[node]['avg.S'] = branches[node]['S']  / (branches[node]['num.genes.full'] + branches[node]['num.genes.partial']);
                    branches[node]['avg.N'] = branches[node]['N']  / (branches[node]['num.genes.full'] + branches[node]['num.genes.partial']);
                # If the branch is not the root and appears in some genes, then compute the averages.

                    if branches[node]['ES']:
                        branches[node]['dS'] = branches[node]['S'] / branches[node]['ES'];
                        
                    if branches[node]['EN']:
                        branches[node]['dN'] = branches[node]['N'] / branches[node]['EN'];

                    if branches[node]['dS']:
                        branches[node]['dNdS'] = branches[node]['dN'] / branches[node]['dS'];
                    else:
                        branches[node]['dNdS'] = "NaN";

                    # new_branches[branch]['avg.dS'] = new_branches[branch]['avg.S'] / new_branches[branch]['avg.ES']
                    # new_branches[branch]['avg.dN'] = new_branches[branch]['avg.N'] / new_branches[branch]['avg.EN']
                    # new_branches[branch]['avg.dNdS'] = new_branches[branch]['avg.dN'] / new_branches[branch]['avg.dS']
            # After averaging ES, EN, S, and N, use them to calculate dS, dN, and dN/dS across all alignments

            branches[node]["clade"] = ";".join(sorted(list(branches[node]["clade"]), key=str.casefold));
            # Sorts the speceies in the clade, while ignoring the case of the letters.

            outline = [];
            for h in headers:
                outline.append(str(branches[node][h]));
            
            outline = ",".join(outline);
            outfile.write(outline + "\n");
            # Compile the output line and write it out to the file.

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    # Output the tree table with the new columns
    ##########################

    CORE.endProg(globs);

#############################################################################

