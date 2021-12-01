# Branch specific rate calculations of dN, dS, and dN/dS while accounting for gene-tree discordance

## 1. Run SLAC with gene trees

To run this method, first generate the SLAC commands and run the SLAC model:

```
python hyphy_gen.py -i [directory with codon alignments in FASTA format] -m slac -genetrees [directory with gene trees] -o [output directory]
```

## 2. Convert SLAC output to .csv

Then, parse the output to .csv format:

```
python hyphy_to_csv.py -i [directory containing hyphy json files (-o from hyphy_gen.py)] -m slac -o [output file name]
```

In addition to the single .csv file for all loci, this will create a directory within the directory defined by `-i` called csv. This contains one .csv file per locus with values parsed out and formatted by branch. This will be used to compute branch specific rates.

## 3. Calculate branch rates

For now, first move to the branch-rates directory:

```
cd branch-rates
```

and then run:

```
python branch_rates.py -i [tree file in Newick format] -r [directory containing .csv files from a SLAC run] -o [output directory]
```

The rates directory (`-r`) will be the csv directory created by hyphy_to_csv.py within the directory with the SLAC .json files.

Within the output directory, a .csv file will be created based on the input tree with columns for the calculated rates. 

## 4. Reading the results into R

Since branch_rates and R label their branches differently, we need to join them up. With `get_tree_info.r` you can create another .csv file with branches labeled as R's ape package understands them, or you can read them directly within an R script.

### STEP 1: Starting an R script for branch analysis

First, load the required packages, including the `get_tree_info.r` script:

```r
library(ape)
# To read the Newick tree into R

library(dplyr)
# To join the two tables together

source("get_tree_info.r")
# Only needed for Option 2, to read the tree table directly in R
```

Next, read your tree into R with the ape function `read.tree()`:

```
tree = read.tree([file with species tree in Newick format])
```

### STEP 2 (OPTION 1): create another table file with ape branch labels

Outside of R, run this command. This is the only command not run in R in this walkthrough:

```
Rscript get_tree_info.r -t [file with species tree in Newick format] -o [output .csv file name]
```

Then, in your R script read this table:

```r
tree_info = read.csv([.csv file from get_tree_info.r], header=T)
```

### STEP 2 (OPTION 2): read the tree info table directly into your R script

Pass the tree you read earlier into the treeToDF function, and unpack the table:

```r
tree_to_df_list = treeToDF(tree)
tree_info = tree-to_df_list[["info"]]
```

### STEP 3: Join the R tree info table and the branch_rates table:

Either way, you should end up with a tree table called `tree_info` that contains the branch labels for your tree as ape labels them. Now we have to join this with the information from branch_rates.

First, read the branch_rates table into your script:

```r
tree_rates = read.csv([.csv file from branch_rates.py], header=T)
```

Then merge the two tables by clade, discarding any duplicate columns:

```r
uniq_info_cols = names(tree_info)[!(names(tree_info) %in% names(tree_rates))] 
# Get non common names

uniq_info_cols = c(uniq_info_cols, "clade") # appending key parameter
# Get a list of columns from the tree_info df to join to the tree rates df

tree_rates = tree_rates %>% left_join((tree_info %>% select(uniq_info_cols)), by="clade")
# Select the columns from tree_info and join to tree_rates, merging by clade
# https://stackoverflow.com/a/61628157
```

Finally, reorder the rows so they are in the order that ape labels branches:

```r
tree_rates = tree_rates[order(tree_rates$node), ]
# Re-order the rows by the R node
```

You are now free to analyze rates by branch!

## Directory overview

| Directory | Description | 
| ------ | ----------- |
| lib/ | General helper functions and code for branch calculations |
| branch_rates.py | Interface script for method being developed to calculate rates per branch in a species tree while accounting for discordance |
