# Branch specific rate calculations of dN, dS, and dN/dS while accounting for gene-tree discordance

## Walkthrough

### 1. Run SLAC with gene trees

To run this method, first generate the SLAC commands and run the SLAC model:

```
python hyphy_gen.py -i [directory with codon alignments in FASTA format] -m slac -genetrees [directory with gene trees] -o [output directory]
```

### 2. Convert SLAC output to .csv

Then, parse the output to .csv format:

```
python hyphy_to_csv.py -i [directory containing hyphy json files (-o from hyphy_gen.py)] -m slac -o [output file name]
```

### 3. Convert species tree from Newick to table format

Next, before you run the branch_rates.py script, you will need to convert your species tree from Newick format to .csv format. Do this using the get_tree_info.r script:

```
Rscript get_tree_info.r -t [file with species tree in Newick format] -o [output .csv file name]
```

In addition to the single .csv file for all loci, this will create a directory within the directory defined by `-i` called csv. This contains one .csv file per locus with values parsed out and formatted by branch. This will be used to compute branch specific rates.

### 4. Calculate branch rates

For now, first move to the branch-rates directory:

```
cd branch-rates
```

and then run:

```
python branch_rates.py -i [tree in table format (-o from get_tree_info.r)] -r [directory containing .csv files from a SLAC run] -o [output directory]
```

The rates directory (`-r`) will be the csv directory created by hyphy_to_csv.py within the directory with the SLAC .json files.

Within the output directory, a .csv file will be created based on the input tree table that will have additional columns for the calculated rates. This can easily be read into R along with the Newick formatted tree for analysis.

## Directory overview

| Directory | Description | 
| ------ | ----------- |
| branch-rates/ | Method being developed to calculate rates per branch in a species tree while accounting for discordance |
| dev/ | Stuff we're still working on |
| hyphy-analyses/ | A subtree of the [HyPhy analyses](https://github.com/veg/hyphy-analyses) repository for easy access to some of the models |
| lib/ | General helper functions and model specific code |
| fix_json.sh/ | A script to convert some values in HyPhy's output that are incompatible with Python's json reader |
| get_tree_info.r/ | Converts a Newick formatted tree to a table in .csv format with nodes labeled as ape reads them |
| hyphy_gen.py | Generates HyPhy commands given alignments, trees, and a model |
| hypyh_to_csv.py | Converts and combines Hyphy json output to .csv format for easy analysis |
| multiple_test_correction.r | A script to correct for multiple tests in the aBSREL model |