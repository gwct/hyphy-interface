# HyPhy interface
### Scripts that generate commands for various HyPhy models and parse its output

## About
[HyPhy](https://www.hyphy.org/) is a program to estimate rates of molecular evolution in coding sequences along a phylogeny. It has various models to perform tests for selection. In many cases, we might want to run these models on many hundreds or thousands of genes, each with its own tree. These scripts endeavor to make this an easy process by generating commands and SLURM submission scripts to run HyPhy on individual genes, and then to subsequently parse the output of HyPhy.

A subtree of various [HyPhy analyses](https://github.com/veg/hyphy-analyses) is included for some models (FitMG94).


#### Generating commands

Generally, to generate commands for a model:

```
python hyphy_gen.py -i [directory with codon alignments in FASTA format] -m [HyPhy model] -genetrees [directory with gene trees] -o [output directory]
```

File must have a matching ID between the alignments and the gene trees. use `-s` to define how to split file names.

Specific models may have more options available. Use `-h` to see the full list of options.

This command will generate both a shell script with one HyPhy command per line, and a SLURM submission script that uses [GNU Parallel](https://www.gnu.org/software/parallel/) to run the commands in the shell script in parallel.

#### Parsing output

To parse output into a .csv file that can be easily read into statistical software:

```
python hyphy_to_csv.py -i [directory containing hyphy json files (-o from hyphy_gen.py)] -m [HyPhy model] -o [output file name]
```

#### Models supported

Current [HyPhy models](https://www.hyphy.org/methods/) supported are mg94, mg94-local, fel, busted, fubar, absrel, slac, and relax

## Methods for estimating branch-specific rates while accounting for discordance

Oftentimes the phylogeny for an individual gene will not match the species tree. In these cases, forcing the gene's sequence onto the species tree will result in the mis-mapping of substitutions and could lead to spurious results when estimating rates and testing for selection. Here we implement a method for averaging rates for each branch in a species tree based on gene tree topologies.

Briefly, instead of using a single species tree for all genes, or only estimating rates when the gene tree matches the species tree, we estimate rates using each gene's tree. Then, for every branch in the species tree we ask whether that branch exists in the gene tree or not. If not, we skip that branch for that gene. If so, we look at the number of sites and number of substitutions in that gene on that branch. Over all genes, we sum these values and calculate dN and dS from these sums.

Notably this only works for methods where the number of sites and number of substitutions are reported per-branch for every gene (like [SLAC](https://www.hyphy.org/methods/selection-methods/#slac) or PAML's M1 model).

This method is being developed in the branch-rates folder. 

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

## Acknowledgements

This material is based upon work supported by the National Science Foundation under Grant Number (DEB 1754096).

Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.