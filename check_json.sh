#!/bin/bash
#PURPOSE: Check if any Python incompatible strings (e.g. inf, null, -nan) are in any .json files in a directory

if [ "$#" -ne 1 ]; then
    echo "Error: an input directory command line argument is required"
    exit 1
fi

indir=$1
echo "Input directory: ${indir}"
for file in $indir/*.json; do
    #echo $file;
    grep -l 'inf,' $file;
    grep -l 'null' $file;
    grep -l '\-nan' $file;
done