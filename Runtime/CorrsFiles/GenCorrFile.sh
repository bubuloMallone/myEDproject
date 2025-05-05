#!/bin/bash
# Script to generate various files for a range of Nmax values

Ns=9
for (( Nmax=2; Nmax<=Ns; Nmax++ ))
do
    # Run the AWK script with current Nmax and Ns as arguments and redirect output to a file
    awk -f corrsFile.awk "$Ns" "$Nmax" > "corrsFile_${Nmax}.${Ns}"
done