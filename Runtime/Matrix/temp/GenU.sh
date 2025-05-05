#!/bin/bash
# Script to generate .U files for a range of Nmax values

for Nmax in 9 12 16  # Adjust the range as needed
do
    awk -f U.awk $Nmax > "Bosons_${Nmax}.U"
done
