#!/bin/bash
# Script to generate .P files for a range of Nmax values

for Nmax in 9  # Adjust the range as needed
do
    for n in $(seq 0 $Nmax)
    do
        awk -f P.awk $Nmax $n > "Bosons_${Nmax}.P_${n}"
    done
done