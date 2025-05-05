#!/bin/bash
# Script to generate various files for a range of Nmax values

awk_dir="awk"

# for Nmax in 7 8 9 10 11 12 # Adjust the range as needed

for (( Nmax=2; Nmax<=4; Nmax++ ))
do
    # # Generate .onsite files
    # awk -f ${awk_dir}/onsite.awk $Nmax > "Bosons_${Nmax}.n"
    
    # Generate .bdag files
    awk -f ${awk_dir}/bdag.awk $Nmax > "Bosons_${Nmax}.bdag"
    
    # Generate .b files
    awk -f ${awk_dir}/b.awk $Nmax > "Bosons_${Nmax}.b"

    # # Generate .bdagb files
    # awk -f ${awk_dir}/bdagb.awk $Nmax > "Bosons_${Nmax}.bdagb"

    # # Generate .nn files
    # awk -f ${awk_dir}/nn.awk $Nmax > "Bosons_${Nmax}.nn"

    # # Generate .Sz files
    # awk -f ${awk_dir}/Sz.awk $Nmax > "Bosons_${Nmax}.Sz"

    # # Generate .charge files
    # awk -f ${awk_dir}/charge.awk $Nmax > "Bosons_${Nmax}.charge"

    # # Generate .U files
    # awk -f ${awk_dir}/U.awk $Nmax > "Bosons_${Nmax}.U"

    # # Generate .hop files
    # awk -f ${awk_dir}/boson_hop.awk $Nmax > "Bosons_${Nmax}.hop"

    # # Generate .P_n files    
    # for n in $(seq 0 $Nmax)
    # do
    #     awk -f ${awk_dir}/P.awk $Nmax $n > "Bosons_${Nmax}.P_${n}"
    # done

done
