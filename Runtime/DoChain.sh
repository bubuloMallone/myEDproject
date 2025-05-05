#!/bin/bash

thr=4
Ns=12
Nmax=2
Ldim=$((Nmax + 1))  

t=-1

# Possible values for Nb
Nb_values=( $Ns $(($Ns-1)) $(($Ns+1)) $(($Ns-2)) $(($Ns+2)) $(($Ns-3)) $(($Ns+3)) $(($Ns-4)) $(($Ns+4)))

# Determine the required Nbits based on Ldim
Nbits=1
power=2  # Start with 2^1
while [[ $power -lt $Ldim ]]; do
    Nbits=$((Nbits + 1))
    power=$((power * 2))
done

# Select the appropriate executable based on Nbits
if [ $Nbits -eq 2 ]; then
    executable="/Users/pietro/SourceCode/sourcecodeed/ED_Suite/EDEngine_2_clang++"
elif [ $Nbits -eq 3 ]; then
    executable="/Users/pietro/SourceCode/sourcecodeed/ED_Suite/EDEngine_3_clang++"
elif [ $Nbits -eq 4 ]; then
    executable="/Users/pietro/SourceCode/sourcecodeed/ED_Suite/EDEngine_4_clang++"
else
    echo "Error: Unsupported number of bits, Nbits=$Nbits"
    exit 1
fi

# Simulation parameters
U=0.0
max=10.0
step=0.1

# Begin simulation loop
while (( $(bc <<< "$U <= $max") == 1 )); do
# Properly format U to ensure the decimal point and zero are included
formatted_U=$(printf "%.2f" "$U")
for r in $(seq 0 $((Ns - 1))); do
for Nb in "${Nb_values[@]}"; do
$executable <<-EOFm
geometry="GeneralLattice";
sites=$Ns;
StatesPerSite=$Ldim;
Nbits=$Nbits;
Sz=0;
SzFile="Matrix/Bosons_$Nmax.Sz";
charge=$Nb;
ChargeFile="Matrix/Bosons_$Nmax.charge";
representation="$r";
LatticeDescriptionFile="ChainFiles/chain.$Ns";
H1Filename="Matrix/Bosons_$Nmax.U";
H1Factor=$formatted_U;
H2Filename="Matrix/Bosons_$Nmax.hop";
H2Factor=$t;
Nevals=1;
initthreads=$thr;
threads=$thr;
xdumpgs=1;
energies=1;
xgsfile="Gs/Bosons_$Nmax.chain.$Ns.t.$t.U.$formatted_U.Nb.$Nb.k.$r.gs";
OutFile="Out/Energies.Bosons_$Nmax.chain.$Ns.t.$t.U.$formatted_U.Nb.$Nb.k.$r.dat";
xcheckpointing=1;
xcheckpointpath="Cp/Bosons_$Nmax.chain.$Ns.t.$t.U.$formatted_U.Nb.$Nb.k.$r.A";
xmoduloiters=10;
EOFm
done
done
U=$(bc <<< "$U + $step")
done