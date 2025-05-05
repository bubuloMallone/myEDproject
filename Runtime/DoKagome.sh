#!/bin/bash

thr=4
Ns=9
t=-1
Nmax=2
Ldim=$((Nmax + 1))
V=0  

tot_time=100.0
Dt=0.1 

# Simulation parameters
U=1.0
max=2.0
step=0.05


# Select the appropriate kset
if [ $Ns -eq 9 ]; then
    Nn="${Ns}D6"
    kset="Gamma.D6.A1 Gamma.D6.A2 Gamma.D6.B1 Gamma.D6.B2 Gamma.D6.E1 Gamma.D6.E2 K.D3.A1 K.D3.A2 K.D3.E"
elif [ $Ns -eq 12 ]; then
    Nn="${Ns}D6"
    kset="Gamma.D6.A1 Gamma.D6.A2 Gamma.D6.B1 Gamma.D6.B2 Gamma.D6.E1 Gamma.D6.E2 M.D2.A1 M.D2.A2 M.D2.B1 M.D2.B2"
elif [ $Ns -eq 21 ]; then
    Nn="${Ns}C6"
    kset="Gamma.C6.A Gamma.C6.B Gamma.C6.E1a Gamma.C6.E1b Gamma.C6.E2a Gamma.C6.E2b 0.C1.A"
elif [ $Ns -eq 24 ]; then
    Nn="${Ns}D2"
    kset="Gamma.D2.A1 Gamma.D2.A2 Gamma.D2.B1 Gamma.D2.B2 M.D2.A1 M.D2.A2 M.D2.B1 M.D2.B2 X.D1.A X.D1.B 0.C1.A"
elif [ $Ns -eq 27 ]; then
    Nn="${Ns}D6"
    kset="Gamma.D6.A1 Gamma.D6.A2 Gamma.D6.B1 Gamma.D6.B2 Gamma.D6.E1 Gamma.D6.E2 K.D3.A1 K.D3.A2 K.D3.E Y.D1.A Y.D1.B"
else
    echo "Error: Unsupported number of sites, Ns=$Ns"
    exit 1
fi

# Possible values for Nb
Nb_values=($Ns $(($Ns-1)) $(($Ns+1)) $(($Ns-2)) $(($Ns+2)) $(($Ns-3)) $(($Ns+3)) $(($Ns-4)) $(($Ns+4)))

# for (( Nmax=2; Nmax<=Ns; Nmax++ )); do
# Ldim=$((Nmax + 1))

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

# Begin simulation loop
while (( $(bc <<< "$U <= $max") == 1 )); do
# Properly format U to ensure the decimal point and zero are included
formatted_U=$(printf "%.2f" "$U")
for r in $kset; do
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
LatticeDescriptionFile="LatticeFiles/kagome.BoseHubbard.$Nn.lat";
H1aFilename="Matrix/Bosons_$Nmax.U";
H1aFactor=$formatted_U;
H2aFilename="Matrix/Bosons_$Nmax.hop";
H2aFactor=$t;
H2bFilename="Matrix/Bosons_$Nmax.nn";
H2bFactor=$V;
Nevals=1;
xNevecs=1;
initthreads=$thr;
threads=$thr;
xPiecewiseConstantTimeEvolution=1;
xTimeEvolution=1;
xTEFunctionFile="TimeEvolutionFiles/TimeEvoFile_Umax.1.0.t.$tot_time.Dt.$Dt.txt";
xMeasurementProtocol="CorrsFiles/corrsFile.$Ns";
xStartingVectorFilename="Gs/Bosons_$Nmax.kagome.$Nn.t.$t.U.$formatted_U.N.$Nb.k.$r.gs";
xOutFile="DummyOut.dat";
xTEOutFile="OutTime/TE.Bosons_$Nmax.kagome.$Nn.t.$t.U.$formatted_U.T.$tot_time.DT.$Dt.dat";
xenergies=1;
xOneBodyCorrelations=1;
xTwoBodyCorrelations=1;
xCorrelationsFileName="CorrsFiles/corrsFile_$Nmax.$Ns";
xdumpgs=1;
xgsfile="Gs/Bosons_$Nmax.kagome.$Nn.t.$t.U.$formatted_U.N.$Nb.k.$r.gs";
OutFile="OutKagome/Energies.Bosons_$Nmax.kagome.$Nn.t.$t.U.$formatted_U.Nb.$Nb.rep.$r.dat";
xcheckpointing=1;
xcheckpointpath="Cp/Bosons_$Nmax.kagome.$Nn.t.$t.U.$formatted_U.N.$Ns.k.0.A";
xmoduloiters=10;
EOFm
done
done
U=$(bc <<< "$U + $step")
done

# done