#!/bin/bash

thr=4
Ns=26
Ldim=2

# Simulation parameters
hx_min=0.00
hx_max=1.00
stepx=0.05

hz_min=0.00
hz_max=0.00
stepz=0.10

# Define the range of hx and hz values
hx_values=$(seq $hx_min $stepx $hx_max)  # Generates a sequence of hx values
hz_values=$(seq $hz_min $stepz $hz_max)  # Generates a sequence of hz values

tot_time=0.001
Dt=0.000001

# Select the appropriate kset
if [ $Ns -eq 10 ]; then
    kset="Gamma.C4.A Gamma.C4.B Gamma.C4.Ea Gamma.C4.Eb 0.C1.A"
elif [ $Ns -eq 16 ]; then
    kset="Gamma.D4.A1 Gamma.D4.A2 Gamma.D4.B1 Gamma.D4.B2 Gamma.D4.E M.D4.A1 M.D4.A2 M.D4.B1 M.D4.B2 M.D4.E Sigma.D1.A Sigma.D1.B X.D2.A1 X.D2.A2 X.D2.B1 X.D2.B2"
elif [ $Ns -eq 18 ]; then
    kset="Gamma.D4.A1 Gamma.D4.A2 Gamma.D4.B1 Gamma.D4.B2 Gamma.D4.E Delta.D1.A Delta.D1.B Sigma.D1.A Sigma.D1.B"
elif [ $Ns -eq 20 ]; then
    kset="Gamma.C4.A Gamma.C4.B Gamma.C4.Ea Gamma.C4.Eb M.D4.B1 M.D4.B2 M.D4.E 0.C1.A 1.C1.A"
elif [ $Ns -eq 26 ]; then
    kset="Gamma.C4.A Gamma.C4.B Gamma.C4.Ea Gamma.C4.Eb 0.C1.A 1.C1.A 2.C1.A"
else
    echo "Error: Unsupported number of sites, Ns=$Ns"
    exit 1
fi

# Determine the required Nbits based on Ldim
Nbits=1
power=2  # Start with 2^1
while [[ $power -lt $Ldim ]]; do
    Nbits=$((Nbits + 1))
    power=$((power * 2))
done

# Select the appropriate executable based on Nbits
if [ $Nbits -eq 1 ]; then
    executable="/Users/pietro/SourceCode/sourcecodeed/ED_Suite/EDEngine_1_clang++"
elif [ $Nbits -eq 2 ]; then
    executable="/Users/pietro/SourceCode/sourcecodeed/ED_Suite/EDEngine_2_clang++"
elif [ $Nbits -eq 3 ]; then
    executable="/Users/pietro/SourceCode/sourcecodeed/ED_Suite/EDEngine_3_clang++"
elif [ $Nbits -eq 4 ]; then
    executable="/Users/pietro/SourceCode/sourcecodeed/ED_Suite/EDEngine_4_clang++"
else
    echo "Error: Unsupported number of bits, Nbits=$Nbits"
    exit 1
fi

echo "Running executable: $executable"

# Begin simulation loop
for hx in $hx_values; do
hx_form=$(printf "%.2f" "$hx")  # Format hx to two decimal places
for hz in $hz_values; do
hz_form=$(printf "%.2f" "$hz")  # Format hz to two decimal places
for r in $kset; do
echo "Running simulation with hx=$hx_form hz=$hz_form representation=$r"
$executable <<-EOFm
geometry="GeneralLattice";
sites=$Ns;
StatesPerSite=$Ldim;
Nbits=$Nbits;
xSz=0.5;
SzFile="Matrix/SpinHalf.Sz_nc";
xcharge=$Ns;
xChargeFile="Matrix/SpinHalf.charge";
representation="$r";
LatticeDescriptionFile="LatticeFiles/TC.mixedfields.$Ns.lat";
hxFilename="Matrix/Spin.X";
hxFactor=-$hx_form;
hyFilename="Matrix/Spin.Y";
hyFactor=0.0;
hzFilename="Matrix/Spin.Z";
hzFactor=-$hz_form;
PAFilename="Matrix/Spin.star";
PAFactor=-1;
PBFilename="Matrix/Spin.plaq";
PBFactor=-1;
Nevals=1;
xNevecs=1;
xDetermineRawDimension=1;
xmaxiterations=100000;
initthreads=$thr;
threads=$thr;
xPiecewiseConstantTimeEvolution=1;
xTimeEvolution=1;
xTEFunctionFile="TimeEvolutionFiles/TimeEvoFile_Umax.1.0.t.$tot_time.Dt.$Dt.txt";
xMeasurementProtocol="CorrsFiles/corrsFile.$Ns";
xStartingVectorFilename="Gs/GS.dummy";
xOutFile="DummyOut.dat";
xTEOutFile="OutTime/TE.dummy";
xenergies=1;
xDensityMatrix=1;
xPrintDensityMatrix=1;
xOneBodyCorrelations=1;
xTwoBodyCorrelations=1;
xCorrelationsFileName="CorrsFiles/corrsFile.$Ns";
xdumpgs=1;
xgsfile="Gs/Gs.dummy";
OutFile="OutToric/Energies.Toric_$Ns.square.hx.$hx_form.hz.$hz_form.rep.$r.dat";
xcheckpointing=1;
xcheckpointpath="Cp/Toric_$Ns.square.hx.$hx_form.hz.$hz_form.rep.$r.dat";
xmoduloiters=10;
EOFm
done
done
done






