#!/bin/bash

# Set the OpenMP flags for this session
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"

# Simulation parameters
thr=2
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=$thr

# Display confirmation of OpenMP variables
echo "OpenMP environment variables:"
echo "OPENBLAS_NUM_THREADS = $OPENBLAS_NUM_THREADS"
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo "LDFLAGS = $LDFLAGS"
echo "CPPFLAGS = $CPPFLAGS"

Ns=16
Ldim=2
hx_min=0.03
hx_max=0.03
stepx=0.01
hx_values=$(seq $hx_min $stepx $hx_max)
hy=0.00
hy_form=$(printf "%.2f" "$hy")
Nevals=100

# time evolution parameters
# tot_time=0.001
# Dt=0.000001

# Begin simulation loop
for hx in $hx_values; do
hx_form=$(printf "%.2f" "$hx") 
hz_form=$(printf "%.2f" "$hx")
echo "Running simulation with hx=$hx_form hz=$hz_form"
python3 /Users/pietro/myEDproject/EDroutine/main.py <<-EOFm
threads=$thr;
sparseMethod=1;
LatticeDescriptionFile="LatticeFiles/TC.mixedfields.$Ns.lat";
hxFilename="Matrix/Spin.X";
hxFactor=-$hx_form;
hzFilename="Matrix/Spin.Z";
hzFactor=-$hz_form;
hyFilename="Matrix/Spin.Y";
hyFactor=$hy_form;
PAFilename="Matrix/Spin.star";
PAFactor=-1.0;
PBFilename="Matrix/Spin.plaq";
PBFactor=-1.0;
Nevals=$Nevals;
OutFile="OutToric/Energies.Toric_$Ns.square.hx.$hx_form.hz.$hz_form.dat";
# OutFile="OutToric/Energies.DUMMY.dat";
# dumpgs=1;
# dumpgsFile="gsToric/GS.Toric_$Ns.square.hx.$hx_form.hz.$hz_form.dat";
EOFm
done