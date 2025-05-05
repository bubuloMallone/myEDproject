#!/bin/csh
@ thr = 1
@ Ns = 9
@ Nmax = 3
@ Ldim = $Nmax + 1
set Nn = "${Ns}D6"
set kset = "Gamma.D6.A1 " # Gamma.D6.A2 Gamma.D6.B1 Gamma.D6.B2 Gamma.D6.E1 Gamma.D6.E2 K.D3.A1 K.D3.A2 K.D3.E M.D2.A1 M.D2.A2 M.D2.B1 M.D2.B2 X.D1.A X.D1.B"
set t = -1
set V = 0

# Determine the required Nbits based on Ldim
@ Nbits = 1
@ power = 2  # Start with 2^1
while ($power < $Ldim)
    @ Nbits++
    @ power = $power * 2
end

# Make sure the number of bits is large enough to contain the states of the local Hilbert space

# Initial value
set U = 0.0
# Upper limit
set max = 1.0
# Step size
set step = 0.1

while (`echo "$U <= $max" | bc` == 1)
# Increment U by step
set U = `echo "$U + $step" | bc`
foreach r ( $kset )
# @ Nm1 = $Ns - 1
# @ Np1 = $Ns + 1

foreach Nb ($Ns) # $Nm1 $Np1 )

/Users/pietro/SourceCode/sourcecodeed/ED_Suite/EDEngine_1_clang++ <<EOFm
geometry="GeneralLattice";
sites=$Ns;
StatesPerSite=$Ldim;
Nbits=$Nbits;
Sz=0;
SzFile="Matrix/Bosons_$Nmax.Sz";
charge=$Nb;
ChargeFile="Matrix/Bosons_$Nmax.charge";
representation="$r";
LatticeDescriptionFile="LatticeFiles/triangular.BoseHubbard.$Nn.lat";
H1aFilename="Matrix/Bosons_$Nmax.U";
H1aFactor=$U;
H2aFilename="Matrix/Bosons_$Nmax.hop";
H2aFactor=$t;
H2bFilename="Matrix/Bosons_$Nmax.nn";
H2bFactor=$V;
Nevals=1;
initthreads=$thr;
threads=$thr;
xdumpgs=1;
energies=1;
xgsfile="Gs/Bosons_$Nmax.triangular.$Nn.t.1.U.$U.N.$Nb.k.$r.gs";
OutFile="Out/Energies.Bosons_$Nmax.triangular.$Nn.t.$t.U.$U.Nb.$Nb.rep.$r.dat";
xcheckpointing=1;
xcheckpointpath="Cp/Bosons_$Nmax.triangular.$Nn.t.1.U.$U.N.$Ns.k.0.A";
xmoduloiters=10;
EOFm
end
end
end
