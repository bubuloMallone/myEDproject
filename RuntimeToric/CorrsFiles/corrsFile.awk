BEGIN {
    Ns=ARGV[1];
    Nmax=ARGV[2];
    ARGC=1;
    
    for (site = 0; site < Ns; site++) {
        print "ONEBODYCORRELATION";
        print "Matrix/Bosons_" Nmax ".n";
        print site;
        print "Corr1Body_bdagb";
    }

    for (site1 = 0; site1 < Ns; site1++) {
        for (site2 = 0; site2 < Ns; site2++) {
            if (site1 != site2){
                print "TWOBODYCORRELATION";
                print "Matrix/Bosons_" Nmax ".bdagb";
                print site1 " " site2;
                print "Corr2Body_bdagb";
            }
        }
    }
}
