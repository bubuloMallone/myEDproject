BEGIN {
    Nmax = ARGV[1];
    for (n = 0; n <= Nmax; n++) {
        print n, n; 
    }
}
