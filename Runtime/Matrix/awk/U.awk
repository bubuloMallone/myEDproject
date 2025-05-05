BEGIN {
    Nmax = ARGV[1];
    print Nmax + 1;

    for (n = 0; n <= Nmax; n++) {
        interaction = int(n * (n - 1) / 2);
        print n, n, interaction;
    }
}

