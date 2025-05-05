BEGIN {
    local_dim=2;
    print local_dim*local_dim;

    for (i = 0; i < local_dim; i++) {
        print i, (local_dim-1)-i, 1;
    }
}