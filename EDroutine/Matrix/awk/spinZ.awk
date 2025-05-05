BEGIN {
    local_dim=2;
    print local_dim*local_dim;

    for (i = 0; i < local_dim; i++) {
        if(i == 0){
            print i, i, 1;
        } 
        else {
            print i, i, -1;
        }
    }
}