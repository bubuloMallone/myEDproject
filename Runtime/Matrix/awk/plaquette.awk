BEGIN {
 local_dim=2;
 n_spins=4;
 # tot_dim=pow(local_dim,n_spins);
 tot_dim=local_dim ** n_spins;
 print (tot_dim)*(tot_dim);

 for(i=0;i<tot_dim;++i) {
    # Convert i to binary and count 1s
    num_ones = 0;
    x = i;
    while (x > 0) {
        num_ones += x % 2;  # Check if the last bit is 1
        x = int(x / 2);  # Remove the last bit by dividing by 2
    }
    if(num_ones % 2 == 0){
        value = 1;
    } 
    else {
        value = -1;
    }
    print i, i, value;
 }
}