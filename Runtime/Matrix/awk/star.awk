BEGIN {
 local_dim=2;
 n_spins=4;
 # tot_dim=pow(local_dim,n_spins);
 tot_dim=local_dim ** n_spins;
 print (tot_dim)*(tot_dim);

 for(i=0;i<tot_dim;++i) {
    print i,(tot_dim-1)-i,1;
 }
}