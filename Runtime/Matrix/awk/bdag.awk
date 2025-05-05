BEGIN {
 Nmax=ARGV[1];
 dim=Nmax+1;
 ARGC=1;
 print dim;

 for(ni=0;ni<Nmax;++ni) {
  # for(nj=0;nj<=Nmax;++nj) { 
   # bdag
   # if(nj==(ni+1)) {
    amplitude=sqrt(ni+1);
    print (ni+1),ni,amplitude; 
   # }
  # }
 }
}