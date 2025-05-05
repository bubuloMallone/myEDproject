BEGIN {
 Nmax=ARGV[1];
 dim=Nmax+1;
 ARGC=1;
 print (Nmax+1)*(Nmax+1);

 for(ni=0;ni<=Nmax;++ni) {
  for(nj=0;nj<=Nmax;++nj) { 
   # n n aka density density 
    amplitude=ni*nj;
    print ni*dim+nj,ni*dim+nj,amplitude; 
  }
 }
}
