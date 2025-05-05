BEGIN {
 Nmax=ARGV[1];
 dim=Nmax+1;
 ARGC=1;
 print (Nmax+1)*(Nmax+1);

 for(ni=0;ni<=Nmax;++ni) {
  for(nj=0;nj<=Nmax;++nj) { 
   # bdag b
   if(ni<Nmax && nj>0) {
    amplitude=sqrt(ni+1)*sqrt(nj);
    print ni*dim+nj,(ni+1)*dim+(nj-1),amplitude; 
   }
  }
 }
}
