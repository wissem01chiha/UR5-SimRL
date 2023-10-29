function k = Kglob(l,n,E,I)
%matrice de rigidite assemblee pour le premeier mailllage (regulier)
 k=zeros(2*n+2,2*n+2);
for  i=1:n
 ke=kele((i-1)*l/n,i*l/n,E,I);
s=zeros(4,2*n+2);
  s(1,2*i-1)=1;
  s(2,2*i)=1;
  s(3,2*i+1)=1;  
 s(4,2*i+2)=1;
 ki=s'*ke*s;
 k=k+ki;
 
end

