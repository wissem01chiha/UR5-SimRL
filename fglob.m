function f = fglob(n,q0,l)
 %vecteur de force asssembler  
 %premier maillage regulier 
 f=zeros(2*n+2,1);
 for i=1:n
 fe=fele((i-1)*l/n,i*l/n,q0,l);
 
 
 f(2*i-1,1)=fe(1,1)+f(2*i-1,1);
 f(2*i,1)=fe(2,1)+ f(2*i);
 f(2*i+1,1)=fe(3,1)+f(2*i+1,1);
 f(2*i+2,1)=fe(4,1)+f(2*i+2,1);
 
 
end

