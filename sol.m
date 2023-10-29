function   u= sol(E,I,q0,F0,l,n )
k=Kglob(l,n,E,I);
 k(1,1)=1;
 k(2,2)=1;
 k(1,2)=0;
 k(2,1)=0;
for j=1:2
for i=3:2*n+2
    k(i,j)=0;
    k(j,i)=0;
end 
end
f=fglob(n,q0,l);
f(1,1)=0;
f(2,1)=0;
Q=zeros(2*n+2,1);
Q(2*n+1,1)=F0;
f(1,1)=0;
f(1,2)=0;
u=k\(f+Q);
 dep=zeros(n+1,1);
 for i=1:n+1
     dep(i,1)=u(2*i-1,1);
 end 
x2=0:0.01:l;
x=0:l/n:l;
   
px=polyfit(x,dep',n);
pval=polyval(px,x2);
 %solution theorique exacte 
 %uth=
plot(x2,pval ),grid 


  
  
  
 
 




















 
    
    
 
 
 
 
   

 
 
 
end

