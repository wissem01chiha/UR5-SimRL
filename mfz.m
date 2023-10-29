function  Mfz = mfz(E,I,q0,F0,l,n)
u=sol(E,I,q0,F0,l,n);
 dep=zeros(n+1,1);
 for i=1:n+1
     dep(i,1)=u(2*i-1,1);
 end 
 x2=0:0.001:l;
x=0:l/n:l;
 
  d=polyfit(x,dep',n+2);
 
 %solution theorique exacte 
   Mz=F0*(l-x2)+(q0*(l.^2)/6)*(2-(3/l).*x2+(x2./l).^3);
   %solution pratique 
   Mfz=zeros(n+1,1);
 for i=1:n+1
     Mfz(i,1)=E*I*polyval(polyder(polyder(d)),(i-1)*l/n);
 end
    mf=polyfit(x,Mfz',n+2);
 plot( x2,polyval(mf',x2),x2,Mz ),grid
 
 
 










end

