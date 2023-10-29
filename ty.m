function  Ty= ty(E,I,q0,F0,n,l )
 %effort tranchant 
 x=0:l/(n):l;
 tt=-1*polyder(mfz(E,I,q0,F0,l,n));
 
Tth=F0+(q0*l/2).*(1-(x/l).^2);
plot(x,Tth,x,polyval(tt,x)),grid 
end

