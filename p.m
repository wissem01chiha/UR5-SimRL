 
function  y=p(i,he)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
switch(i)
     case (1)
         x=0:1/100:he;
         y=1-3*(x.^2)/(he^2)+(2*x.^3)/(he^3);
         plot(x,y) ,grid
         title("fonction de forme psi1")
     case(2)
         x=0:1/100:he;
         y=x-2*(x.^2)/he+(x.^3)/(he^2);
         plot(x,y),grid
          title("fonction de forme psi2")
     case(3)
         x=0:1/100:he;
         y=3*x.^2/(he^2)-2*(x.^3)/he^3;
         plot(x,y) , grid
          title("fonction de forme psi3")
     case(4)
         x=0:1/100:he;
         y=-1*x.^2/he+(x.^3)/(he^2);
         plot(x,y) ,grid
          title("fonction de forme psi4")
     
 end

