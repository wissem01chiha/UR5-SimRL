function   f= derive(k,he)
 
switch(k)
     
    case(1)
         x=0:0.1:he;
        f=  polyder(polyder([2/(he.^3) -3/(he.^3) 0 1 ]));
        y=polyval(f,x);
        %plot(x,y),grid
    case(2)
       x=0:0.1:he;
        f=polyder(polyder([1/he.^2 -2/he 1 0 ]));
           y=polyval(f,x);
       % plot(x,y),grid
    case(3)
           x=0:0.1:he;
        f=polyder(polyder([-2/he.^3 3/he.^2 0  0 ]));
           y=polyval(f,x);
        %plot(x,y),grid
    case(4)
           x=0:0.1:he;
        f=polyder(polyder([1/he.^2  -1/he  0 0 ]));
           y=polyval(f,x);
      %plot(x,y),grid

end


end

