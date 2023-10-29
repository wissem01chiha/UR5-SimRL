function fe = fele(x1,x2,q0,L)
he=x2-x1;
%fonction qui donne le vecteur de force emlenetaire  de la poutre 
 
%modelisation de la charge repartie q(x)
q=[q0/L 0];
 
 f1=diff(polyval(polyint(conv([2/he.^3 -3/he.^3 0 1 ],q)),[0 he]));
 f2=diff(polyval(polyint(conv([1/he.^2 -2/he 1 0 ],q)),[0 he]));
  f3=diff(polyval(polyint(conv([-2/he.^3 3/he.^2 0 0],q)),[0 he]));
  f4=diff(polyval(polyint(conv([1/he.^2 -1/he 0 0],q)),[0 he]));
  
  fe=[ f1 f2 f3 f4 ]';
 
 
 











 
end

