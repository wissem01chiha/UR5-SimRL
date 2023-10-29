function  ke = kele(x1,x2,E,I)
  %definition du longeur de l'elemnt 
he=x2-x1;
 
%matrice ke

 ke=zeros(4,4);
 for i=1:4
     for j=1:4
         ke(i,j)=E*I*diff(polyval(polyint(conv(derive(i,he),derive(j,he))),[0 he]));
     end
 end

end

