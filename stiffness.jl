"""

def kele(x1,x2,E,Igz):
    k=np.zeros((4,4))
    he=x2-x1
    x=np.linspace(0,he,100)
    for i in range(0,4):
        for j in range(i,4):
          #print(d2_psi(x1,x2,i+1))
            k[i][j]=integrate.trapz((E*Igz)*d2_psi(x1,x2,i+1)*d2_psi(x1,x2,j+1),x)
            k[j][i]=k[i][j]   
    return k

"""
function kele(x1::double,x2::double::E::double,IGz::double)::Vector{float}



    
    
end