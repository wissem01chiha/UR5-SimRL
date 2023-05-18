      subroutine FORCE(Fext,ndpnu,ndpnv,ndlu,ndlv)
      implicit none
      include 'donnees.inc'
      include 'maillage.inc'      
      integer i,j,k,l,m,n,ndpnu,ndpnv,ndl,icontv,ndlu,ndlv,icontu
      real*8 FU(1000),FV(1000),Fext(*)
      
      

      FU=0.
      FV=0.
      call DZERO(Fext,ndlu+ndlv) 
      do i=1,NLOAD
       if(TYPLOAD(i).eq.'FORCE')then
        if(LOAD(2,i).eq.1)then
         icontu=(LOAD(1,i)-1)*ndpnu
         FU(icontu+1)=LOADV(i)
        elseif(LOAD(2,i).eq.2)then
         icontv=(LOAD(1,i)-1)*ndpnv
         FV(icontv+1)=LOADV(i)
        elseif(LOAD(2,i).eq.3)then
         icontv=(LOAD(1,i)-1)*ndpnv
         FV(icontv+2)=LOADV(i)        
        endif
       endif 
      enddo
      do i=1,ndlu
       Fext(i)=FU(i)
      enddo
      do i=1,ndlv
       Fext(i+ndlu)=FV(i)
      enddo        
c         stop

c        do i=1,NLOAD
c         write(*,1)i,TYPLOAD(i),LOAD(1,i),LOAD(2,i),LOADV(i)
c        enddo
c        write(*,1)ndpnu,ndpnv,ndlu,ndlv
      


      
      
1     FORMAT(162(G15.8,1X))      
      end
      

      
      
      
      
      SUBROUTINE ASSVEC(F,FELE,IM,NBN,NDPN,NDL)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %     ASSEMBLAGE DU VECTEUR ELEMENTAIRE FELE DANS F               %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'maillage.inc'
      DIMENSION F(1),FELE(1),IGLO(NBN*NDPN) 
      integer IM

COMM  ......//////// TABLE CORRESPONDANCE DES DDL  ////////......
      DO 5 K=1,NBN
      DO 5 KK=1,NDPN
      IGLO(NDPN*(K-1)+KK)=NDPN*(MAIL(IM,K)-1)+KK
 5    CONTINUE
 
COMM .....////// ASSEMBLAGE DU VECTEUR  ELEMENTAIRE //////.............
      N=NDPN*NBN        
      DO 15 I=1,N
15    F(IGLO(I))=F(IGLO(I))+FELE(I) 
      END
      
            
