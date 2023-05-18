      
      subroutine SERIES(X0u,X0v,Xiu,Xiv,Fnli,Ci,ndlu,ndlv,ndl)
      implicit none
      include 'donnees.inc'
      include 'maillage.inc'
c       include 'varstok.inc'
      real*8 MTG(1000*1000),Fext(1000),Fnl(1000),Xi(1000),X1(1000),
     & AUX(1000*1000),X0u(*),X0v(*),Xiu(ndlu,*),Xiv(ndlv,*),Fnli(ndl,*),
     & Ci(*)
      integer nintu,nintv,ndpnu,ndpnv,ndlu,ndlv,ndlue,ndlve,nelu,nelv,
     & NLIU,NLIV,LIV(2,1000),VLIV(1000),LIU(2,1000),VLIU(1000),nel,inc,
     & i,j,k,l,iordre,ndl,deco,Ichoix
      
      
      ndpnu=1
      ndpnv=2
      nintu=2
      nintv=3
      nelu=2
      nelv=2
      nel=2
      ndlue=ndpnu*nelu
      ndlve=ndpnv*nelv
 
      NLIU=0
      NLIV=0
      MTG=0.
      Fnl=0.
      
      
      do i=1,NLIAI
       
       if    (LIAI(2,i).eq.1 )then
        NLIU=NLIU+1
        LIU(1,NLIU)=LIAI(1,i)
        LIU(2,NLIU)=LIAI(2,i)
        VLIU(NLIU)=LIAIV(i)
        
       elseif(LIAI(2,i).eq.2 .or. LIAI(2,i).eq.3)then
        NLIV=NLIV+1
        LIV(1,NLIV)=LIAI(1,i)
        LIV(2,NLIV)=LIAI(2,i)-1   
        VLIV(NLIV)=LIAIV(i)
       endif
      enddo !i=1,NLIAI
  
      call RTG(MTG,X0u,X0v,nel,ndpnu,ndpnv,ndlu,ndlv,ndl) 
      call FORCE(Fext,ndpnu,ndpnv,ndlu,ndlv)     
      call CALIMPNNL(MTG,LIU,VLIU,Fext,NLIU,ndpnu,0,ndl)
      call CLTERMUNIT(MTG,LIU,NLIU,ndpnu,0,ndl)
      call CALIMPNNL(MTG,LIV,VLIV,Fext,NLIV,ndpnv,ndlu,ndl)
      call CLTERMUNIT(MTG,LIV,NLIV,ndpnv,ndlu,ndl)  
      
      CALL SHIFTD(Fext,X1,ndl)
c       call EcrireM(MTG,ndl,ndl)
c       stop 'ffff'
      
      deco=1
      AUX=MTG
      call RESOL_LU(MTG,X1,ndl,deco)
      
      
      CALL SHIFTD(X1(1),Xiu(1,1),ndlu)
      CALL SHIFTD(X1(1+ndlu),Xiv(1,1),ndlv)    
   
      CALL CVORD1(Ci(1),X1,ndl,pilot) 
       
      do iordre=2,NORDRE
       call FNLG(Fnl,Xiu,Xiv,X0u,X0v,nel,ndpnu,ndpnv,iordre,ndlu,ndlv)
       CALL SHIFTD(Fnl,Fnli(1,iordre),ndl)       
       call CALIMP0(LIU,Fnl,NLIU,ndpnu,0,ndl)
       call CALIMP0(LIV,Fnl,NLIV,ndpnv,ndlu,ndl)
       CALL SHIFTD(Fnl,Xi,ndl)
       
       deco=1
       MTG=AUX
       call RESOL_LU(MTG,Xi,ndl,deco) !(MTG,X1,ndl,deco)
       
       CALL SHIFTD(Xi(1)     ,Xiu(1,iordre),ndlu)
       CALL SHIFTD(Xi(1+ndlu),Xiv(1,iordre),ndlv)  
       
       CALL CVORDP(Ci(1),X1,Xi,ndl,iordre,pilot) 
      
       
      enddo  !iordre=2,NORDRE
      
      
c 	   do iordre=1,NORDRE
c 	    write(*,*)iordre,Ci(iordre)
c 	   enddo  
c 	   pause 'hhhyy'      
      
      
2     FORMAT(500(G15.8,1X))
1     FORMAT(162(G15.8,1X))        
4     FORMAT(1X,A11,I3,1X,G10.3,2(G15.6)) 
3     FORMAT(1X,A7,G10.3)       
      end
      
      
      SUBROUTINE CVORD1(C,V1,NDL,Ichoix)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 % 
COMM %  ON CALCULE LE COEFFICIENT  C(1) ET LE VECTEUR V(1)             %
COMM %                                                                 %    
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(*),V1(*)
COMM  il faut eliminer XL   Ichoix=1 :force  2: deplacement


      IF (Ichoix.EQ.1)C(1)=1.D0  
      IF (Ichoix.EQ.2)C(1)=DSQRT(1.D0/(1.D0+PROSCA(V1(1),V1(1),NDL)))
      IF (Ichoix.EQ.3)C(1)=DSQRT(1.D0/(PROSCA(V1(1),V1(1),NDL))) 

      
      
      DO 1 I=1,NDL
 1    V1(I)=V1(I)*C(1)
 
   
      END      
      
      SUBROUTINE CVORDP(C,V1,Vi,NDL,IORDRE,Ichoix)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 % 
COMM %  ON CALCULE LE COEFFICIENT  C(p) ET ON TERMINE LE VECTEUR V(p)  %
COMM %                                                                 %    
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(*),V1(*),Vi(*)
COMM
COMM
      IF (Ichoix.EQ.1) C(IORDRE)=0.D0 
      IF (Ichoix.EQ.2) C(IORDRE)=-PROSCA(Vi(1),V1(1),NDL)*C(1)
      IF (Ichoix.EQ.3) C(IORDRE)=-PROSCA(Vi(1),V1(1),NDL)*C(1)

COMM  
COMM ..../////   calcul de V(p) //////...................................
COMM	
     
      DO 20 I=1,NDL
 20   Vi(I)=Vi(I)+V1(I)*C(IORDRE)/C(1)
      END       
