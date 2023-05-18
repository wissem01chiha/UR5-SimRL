      SUBROUTINE ASSMATNXN(MATEL,MAT,IM,NBN,NDPN,NDL)
!                          

COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %     ASSEMBLAGE DE LA MATRICE ELEMENTAIRE RGL DANS RIG           %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'maillage.inc'
      DOUBLE PRECISION MATEL(NBN*NDPN,NBN*NDPN),MAT(NDL,1)
      integer I,J,Jdln,Jnoe,Inoe,Idln,IM,NBN,NDPN
      
      Jnoe=0
      Jdln=0   
  
 

      DO J=1,NBN*NDPN 
      
       if(mod(J-1,NDPN).eq.0)Jnoe=Jnoe+1
       if(mod(J-1,NDPN).eq.0)Jdln=0
       Jdln=Jdln+1
      

       
       Inoe=0
       Idln=0       
       DO  I=1,NBN*NDPN
        if(mod(I-1,NDPN).eq.0)Inoe=Inoe+1
        if(mod(I-1,NDPN).eq.0)Idln=0 
c         write(*,*)I,J,Inoe,Idln,Jnoe,Jdln
        
        Idln=Idln+1       

        JG=NDPN*(MAIL(IM,Jnoe)-1)+Jdln
        IG=NDPN*(MAIL(IM,Inoe)-1)+Idln
c         write(*,*)IM,IG,JG
        
        MAT(IG,JG)=MAT(IG,JG)+MATEL(I,J)
 
       ENDDO      
      ENDDO
  
          
          
 1     FORMAT(162(G15.8,1X)) 
      END  
      
      
      SUBROUTINE ASSMATNXM(MATEL,MAT,IM,NBN,NDPN1,NDPN2,N1,N2)
!                          

COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %     ASSEMBLAGE DE LA MATRICE ELEMENTAIRE RGL DANS RIG           %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'maillage.inc'
      DOUBLE PRECISION MATEL(NBN*NDPN1,NBN*NDPN2),MAT(N1,1)
      integer I,J,Jdln,Jnoe,Inoe,Idln,IM,NBN,NDPN
      
      Jnoe=0
      Jdln=0          
      DO J=1,NBN*NDPN2 
      
       if(mod(J-1,NDPN2).eq.0)Jnoe=Jnoe+1
       if(mod(J-1,NDPN2).eq.0)Jdln=0
       Jdln=Jdln+1
      

       
       Inoe=0
       Idln=0       
       DO  I=1,NBN*NDPN1
        if(mod(I-1,NDPN1).eq.0)Inoe=Inoe+1
        if(mod(I-1,NDPN1).eq.0)Idln=0 
c         write(*,*)I,J,Inoe,Idln,Jnoe,Jdln
        
        Idln=Idln+1       

        JG=NDPN2*(MAIL(IM,Jnoe)-1)+Jdln
        IG=NDPN1*(MAIL(IM,Inoe)-1)+Idln
c         write(*,*)IM,IG,JG
        
        MAT(IG,JG)=MAT(IG,JG)+MATEL(I,J)
 
       ENDDO      
      ENDDO
c       pause 'ggg'
          
          
 1     FORMAT(162(G15.8,1X)) 
      END       
      
      SUBROUTINE ASSKG(MTG,MUU,MVV,MUV,NDLU,NDLV,NDL)
!                          

COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %     ASSEMBLAGE DE LA MATRICE ELEMENTAIRE RGL DANS RIG           %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'maillage.inc'
      INTEGER I,J,NDL 
      DOUBLE PRECISION MUU(NDLU,1),MVV(NDLV,1),MUV(NDLU,1),
     & MTG(NDL,1),FGU(NDL),
     & FGV(NDL)
      
      
      MTG=0.
      DO I=1,NDL
       DO J=1,NDL
       
        IF    (J.le.NDLU.and.I.le.NDLU)then
         MTG(I,J)=MUU(I,J)
        elseif(J.gt.NDLU.and.I.le.NDLU)then
         MTG(I,J)=MUV(I,J-NDLU)
        elseif(J.le.NDLU.and.I.gt.NDLU)then
         MTG(I,J)=MUV(J,I-NDLU)       
        elseif(J.gt.NDLU.and.I.gt.NDLU)then
         MTG(I,J)=MVV(I-NDLU,J-NDLU)       
        endif
       ENDDO
       
       
      ENDDO
c       DO I=1,NDL
c       write(*,*)I
c        write(19,2)(MTG(I,J),J=1,NDL)
c       ENDDO
      
      

      
          
 1     FORMAT(162(G15.8,1X)) 
 2     FORMAT(1000(G15.8,1X)) 
      END        
      
