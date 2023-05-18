      SUBROUTINE GAU_TRI1(ITYP,GAUSS,POIDS)   !triangulaire 1 pg
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,8),POIDS(8),RQ4(1)
COMM   
      DATA RQ4/0.333333333/
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,1
        DO 1 J=1,1
        DO 1 K=1,1
        N=I+(J-1)*2+(K-1)*4
        GAUSS(1,N)=RQ4(1)
        GAUSS(2,N)=RQ4(1)
        GAUSS(3,N)=0
        POIDS(N)=.5D0
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      END        
      
      
      
      
      
      
      
      
      SUBROUTINE GRIGCO_4pts(ITYP,GAUSS,POIDS)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(2)
COMM   
      DATA RQ4/-0.577350269189626,0.577350269189626/
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,2
        DO 1 J=1,2
        DO 1 K=1,1
        N=I+(J-1)*2+(K-1)*4
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=RQ4(J)
        GAUSS(3,N)=0
        POIDS(N)=1.D0
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      END   
      
      SUBROUTINE GRIGCO_1pts(ITYP,GAUSS,POIDS)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(2)
COMM   
      DATA RQ4/-0,0/
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,1
        DO 1 J=1,1
        DO 1 K=1,1
        N=I+(J-1)*2+(K-1)*4
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=RQ4(J)
        GAUSS(3,N)=0
        POIDS(N)=4.D0
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      END        

      SUBROUTINE FORMQ4(KSI,ETA,N,NKSI,NETA)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  FONCTIONS DE FORME ET LEURS DERIVEES  AU POINT  KSI,ETA        % 
COMM %           -- QUADRANGLE A HUIT NOEUDS : Q4 --                   %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION KSI,ETA,N(*),NKSI(*),NETA(*)
!        write(*,1)KSI,ETA
COMM....................../////////////////////////....................      
      X1PK=1.D0+KSI !
      X1PE=1.D0+ETA
      X1MK=1.D0-KSI
      X1ME=1.D0-ETA

      
COMM ...........//// Fonctions d' interpolation ////////................
      N(1)=0.25D0*X1MK*X1ME
      N(2)=0.25D0*X1PK*X1ME
      N(3)=0.25D0*X1PK*X1PE
      N(4)=0.25D0*X1MK*X1PE
   
COMM .....//// Derivee des fonctions d' interpolation (KSI,ETA) //////....                        
      NKSI(1)=-0.25D0 *X1ME 
      NKSI(2)= 0.25D0 *X1ME 
      NKSI(3)= 0.25D0 *X1PE 
      NKSI(4)=-0.25D0 *X1PE    
                                      
      NETA(1)=-0.25D0 *X1MK 
      NETA(2)=-0.25D0 *X1PK 
      NETA(3)= 0.25D0 *X1PK 
      NETA(4)= 0.25D0 *X1MK 
1     FORMAT(10(G15.8,1X))
      END
COMM ...........................////////////...........................
COMM ...//////  FIN DU MODULE REGROUPANT LES SUBROUTINES DE  //////....
COMM ...//////  CALCULS  DE RIGIDITE TANGENTE                //////....
COMM ...........................////////////...........................




      SUBROUTINE FORMQ8(KSI,ETA,N,NKSI,NETA)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  FONCTIONS DE FORME ET LEURS DERIVEES  AU POINT  KSI,ETA        % 
COMM %           -- QUADRANGLE A HUIT NOEUDS : Q8 --                   %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION KSI,ETA,N(*),NKSI(*),NETA(*)
COMM....................../////////////////////////....................      
      X1PK=1.D0+KSI
      X1PE=1.D0+ETA
      X1MK=1.D0-KSI
      X1ME=1.D0-ETA
      X1MKC=1.D0-KSI*KSI
      X1MEC=1.D0-ETA*ETA
      X2KPE=2.D0*KSI+ETA
      X2KME=2.D0*KSI-ETA
      XKP2E=KSI+2.D0*ETA
      XKM2E=KSI-2.D0*ETA
COMM ...........//// Fonctions d' interpolation ////////................
      N(1)=-0.25D0*X1MK*X1ME*(1.D0+KSI+ETA)
      N(2)=-0.25D0*X1PK*X1ME*(1.D0-KSI+ETA)
      N(3)=-0.25D0*X1PK*X1PE*(1.D0-KSI-ETA)
      N(4)=-0.25D0*X1MK*X1PE*(1.D0+KSI-ETA)
      N(5)= 0.5D0*X1MKC*X1ME
      N(6)= 0.5D0*X1PK*X1MEC
      N(7)= 0.5D0*X1MKC*X1PE
      N(8)= 0.5D0*X1MK*X1MEC      
COMM .....//// Derivee des fonctions d' interpolation (KSI,ETA) //////....                        
      NKSI(1)=0.25D0 *X1ME *X2KPE
      NKSI(2)=0.25D0 *X1ME *X2KME
      NKSI(3)=0.25D0 *X1PE *X2KPE
      NKSI(4)=0.25D0 *X1PE *X2KME      
      NKSI(5)=-X1ME*KSI
      NKSI(6)= 0.5D0*X1MEC      
      NKSI(7)=-X1PE*KSI
      NKSI(8)=-0.5D0*X1MEC                                     
      NETA(1)= 0.25D0 *X1MK *XKP2E
      NETA(2)=-0.25D0 *X1PK *XKM2E
      NETA(3)= 0.25D0 *X1PK *XKP2E
      NETA(4)=-0.25D0 *X1MK *XKM2E
      NETA(5)=-0.5D0*X1MKC
      NETA(6)=-X1PK*ETA      
      NETA(7)= 0.5D0*X1MKC 
      NETA(8)=-X1MK*ETA       
      END
COMM ...........................////////////...........................
COMM ...//////  FIN DU MODULE REGROUPANT LES SUBROUTINES DE  //////....
COMM ...//////  CALCULS  DE RIGIDITE TANGENTE                //////....
COMM ...........................////////////...........................


      SUBROUTINE GRIGCO_2x2(ITYP,GAUSS,POIDS)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(2),Wei(2),GAU(2)
COMM   
      DATA RQ4/-0.577350269189626, 0.577350269189626/
      DATA Wei/1.,1./
      DATA GAU/-0.577350269189626, 0.577350269189626/
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,2
        DO 1 J=1,1
        DO 1 K=1,2
        N=I+(J-1)*2+(K-1)*2
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=GAU(K)
        GAUSS(3,N)=0.
        POIDS(N)=Wei(I)
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      END  
      
      SUBROUTINE GRIGCO_1x1(ITYP,GAUSS,POIDS,NBINT)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(1),Wei(1),GAU(1)
      INTEGER NBINT
COMM   
      DATA RQ4/0. /
      DATA Wei/2./
      DATA GAU/0./
      NBINT=1
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,1
        DO 1 J=1,1
        DO 1 K=1,1

        N=I+(J-1)*2+(K-1)*2       
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=0.
        GAUSS(3,N)=0.
        POIDS(N)=Wei(I)      
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
 11   FORMAT(162(G15.8,1X))       
      END   
      
      SUBROUTINE GRIGCO_2x1(ITYP,GAUSS,POIDS,NBINT)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(2),Wei(2),GAU(2)
      INTEGER NBINT
COMM   
      DATA RQ4/-0.577350269189626, 0.577350269189626/
      DATA Wei/1.,1./
      DATA GAU/-0.577350269189626, 0.577350269189626/
      NBINT=2
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,2
        DO 1 J=1,1
        DO 1 K=1,1

        N=I+(J-1)*2+(K-1)*2       
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=0.
        GAUSS(3,N)=0.
        POIDS(N)=Wei(I)      
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
 11   FORMAT(162(G15.8,1X))       
      END   
      
      
      SUBROUTINE GRIGCO_3x1(ITYP,GAUSS,POIDS,NBINT)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(3),Wei(3),GAU(2)
      INTEGER NBINT
COMM   
      DATA RQ4/-.774596669241,0,.774596669241/
      DATA Wei/.555555555555,.888888888888,.555555555555/
      NBINT=3
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,3
        DO 1 J=1,1
        DO 1 K=1,1

        N=I+(J-1)*2+(K-1)*2       
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=0.
        GAUSS(3,N)=0.
        POIDS(N)=Wei(I)
      
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      
 11   FORMAT(162(G15.8,1X))       
      END  
      
      SUBROUTINE GRIGCO_4x1(ITYP,GAUSS,POIDS,NBINT)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(4),Wei(4),GAU(2)
      INTEGER NBINT
COMM   
      DATA RQ4/-.861136311594053,-.339981043584856,.339981043584856,
     &      .861136311594053/
      DATA Wei/ .347854845137454, .652145154862546,.652145154862546,
     &      .347854845137454/
      NBINT=4
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,4
        DO 1 J=1,1
        DO 1 K=1,1

        N=I+(J-1)*2+(K-1)*2       
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=0.
        GAUSS(3,N)=0.
        POIDS(N)=Wei(I)
      
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      
 11   FORMAT(162(G15.8,1X))       
      END    
     
      SUBROUTINE GRIGCO_8x1(ITYP,GAUSS,POIDS,NBINT)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(8),Wei(8)
      INTEGER NBINT
COMM   
      DATA RQ4/  -.9602898564975362,-.7966664774136268,
     &           -.5255324099163290,-.1834346424956498,  
     &            .1834346424956498, .5255324099163290,
     &            .7966664774136268, .9602898564975362/
     
     
     
      DATA Wei/   .10122853629036970, .22238103445337450,
     &            .31370664587788740, .36268378337836200,
     &            .36268378337836200, .31370664587788740,
     &            .22238103445337450, .10122853629036970/
      NBINT=8
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,8
        DO 1 J=1,1
        DO 1 K=1,1

        N=I+(J-1)*2+(K-1)*2       
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=0.
        GAUSS(3,N)=0.
        POIDS(N)=Wei(I)
      
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      
 11   FORMAT(162(G15.8,1X))       
      END 
      
      SUBROUTINE GRIGCO_16x1(ITYP,GAUSS,POIDS,NBINT)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(16),Wei(16)
      INTEGER NBINT
COMM   
      DATA RQ4/  -.9894009349916499,-.9445750230732326,
     &           -.8656312023878318,-.7554044083550030,
     &           -.6178762444026438,-.4580167776572274,
     &           -.2816035507792589,-.0950125098376374, 
     &            .0950125098376374, .2816035507792589, 
     &            .4580167776572274, .6178762444026438,  
     &            .7554044083550030, .8656312023878318,
     &            .9445750230732326, .9894009349916499/
     
     
     
      DATA Wei/   .02715245941175185, .06225352393864778,
     &            .09515851168249290, .12462897125553390,
     &            .14959598881657330, .16915651939500240,
     &            .18260341504492360, .18945061045506850,
     &            .18945061045506850, .18260341504492360,
     &            .16915651939500240, .14959598881657330,
     &            .12462897125553390, .09515851168249290,
     &            .06225352393864778, .02715245941175185/
      NBINT=16
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,16
        DO 1 J=1,1
        DO 1 K=1,1

        N=I+(J-1)*2+(K-1)*2       
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=0.
        GAUSS(3,N)=0.
        POIDS(N)=Wei(I)
      
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      
 11   FORMAT(162(G15.8,1X))       
      END    
      
      
      SUBROUTINE GRIGCO_5x1(ITYP,GAUSS,POIDS,NBINT)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ4(5),Wei(5),GAU(2)
      INTEGER NBINT
COMM   
      DATA RQ4/-.906179845938664,-.538469310105683,0.,
     &      .538469310105683,.906179845938664/
      DATA Wei/ .236926885056189, .478628670499365,.568888889888889,
     &      .478628670499365,.236926885056189/
     
c       DATA Wei/ .236926885056189, .478628670499367,.501960784313726,
c      &      .478628670499365,.236926885056189/     
     
    
     
      NBINT=5
      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,5
        DO 1 J=1,1
        DO 1 K=1,1

        N=I+(J-1)*2+(K-1)*2       
        GAUSS(1,N)=RQ4(I)
        GAUSS(2,N)=0.
        GAUSS(3,N)=0.
        POIDS(N)=Wei(I)
      
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      
 11   FORMAT(162(G15.8,1X))       
      END         
      
      SUBROUTINE GRIGCO_3x2(ITYP,GAUSS,POIDS)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  EN FONCTION DU TYPE DE L'ELEMENT, ON RENVOIE LES COORDONNEES   %
COMM %  DES POINTS DE GAUSS, ET LES POIDS DE GAUSS    4pts             %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAUSS(3,*),POIDS(*),RQ2(2),RQ3(3),Wei(3),GAU(2)
COMM   
      DATA RQ3/-.774596669241,0,.774596669241/
      DATA Wei/.555555555555,.888888888888,.555555555555/
      DATA RQ2/-0.577350269189626, 0.577350269189626/

      IF (ITYP.EQ.41.OR.ITYP.EQ.83) THEN
        DO 1 I=1,2
        DO 1 J=1,3
        DO 1 K=1,1

        N=J+(I-1)*3+(K-1)*2       
        GAUSS(1,N)=RQ3(J)
        GAUSS(2,N)=RQ2(I)
        GAUSS(3,N)=0.
        POIDS(N)=Wei(J)
      
 1      CONTINUE  
        ELSE
       STOP ' PB :GRIGCOD  '
      ENDIF
      
 11   FORMAT(162(G15.8,1X))       
      END         
      
      SUBROUTINE FORMS3(KSI,N,NK,NKK,le) !elet poutre cubique Hermite
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  FONCTIONS DE FORME ET LEURS DERIVEES  AU POINT  KSI,ETA        % 
COMM %           -- QUADRANGLE A HUIT NOEUDS : Q4 --                   %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION KSI,N(*),NK(*),NKK(*),le,leS8,X1PK,X1MK,X2PK,
     & X2MK,X1PKC,X1MKC,DKSI
COMM....................../////////////////////////....................      
      
 
      X1PK=1.D0+KSI ! 
      X1MK=1.D0-KSI
      X2PK=2.D0+KSI
      X2MK=2.D0-KSI
      X1PKC=(1.D0+KSI)**2
      X1MKC=(1.D0-KSI)**2
      DKSI=2*KSI
      leS8=le/8.

      
COMM ...........//// Fonctions d' interpolation ////////................
      N(1)= 0.25D0*X1MKC*X2PK
      N(2)= leS8  *X1MKC*X1PK
      N(3)= 0.25D0*X1PKC*X2MK
      N(4)=-leS8  *X1PKC*X1MK
      

      
   
COMM .....//// Derivee des fonctions d' interpolation (KSI,ETA) //////....                         
                                      
      NK(1)= 0.25D0*X1MKC-0.5*X1MK*X2PK
      NK(2)=leS8*X1MKC-2*leS8*X1MK*X1PK  
      NK(3)= 0.5*X1PK*X2MK  -0.25*X1PKC
      NK(4)=leS8*X1PKC-2*leS8*X1PK*X1MK
      
   
      
      NKK(1)=0.5*X2PK+KSI-1
      NKK(2)=2.*leS8*X1PK-4.*leS8*X1MK
      NKK(3)=-KSI+0.5*X2MK-1 
      NKK(4)=4.*leS8*X1PK-2.*leS8 *X1MK   
         
1     FORMAT(162(G15.8,1X))     
      END
COMM ...........................////////////...........................
COMM ...//////  FIN DU MODULE REGROUPANT LES SUBROUTINES DE  //////....
COMM ...//////  CALCULS  DE RIGIDITE TANGENTE                //////....
COMM ...........................////////////...........................

      SUBROUTINE FORMS1(KSI,N,NKSI) !elet poutre lineair
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  FONCTIONS DE FORME ET LEURS DERIVEES  AU POINT  KSI,ETA        % 
COMM %           -- QUADRANGLE A HUIT NOEUDS : Q4 --                   %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION KSI,N(*),NKSI(*),le
COMM....................../////////////////////////....................      
      X1PK=1.D0+KSI !
      X1MK=1.D0-KSI
      X2PK=2.D0+KSI
      X2MK=2.D0-KSI
      X1PKC=1.D0+KSI*KSI
      X1MKC=1.D0-KSI*KSI
      DKSI=2*KSI
      leS8=le/8


      
COMM ...........//// Fonctions d' interpolation ////////................
      N(1)=0.5D0*X1MK
      N(2)=0.5D0*X1PK
   
COMM .....//// Derivee des fonctions d' interpolation (KSI,ETA) //////....                         
                                      
      NKSI(1)=-0.5D0   
      NKSI(2)= 0.5D0    
   
      END


      
      SUBROUTINE FORMS2(KSI,N,NKSI) !elet poutre quadratique
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %  FONCTIONS DE FORME ET LEURS DERIVEES  AU POINT  KSI,ETA        % 
COMM %           -- QUADRANGLE A HUIT NOEUDS : Q4 --                   %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION KSI,N(*),NKSI(*),le
COMM....................../////////////////////////....................      
      X1PK=1.D0+KSI !
      X1MK=1.D0-KSI
      X2PK=2.D0+KSI
      X2MK=2.D0-KSI
      X1PKC=1.D0+KSI*KSI
      X1MKC=1.D0-KSI*KSI
      DKSI=2*KSI
      
COMM ...........//// Fonctions d' interpolation ////////................
      N(1)=-0.5D0*X1MK*KSI
      N(2)=0.5D0*X1PK*KSI
      N(3)=X1MK*X1PK
   
COMM .....//// Derivee des fonctions d' interpolation (KSI,ETA) //////....                         
                                      
      NKSI(1)=KSI-0.5D0   
      NKSI(2)=KSI+0.5D0  
      NKSI(3)=-DKSI
   
      END      
