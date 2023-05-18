      SUBROUTINE MATVEC(A,B,C,N1,N2)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %  EFFECTUE LE PRODUIT  C=A.B   A(N1,N2) ; B(N2) ; C(N1)         %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(N1,N2),B(N2),C(N1)
      integer N1,N2
      DO  I=1,N1
      C(I)=0.D0
      DO  J=1,N2
      C(I)=C(I)+A(I,J)*B(J)
 2    ENDDO
 1    ENDDO      
      END
      
      SUBROUTINE TMATVEC(A,B,C,N1,N2)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                         T                                      %
COMM %  EFFECTUE LE PRODUIT  C= A.B   A(N1,N2) ; B(N1) ; C(N2)        %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(N1,N2),B(N1),C(N2)
      DO 1 I=1,N2
      C(I)=0.D0
      DO 1 J=1,N1
 1    C(I)=C(I)+A(J,I)*B(J)      
      END
      
      SUBROUTINE MATMAT(A,B,C,N1,N2,N3)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %  EFFECTUE LE PRODUIT  C=A.B   A(N1,N2) ; B(N2,N3) ; C(N1,N3)   %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(N1,N2),B(N2,N3),C(N1,N3)
      INTEGER I,J,K,N1,N2,N3
      C=0
      DO 1 I=1,N1
      DO 1 J=1,N3
      C(I,J)=0.D0
      DO 1 K=1,N2
 1    C(I,J)=C(I,J)+A(I,K)*B(K,J)      
      END
      
      SUBROUTINE TMATMAT(A,B,C,N1,N2,N3)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %  EFFECTUE LE PRODUIT  C=TA.B   A(N1,N2) ; B(N2,N3) ; C(N1,N3)   %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(N1,N2),B(N2,N3),C(N1,N3)
      INTEGER I,J,K,N1,N2,N3
      C=0
      DO 1 I=1,N1
      DO 1 J=1,N3
      C(I,J)=0.D0
      DO 1 K=1,N2
 1    C(I,J)=C(I,J)+A(K,I)*B(K,J)      
      END      
        
      SUBROUTINE TVECVEC(A,B,C,N1,N2)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %  EFFECTUE LE PRODUIT  C=TA.B   A(N1) ; B(N2) ; C(N1,N2)   %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(1),B(1),C(N1,1)
      INTEGER I,J,K,N1,N2
      DO 1 I=1,N1
      DO 1 J=1,N2
 1    C(I,J)=A(I)*B(J)      
      END  
      
      SUBROUTINE TVECVEC_add(A,B,C,N1,N2)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %  EFFECTUE LE PRODUIT  C=TA.B   A(N1) ; B(N2) ; C(N1,N2)   %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(*),B(*),C(N1,N2)
      INTEGER I,J,K,N1,N2
      DO 1 I=1,N1
      DO 1 J=1,N2
 1    C(I,J)=C(I,J)+A(I)*B(J)      
      END      
      
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %                EVE   Version 2.1  (Novembre 94)                 %
COMM %                                                                 %
COMM %       //////  MODULE REGROUPANT LES SUBROUTINES   //////        %
COMM %       //////             UTILITAIRES              //////        %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      SUBROUTINE TBDB(B,D,DB,ARGL,ND,NDDLE,ITS)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %    CALCUL LE PRODUIT TB*D*B ---> STOCKE DANS ARGL              %
COMM %     ARGL : MATRICE TRIANGULAIRE SUPERIEURE STOCKEE PAR LIGNE   %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARGL(ITS),B(ND,NDDLE),D(ND,ND),DB(ND,NDDLE)
COMM
      DO 1 I=1,ND  
      DO 1 J=1,NDDLE
      DB(I,J)=0.D0
      DO 1 K=1,ND
 1    DB(I,J)=DB(I,J)+D(I,K)*B(K,J)
 
      L=1
      DO 6 I=1,NDDLE 
      DO 7 J=I,NDDLE
      ARGL(L)=0.D0
        DO 8 K=1,ND
        ARGL(L)=ARGL(L)+B(K,I)*DB(K,J)
 8      CONTINUE
        L=L+1
 7    CONTINUE
 6    CONTINUE
      END 
      
      SUBROUTINE TNN(N,ARGL,ND,NDDLE,ITS)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %    CALCUL LE PRODUIT TN*N   ---> STOCKE DANS ARGL              %
COMM %     ARGL : MATRICE TRIANGULAIRE SUPERIEURE STOCKEE PAR LIGNE   %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION ARGL(ITS),N(ND,NDDLE)
COMM
      L=1
      DO 6 I=1,NDDLE 
      DO 7 J=I,NDDLE
      ARGL(L)=0.D0
        DO 8 K=1,ND
        ARGL(L)=ARGL(L)+N(K,I)*N(K,J)
 8      CONTINUE
        L=L+1
 7    CONTINUE
 6    CONTINUE
      END 

      
      SUBROUTINE MATLIG(A,B,N)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %  ON STOCKE LA MATRICE CARRE SYMETRIQUE A(N,N) SOUS FORME       % 
COMM %  DE LIGNE B( )                                                 %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(N,N),B(1)
      K=1      
      DO 1 I=1,N
      DO 1 J=I,N
      B(K)=A(I,J)
1     K=K+1
      END
      
      SUBROUTINE LIGMAT(B,A,N)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %  ON REFORME LA MATRICE CARRE SYMETRIQUE A(N,N) A PARTIR        % 
COMM %  DU VECTEUR B QUI CONTIENT SES LIGNES                          %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(N,N),B(1)
      K=1      
      DO 1 I=1,N
      DO 1 J=I,N
      A(I,J)=B(K)
      A(J,I)=B(K)
1     K=K+1
      END 
                          
      SUBROUTINE DSOMME(A,B,N)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %                 EFFECTUE LE CALCUL  A=A+B                      %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(1),B(1)
      DO 1 I=1,N
 1    A(I)=A(I)+B(I)
      END
      
      SUBROUTINE DDIFFE(A,B,N)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %                 EFFECTUE LE CALCUL  A=A-B                      %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(1),B(1)
      DO 1 I=1,N
 1    A(I)=A(I)-B(I)
      END
      
      SUBROUTINE DINCRE(A,B,N,C)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %                 EFFECTUE LE CALCUL  B(I) = A(I)*C              %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(1),B(1),C
      DO 1 I=1,N
 1    B(I)=A(I)*C
      END
      
      SUBROUTINE DACTUA(A,B,N,C)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %                 EFFECTUE LE CALCUL  A(I) = A(I) + B(I)*C       %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(1),B(1),C
      INTEGER N
      DO 1 I=1,N
 1    A(I)=A(I)+B(I)*C
      END
            
      SUBROUTINE DZERO(A,N)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %                MISE A ZERO DE A(I),I=1,N                       %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(*)
      DO 1 I=1,N
 1    A(I)=0.D0
      END
      
      SUBROUTINE ZERO(M,k,l)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION M(k,l) 
        do 10 i=1,k
        do 10 j=1,l
10      M(i,j)=0.D0 
      end      
      
      SUBROUTINE COMPOS(V,IN,ID,NDPN,COMP)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %   On recupere la Ieme composante du vecteur V                  % 
COMM %   IN --> Numero du noeud   ID --> de 1 a 6 numero du d.d.l.    %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION V(1),COMP
      COMP=V((IN-1)*NDPN+ID)
      END 
                  
      FUNCTION CALMAX(X,N)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %     CALCUL DE LA COMPOSANTE MAXIMUN DE ABS(X(I)) ,I=1,N        %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION X(N),CALMAX,R
      CALMAX=0.D0
      DO 1 I=1,N
      R=DABS(X(I))
      IF(R.GT.CALMAX) CALMAX=R
 1    CONTINUE
      END
      
      SUBROUTINE NORMER_F(F,NDL,X)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 % 
COMM %    ON CALCULE LA NORME EUCLIDIENNE DE F                         %
COMM %                                                                 %    
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION F(*)
      X=DSQRT(PROSCA(F,F,NDL))
      END
            
    
      subroutine prod_scal(U,V,N,PROSCA)
      implicit none
      integer N,I
      real*8 U(1),V(1),PROSCA
      PROSCA=0.D0
      DO 1 I=1,N
 1    PROSCA=PROSCA+U(I)*V(I)      
      
      
      end

    
      DOUBLE PRECISION FUNCTION PROSCA(U,V,N)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %      CALCUL DU PRODUIT SCALAIRE DES VECTEURS U(N),V(N)         %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION U(*),V(*)
      PROSCA=0.D0
      DO 1 I=1,N
 1    PROSCA=PROSCA+U(I)*V(I)
      END

      DOUBLE PRECISION FUNCTION PROSCASI(U,V,N)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %      CALCUL DU PRODUIT SCALAIRE DES VECTEURS U(N),V(N)         %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(N),V(N)

      DO 1 I=1,N
 1    PROSCASI=PROSCASI+U(I)*V(I)
      END

      SUBROUTINE SHIFTD(X,Y,LX)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %       ON MET LE TABLEAU X DANS Y SUR UNE LONGEUR LX            %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION X(*),Y(*)
      DO 1 IA=1,LX
    1 Y(IA)=X(IA)
      END
      
      SUBROUTINE MATVECSI(A,B,C,N1,N2)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %  EFFECTUE LE PRODUIT  C=C+A.B   A(N1,N2) ; B(N2) ; C(N1)       %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(N1,N2),B(N2),C(N1)
      DO 1 I=1,N1
      DO 1 J=1,N2
 1    C(I)=C(I)+A(I,J)*B(J)      
      END  
      SUBROUTINE TMATVECSI(A,B,C,N1,N2)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                            T                                   %
COMM %  EFFECTUE LE PRODUIT  C=C+  A.B   A(N1,N2) ; B(N1) ; C(N2)     %  
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(N1,N2),B(N1),C(N2)
      DO 1 I=1,N2
      DO 1 J=1,N1
 1    C(I)=C(I)+A(J,I)*B(J)      
      END    
      
      SUBROUTINE TMATVEC_add(A,B,C,N1,N2,COEF)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                              T                                      %
COMM %  EFFECTUE LE PRODUIT  C = C + A.B * COEF   A(N1,N2) ; B(N1) ; C(N2) %  
COMM %                                                                     %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(N1,N2),B(N1),C(N2),COEF
      DO I=1,N2
c      
      C_inter=0.D0
      DO J=1,N1
      C_inter=C_inter+A(J,I)*B(J)
      ENDDO
c      
      C(I)=C(I)+C_inter*COEF 
      ENDDO   
      END
      
      SUBROUTINE TB_A_add(B,A,C,n,l,m,COEF)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                              T                                            %
COMM %  EFFECTUE LE PRODUIT  C = C + A.B * COEF   A(N2,N1) ; B(N2,N3) ; C(N1,N3) %  
COMM %                                                                           %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION A(l,m),B(l,n),C(n,m),COEF
      DO  I=1,n
       DO  J=1,m
        C_inter=0.D0
        DO  K=1,l
         C_inter=C_inter + B(K,I)*A(K,J)
        ENDDO
        C(I,J)=C(I,J) + COEF * C_inter
       ENDDO
      ENDDO
      END 
      
     
      

      
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %      produit de matrices C = C + TB.D.A                         %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      SUBROUTINE TB_D_A_add(B,D,A,n,l,m,C,COEF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION B(l,n),D(l,l),A(l,m),C(n,m),T(l,m)
      DO  I=1,l
       DO  J=1,m
        T(I,J)=0.D0
        DO  K=1,l
1        T(I,J)=T(I,J)+D(I,K)*A(K,J)
        ENDDO
       ENDDO
      ENDDO

      DO  I=1,n
       DO  J=1,m
        C_inter=0.D0
        DO K=1,l
3        C_inter=C_inter + B(K,I)*T(K,J)
        ENDDO
2       C(I,J)=C(I,J) + COEF * C_inter
       ENDDO
      ENDDO
      
      END  
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %      produit de matrices C = C + TB.D.A                         %
COMM %         resultat : matrice symmetrique                          %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      SUBROUTINE TB_D_A_addsym(B,D,A,n,l,m,C,COEF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION B(l,n),D(l,l),A(l,m),C(n,m),T(l,m)
      DO 1 I=1,l
      DO 1 J=1,m
      T(I,J)=0.D0
      DO 1 K=1,l
1     T(I,J)=T(I,J)+D(I,K)*A(K,J)
c--------
      DO 2 I=1,n
      DO 2 J=I,m
      C_inter=0.D0
      DO 3 K=1,l
3     C_inter=C_inter + B(K,I)*T(K,J)
c--------
      C(I,J)=C(I,J) + COEF * C_inter
2     C(J,I)=C(I,J)
      END

 
      
      DOUBLE PRECISION FUNCTION SIGNE(a)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %      CALCUL DU PRODUIT SCALAIRE DES VECTEURS U(N),V(N)         %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SIGNE=-1.D0
      IF (a.NE.0.D0) THEN
      SIGNE=a/DABS(a)
      ENDIF

      END      
      
      
COMM .........................////////////.............................
COMM ............////// FIN DU MODULE UTILITAIRE //////................
COMM .........................////////////.............................

COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %                 Ecriteure d'un tableau				            %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine EcrireT(A,N)
	DOUBLE PRECISION A(1)
	do 1 i=1,N
	write(*,*)i,A(i)
1	continue
	end
	
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %                 Ecriteure d'une matrice				            %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine EcrireM(A,N1,N2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DOUBLE PRECISION A(N1,1)
2     FORMAT(500(G15.8,1X))
	  do 1 i=1,N1
	  write(19,2)(A(i,j),j=1,N2)
1	  continue
	  end	

c=======================================================================
      subroutine exec_command( cmd, input, output, ier )
C Execute cmd, avec comme entree le fichier input, et sortie le fichier output
C input et output peuvent etre vide 
      implicit none

      character*(*)  cmd, input, output

      character*200  cmdstring
      integer*4 system
      integer ier

      cmdstring = cmd

      if (input.ne.' ') then
         cmdstring = cmdstring(1:len_trim(cmdstring)) // '<' //input
      endif

      if (output.ne.' ') then
         cmdstring = cmdstring(1:len_trim(cmdstring)) // '>' // output
      endif

C      write(*,*) 'cmd=', cmd
C      write(*,*) 'input=', input
C      write(*,*) 'output=', output
C      write(*,*) 'Debug:', cmdstring(1:len_trim(cmdstring))
      ier = system(cmdstring)
      
C      write(*,*) 'result=', ier
      end
	  
COMM .........................////////////.............................
COMM ............////// FIN DU MODULE UTILITAIRE //////................
COMM .........................////////////.............................


      Subroutine conv_tab(T1,Tc,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION T1(N),Tc(N)
      INTEGER N
      DO i=1,N
      Tc(i)=T1(i)
      ENDDO
      END	  
CBR
CBR
CBR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CBR    Methode d'elimination de Gauss, matrice carre
C  
      SUBROUTINE MATINV(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N),IPIVOT(N),INDEX(N,3),PIVOT(N)

c********************
               DO 20 J=1,N
  20           IPIVOT(J)=0
               DO 550 I=1,N
C*************************
C  LA RECHERCHE DU PIVOT
C*************************
               AMAX=0.
               DO 105 J=1,N
               IF(IPIVOT(J)-1)60,105,60
  60           DO 100 K=1,N
               IF(IPIVOT(K)-1)80,100,740
  80           IF(ABS(AMAX)-ABS(A(J,K)))85,100,100
  85           IROW=J
               ICOLUM=K
               AMAX=A(J,K)
  100          CONTINUE
  105          CONTINUE
               IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1
C**************************************************************
C  CHANGEMENT DES LIGNES POUR METTRE LE PIVOT DANS LA DIAGONALE
C**************************************************************
               IF(IROW-ICOLUM)140,260,140
c***************
C140      DETERM=-DETERM
C****************
  140           CONTINUE
C                SWAP=SM(IROW)
C                SM(IROW)=SM(ICOLUM)
C                SM(ICOLUM)=SWAP
                DO 200 L=1,N
                SWAP=A(IROW,L)
                A(IROW,L)=A(ICOLUM,L)
  200           A(ICOLUM,L)=SWAP
  260           INDEX(I,1)=IROW
                INDEX(I,2)=ICOLUM
                PIVOT(I)=A(ICOLUM,ICOLUM)
c*****************************************
C              DETERM=DETERM*PIVOT(I)
c  DIVISION DE LA LIGNE PIVOT PAR LE PIVOT
C*****************************************
                A(ICOLUM,ICOLUM)=1.
                DO 350 L=1,N
                A(ICOLUM,L)=A(ICOLUM,L)/PIVOT(I)
 350  CONTINUE
c****************************
c REDUCTION DE LA LIGNE PIVOT
C****************************
                DO 551 L1=1,N
                IF(L1-ICOLUM)400,551,400
  400           T=A(L1,ICOLUM)
                A(L1,ICOLUM)=0.
                DO 450 L=1,N
  450           A(L1,L)=A(L1,L)-A(ICOLUM,L)*T
  551           CONTINUE
  550           CONTINUE
C******************************
C  CHANGEMENT DE COLUMNS
C*****************************
                DO 710 I=1,N
                L=N+1-I
                IF(INDEX(L,1)-INDEX(L,2))630,710,630
  630           JROW=INDEX(L,1)
                JCOLUM=INDEX(L,2)
                DO 705 K=1,N
                SWAP=A(K,JROW)
                A(K,JROW)=A(K,JCOLUM)
                A(K,JCOLUM)=SWAP
  705           CONTINUE
  710           CONTINUE
  740           RETURN
                END
                
                
                
       subroutine inverser(A,n)
       implicit double precision (A-H,O-Z)
       double precision A(n,n)
       integer n
       real*8 b(n,n)
       
!        dimension A(10000,10000)
       do 1 i=1,n
       do 2 j=1,n
2      b(i,j)=0.D0
       b(i,i)=1.
1     continue 
       do 3 k=1,n
         DP=1./A(k,k)
         do 4 j=1,n
         A(k,j)=A(k,j)*DP
         b(k,j)=b(k,j)*DP
4        continue 
        do 5 j=1,n
        IF (j.NE.k) THEN
         D=A(j,k)
         do 6 L=1,n
         A(j,L)=A(j,L)-A(k,L)*D
         b(j,L)=b(j,L)-b(k,L)*D 
6        continue 
        ELSE
        ENDIF          
5       continue 
3     continue 
      do 7 i=1,n
      do 7 j=1,n
7      A(i,j)=b(i,j)
      END            
c --------------	  
	  
	  
      REAL FUNCTION FindDet(matrix, n)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, DIMENSION(n,n) :: matrix
     
      REAL :: m, temp
      INTEGER :: i, j, k, l
      LOGICAL :: DetExists = .TRUE.
      l = 1
	!Convert to upper triangular form
      DO k = 1, n-1
      IF (matrix(k,k) == 0) THEN
			DetExists = .FALSE.
			DO i = k+1, n
				IF (matrix(i,k) /= 0) THEN
					DO j = 1, n
						temp = matrix(i,j)
						matrix(i,j)= matrix(k,j)
						matrix(k,j) = temp
					END DO
					DetExists = .TRUE.
					l=-l
					EXIT
				ENDIF
			END DO
			IF (DetExists .EQV. .FALSE.) THEN
				FindDet = 0
				return
			END IF
		ENDIF
		DO j = k+1, n
			m = matrix(j,k)/matrix(k,k)
			DO i = k+1, n
				matrix(j,i) = matrix(j,i) - m*matrix(k,i)
			END DO
		END DO
      END DO
	
	!Calculate determinant by finding product of diagonal elements
      FindDet = l
      DO i = 1, n
		FindDet = FindDet * matrix(i,i)
      END DO
	
      END FUNCTION FindDet
      
      
      
      
