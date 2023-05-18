      subroutine increm()
      
      implicit none
      include 'donnees.inc'
      include 'maillage.inc'
      include 'varstok.inc'
      
      real*8 FNLNORM(50),RAYON,X1(1000),XN(1000),X0(1000),NrmX1,NrmXN,
     & DA,APOSITIF,LAMBDA,Xu(1000),Xv(1000),A,CMAX,APOS,CTARGET,ALEFT,
     & ARIGHT,X(1000),Xray(50*1000)
      integer i,j,k,l,IP,IFNL,ndpnu,ndpnv,ndlu,ndlv,ndl,NPOI,iordre,
     & Ichoix,IPOI
      
 43   FORMAT( 1X,10X,'FNL',1X,I2,1X,' = ',G15.8 )
 54   FORMAT(/2X,'LONGUEUR DU PAS MAN = ',G15.8,2X,'LAMBDA =',
     &     G15.8 /) 
 1    FORMAT(162(G15.8,1X))
 
 
               ndpnu=1
               ndpnv=2
      ndlu=ndpnu*nbnoe
      ndlv=ndpnv*nbnoe
         ndl=ndlu+ndlv
                  IP=1
                  C0=0.
                 X0u=0.
                 X0v=0.
                   X=0.
                  X1=0.
                  Ni=0.
                  N0=0.
                 Xiu=0.
                 Xiv=0.
                  Xu=0.
                  Xv=0.
              LAMBDA=0.
      OPEN(14,FILE='courbeMAN',STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      do while(IP.le.nbincr) 
      WRITE(*,'(   )')
      WRITE(*,'(A14,I4,A10)') '=========|PAS',IP,'|========='
      WRITE(*,'(   )')
      call DZERO(Fnli,ndl*NORDRE)
      CALL SERIES(X0u,X0v,Xiu,Xiv,Fnli,Ci,ndlu,ndlv,ndl)
 
       
       FNLNORM=0.
       DO IFNL=2,NORDRE
       call NORMER_F(Fnli((IFNL-1)*ndl+1),ndl,FNLNORM(IFNL))
       WRITE(6,43) IFNL,FNLNORM(IFNL)
       
c        
       ENDDO
       WRITE(6,*)
       
c ------calcul rayon   

c -------   
       CALL SHIFTD(Xiu(1),X1(1),ndlu)
       CALL SHIFTD(Xiv(1),X1(ndlu+1),ndlv)
       CALL SHIFTD(Xiu((NORDRE-1)*ndlu+1),XN(1),ndlu)
       CALL SHIFTD(Xiv((NORDRE-1)*ndlv+1),XN(ndlu+1),ndlv)
       do iordre=1,NORDRE
        CALL SHIFTD(Xiu((iordre-1)*ndlu+1),Xray((iordre-1)*ndl+1),ndlu)
        CALL SHIFTD(Xiv((iordre-1)*ndlv+1),
     &        Xray((iordre-1)*ndl+ndlu),ndlv)
       enddo
c -------          
       
       call NORMER_F(X1,ndl,NrmX1)
       call NORMER_F(XN,ndl,NrmXN)
       call CALCUL_RAY(CRIT,Xray,ndl,NORDRE,RAYON)

c        RAYON=(CRIT*NrmX1/NrmXN)**(1./(NORDRE-1))
c        if(RAYON.gt.0.1)RAYON=0.1
c -------calcul de la solution   


c -------   
       CALL SHIFTD(X0u,X0,ndlu)
       CALL SHIFTD(X0v,X0(ndlu+1),ndlv)
c -------          
       NPOI=1

       CALL SENSCO(X1,X,IP,APOSITIF,RAYON,NPOI,DA,ndl)
       
       IF (pilot.EQ.1)THEN
        CMAX=1.
        CALL VERIF_CHARGE(C0,CMAX,IP,APOSITIF,DA,NPOI,DA)
       ELSEIF(pilot .EQ.2)THEN
       CMAX=1.
       APOS    =APOSITIF 
       CTARGET =CMAX
       ALEFT   =0.D0
       ARIGHT  =RAYON
       
       CALL VERIF_CMAX(C0,Ci,CMAX,APOS,CTARGET,ALEFT,ARIGHT,NORDRE,NPOI,
     &   DA) 
        
       ELSE 
        stop 'limite pilotage non implementé'
       ENDIF 
       A=DA !*DBLE(IPOI)
       
       
       CALL SOLPOL(X0v,Xiv,Xv,ndlv,NORDRE,A)
       CALL SOLPOL(X0u,Xiu,Xu,ndlu,NORDRE,A)
       CALL SOLPOL(C0,Ci,LAMBDA,1,NORDRE,A)
    
       CALL VORIEN(Xiu,X(1),APOSITIF,NORDRE,ndlu,RAYON)
       CALL VORIEN(Xiv,X(1+ndlu),APOSITIF,NORDRE,ndlv,RAYON)
c -------       
       CALL SHIFTD(LAMBDA,C0,1)
       CALL SHIFTD(Xu,X0u,ndlu)
       CALL SHIFTD(Xv,X0v,ndlv)
       CALL SHIFTD(Xu,X(1),ndlu)
       CALL SHIFTD(Xv,X(ndlu+1),ndlv)       
c -------     

       call DEFORME(X0u,X0v,Xiu,Xiv,IP,ndpnu,ndpnv,ndlu,ndlv)
       
       if(mod(IP-1,dincsto).eq.0.or.IP.eq.1.or.IP.eq.nbincr)then 
       write(14,1)C0,Xu(nctr),Xv(ndpnv*(nctr-1)+1),
     &    Xv(ndpnv*(nctr-1)+2),
     &    Xv(ndpnv*(nctr-1)+1)* YOUNG*Iqua/long**3 
       endif

       WRITE(*,54) DA,LAMBDA

      
      IF (pilot.EQ.1.OR.pilot.EQ.2)THEN
      
       IF (DABS(LAMBDA-CMAX).LT.1.D-5) stop 'Arret normal'
       
      ELSEIF(pilot.EQ.3)THEN
      
c        IF (VALAT) 
       
      ENDIF       

       IP=IP+1
      
      enddo !do while(IP.le.nbincr)

      end

      
      
      SUBROUTINE VORIEN(V,U,APOSITIF,NOR,NDL,RAYON)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %   ON FORME UN VECTEUR ORIENTE DANS LE SENS DU PARCOURS          %
COMM %   A L'AIDE DE dV/dA EN FIN DE PAS (A=RAYON)                     %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION V(NDL,*),U(*)
C
       IF(APOSITIF.EQ.1.D0) THEN
        A=RAYON
       ELSE
        A=-RAYON
       ENDIF
       
       DO K=1,NDL
        U(K)=V(K,1)
       ENDDO
       
       DO   J=2,NOR 
       
        PA=1.D0
        DO  JJ=1,J-1
        PA=PA*A  
        ENDDO
c         On obtient A**(J-1)

        PA=PA*J
        DO  K=1,NDL
        U(K)=U(K)+V(K,J)*PA
        ENDDO
c         On obtient dU/dA
        ENDDO 
     
       IF(APOSITIF.EQ.1.D0) THEN
       
       ELSE
        DO  K=1,NDL
         U(K)=-U(K)
        ENDDO
       ENDIF
       
      END 
      
      
      
      SUBROUTINE VERIF_CMAX(C0,CI,CMAX,APOS,CTARGET,ALEFT,ARIGHT,
     & NORDRE,N,DA)

COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %																 %
COMM %    LONGUEUR D'ARC IMPOSÉE                                       %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION CI(*)
	
	DOUBLE PRECISION C0,CL,PRECISION,FAR
     &    ,FAL,A,C,FA,CR,DA,RAY1,CMAX
	INTEGER N,ITER_MAX,ITER,NORDRE,IP 
	REAL*8 APOS,CTARGET,ALEFT,ARIGHT
		
	
	RAY1=ARIGHT
	N=1

	
	
C CALCUL DE LAMBDA_MIN, LAMBDA_MAX
	

      ITER=0
      ITER_MAX=20000
	PRECISION = 1.D-6
	
	ALEFT=0
	CALL SOLPOL(C0,CI,CL,1,NORDRE,ALEFT)
	ARIGHT=DA
	CALL SOLPOL(C0,CI,CR,1,NORDRE,ARIGHT)
	FAR=CMAX-CR
	FAL=CMAX-CL

c 	write(*,*)FAL,FAR
	
	COND1=FAR
	
	 
	IF( COND1 .LT. 0)THEN

	DO WHILE(DABS(FAR-FAL).GT. PRECISION .AND. ITER.LT.ITER_MAX)

		
	A = (ALEFT+ARIGHT)/DBLE(2.D0)

	CALL SOLPOL(C0,CI,C,1,NORDRE,A)
	FA=CMAX-C

	IF(FAL*FA .GT. 0.D0) THEN
	ALEFT=A
	ELSE
	ARIGHT=A
	ENDIF

	CALL SOLPOL(C0,CI,CL,1,NORDRE,ALEFT)
	CALL SOLPOL(C0,CI,CR,1,NORDRE,ARIGHT)
	FAR=CMAX-CR
	FAL=CMAX-CL
	ITER=ITER+1

	ENDDO
	DA=A

	ENDIF

        CALL SOLPOL(C0,CI,CR,1,NORDRE,A)

	
	END
	
	
      SUBROUTINE VERIF_CHARGE(C,CMAX,IP,APOSITIF,RAYON,NPOI,DA)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %   ON DETERMINE LE SENS DE PARCOURT ET LE PAS                    %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      IF ((C+DABS(RAYON)).GT.CMAX) RAYON=CMAX-C  
      
       DA=DABS(RAYON)/DBLE(NPOI)
      END
COMM ..................../////////////............................

      SUBROUTINE SOLPOL(V0,V,U,NDIM,NOR,A)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %   ON FORME LA SOLUTIONS POLYNOMIALES --> U=V0+A*V(1)+...+ApV(p) %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DOUBLE PRECISION V0(*),V(NDIM,*),U(*) 
      
      CALL SHIFTD(V0,U,NDIM)
      DO 2  J=1,NOR 
      CALL DACTUA(U,V(1,J),NDIM,A**J)
 2    CONTINUE 
      END
      
      SUBROUTINE SENSCO(V,U,IP,APOSITIF,RAYON,NPOI,DA,NDL)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %   ON DETERMINE LE SENS DE PARCOURT ET LE PAS                    %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*),U(*)
      IF(IP.LT.1) THEN
       APOSITIF=1.D0
       DA=RAYON/DBLE(NPOI)          
      ELSE      
        PROJ=PROSCA(V,U,NDL)
        IF(PROJ.GE.0) THEN
         APOSITIF=1.D0
         DA=RAYON/DBLE(NPOI)       
        ELSE

        APOSITIF=0.D0
        DA=-RAYON/DBLE(NPOI)
        ENDIF
        
      ENDIF
	
      END
      
