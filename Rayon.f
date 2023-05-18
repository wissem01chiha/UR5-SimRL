      
      SUBROUTINE CALCUL_RAY(CRIT,V,NDL,NORDRE,RAYON) 
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %   ON DETERMINE LE RAYON DE CONVERGENCE DE LA SERIE              %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION V(NDL,*)
      INTEGER NDL,NORDRE
	
      XNOR1=DSQRT(PROSCA(V(1,1),V(1,1),NDL))
	
      XNORP=DSQRT(PROSCA(V(1,NORDRE),V(1,NORDRE),NDL))

	
      RAYON=(CRIT*(XNOR1/XNORP))**(1.D0/DBLE(NORDRE-1))    

      END
      SUBROUTINE CALCUL_RAY_2(CRIT,V,NDL,NORDRE,RAYON,SRay) 
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %   ON DETERMINE LE RAYON DE CONVERGENCE DE LA SERIE POUR ORDRE 2 %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION V(NDL,*)
      DOUBLE PRECISION Vi(NORDRE),V1,V2,V3,VN,chimal,yamin,wassat,SRay
      integer NDL,NORDRE
      chimal=1.d-15
      yamin=10

	
      do i=1,NORDRE
      Vi(i)=DSQRT(PROSCA(V(1,i),V(1,i),NDL))
      enddo

       Fyamin=fonc(yamin,NORDRE,CRIT,Vi)
       Fchimal=fonc(chimal,NORDRE,CRIT,Vi)	
	
      ITER=0
      ITER_MAX=20000
      precision = 1.d-4**NORDRE

      DO WHILE(dabs(Fchimal-Fyamin).GT. precision.AND. iter.LT.ITER_MAX)

      wassat=(chimal+yamin)/DBLE(2.d0)
 !     Fwassat=fonc(wassat,NORDRE,CRIT,V1,V2,V3,VN)
      Fwassat=fonc(wassat,NORDRE,CRIT,Vi)
      if(Fwassat*Fchimal .GT. 0) then
      chimal=wassat
      else
      yamin=wassat
      endif
      ITER=ITER+1
      if(ITER .gt. ITER_MAX) stop 'Pb Calcul Rayon'	
      ENDDO

      RAYON=wassat

      END

      DOUBLE PRECISION FUNCTION fonc(RAYON,NORDRE,CRIT,Vi)
	
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      integer NORDRE,i,p
      DOUBLE PRECISION Vi(NORDRE),V1,V2,V3,RAYON,CRIT,som

      p=(NORDRE)
!      fonc = (VN*RAYON**(NORDRE)-CRIT*V2*RAYON**2-CRIT*V3*RAYON**3-
!     *CRIT*V1)
      fonc=Vi(NORDRE)*RAYON**p
      do i=1,NORDRE-1
      fonc=fonc-CRIT*Vi(NORDRE-i)*RAYON**(-i+p)
      enddo
	
	
      END
      

      SUBROUTINE CALCUL_RAY_N(CRIT,V,NDL,NORDRE,RAYON,SRay) 
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                 %
COMM %   ON DETERMINE LE RAYON DE CONVERGENCE DE LA SERIE POUR ORDRE N %
COMM %                                                                 %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION V(NDL*NORDRE)
      DOUBLE PRECISION chimal,yamin,wassat,precis,Fyamin,Fchimal,foncN,
     *SRay,x,chim,yam,Fyam,Fchim
      integer NDL,NORDRE,i,j,ITER_MAX,npas,k
      real*8, allocatable :: VR(:,:)
      allocate(VR(NDL,NORDRE))
      chimal=1.d-10
      yamin=SRay
      Fchimal=0
      Fyamin=0
      Fwassat=0
 !     write(*,*)'calcul Rayon'


      DO j=1,NORDRE
      DO i = 1,NDL
      VR(i,j)=V(i+NDL*(j-1))
      ENDDO
      ENDDO

1     continue	

       npas=1000	
       k=0	   
       do i=1,npas-1
       x=x-SRay/real(npas) !(1.d0-SRay)/real(npas)
       chim=x
       Fchim=foncN(chim,NORDRE,CRIT,VR,NDL)

       write(12,53)chim,Fchim	   
       enddo
       close(12)

	   
      Fchimal=foncN(chimal,NORDRE,CRIT,VR,NDL) 

!      Fchimal=fonc(chimal,NORDRE,CRIT,V1,V2,VN)
      Fyamin=foncN(yamin,NORDRE,CRIT,VR,NDL) 

 !     Fyamin=fonc(yamin,NORDRE,CRIT,V1,V2,VN)
	
      ITER=0
      ITER_MAX=500
      precis = 1.d-3

      DO WHILE(dabs(Fchimal-Fyamin).GT. precis.AND. iter.LT.ITER_MAX)
     
      wassat=(chimal+yamin)/DBLE(2.d0)
      Fwassat=foncN(wassat,NORDRE,CRIT,VR,NDL)  
  !    Fwassat=fonc(wassat,NORDRE,CRIT,V1,V2,VN)
      if(Fwassat*Fchimal .GT. 0) then
      chimal=wassat
      else
      yamin=wassat
      endif
      
      ITER=ITER+1

      If (ITER .gt. ITER_MAX) stop 'PB calcul_Ray'  	
      ENDDO

      RAYON=wassat
 !     write(*,*)'Fin calcul Rayon'
 53    FORMAT(8(G15.8,4X))
      END
	  
      DOUBLE PRECISION FUNCTION foncN(RAYON,NORDRE,CRIT,VR,NDL)	
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DOUBLE PRECISION RAYON,CRIT,som(NDL),VR(NDL,NORDRE),Term1,Term2
      integer NORDRE,NDL,i,j,p
      DIMENSION V(NDL)		
      som=0.d0
      foncN	 =0.d0 
      p= NORDRE
      do i=1,NORDRE
      do j=1,NDL
      som(j)=som(j)+VR(j,i)*(RAYON**(i-NORDRE))
      enddo
      enddo

      Term2=DSQRT(PROSCA(som,som,NDL))
      do i=1,NDL
      V(i)=VR(i,NORDRE)
      enddo	  
      Term1=DSQRT(PROSCA(V,V,NDL))

      foncN = Term1-CRIT*Term2
!       write(*,*)NORDRE,NDL
!       stop
	
      END
