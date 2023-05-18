      subroutine FNLG(Fnl,Xiu,Xiv,X0u,X0v,nel,ndpnu,ndpnv,ORDRE,ndlu,
     &  ndlv)
      implicit none
      include 'maillage.inc' 
c       include 'varstok.inc'      
      real*8  xn(nel),yn(nel),VEGU(nel*ndpnu),VEGV(nel*ndpnv),
     &  Fnlu(1000),Fnlv(1000),Fnl(*),Fenlu(2),Fenlv(4),X0u(*),X0v(*),
     & Xiu(ndlu,*),Xiv(ndlv,*)
      integer i,j,k,l,m,n,nel,ndpn,nbint,ndl,IM,ndpnu,ndpnv,ndlu,ndlv,
     & ORDRE 
1     FORMAT(162(G15.8,1X))    

      Fnlu=0.
      Fnlv=0.
      
      ndl=ndlu+ndlv
      call DZERO(Fnl,ndl) 
      do IM=1,nmai
      Fenlu=0.
      Fenlv=0.
      call FNL_el(Fenlu,Fenlv,Xiu,Xiv,X0u,X0v,IM,nel,ndpnu,ndpnv,
     & ndlu,ndlv,ndl,ORDRE)
      call ASSVEC(Fnlu,Fenlu,IM,nel,ndpnu,ndlu)
      call ASSVEC(Fnlv,Fenlv,IM,nel,ndpnv,ndlv)
      do i=1,ndl
       if(i.le.ndlu)then
        Fnl(i)=Fnlu(i)
       else
        Fnl(i)=Fnlv(i-ndlu)
       endif
      enddo
      
      enddo !IM=1,nmai
      
c        do i=1,ndlu
c         write(*,1)i,Fnlu(i)
c        enddo   
c        pause 'fffhhh'        
      end
      subroutine FNL_el(Fenlu,Fenlv,Xiu,Xiv,X0u,X0v,IM,nel,ndpnu,ndpnv,
     &  ndlu,ndlv,ndl,ORDRE)
      implicit none
      include 'donnees.inc'   
      include 'maillage.inc'
      
      real*8 le,xn(nel),yn(nel),GAUSS(3,16),POIDS(16),TETA1,TETA2,TETA3,
     & Jc,N(2),NK(2),H(4),HK(4),HKK(4),coef1,Fenlu(*),Fenlv(*),VKG0,
     & VKG1,VKG2,VEG0(nel*ndpnv),Xiu(ndlu,*),Xiv(ndlv,*),X0u(*),
     & X0v(*),VEG1(nel*ndpnv),VEG2(nel*ndpnv),UEG(nel*ndpnu),UKG,
     & XG,YG,Ni,terme
      integer i,j,k,l,ll,m,nel,ndpnu,ndpnv,ndpn,nbint,ndle,IM,IG,ii,
     & iordre,ORDRE,sord,ndlu,ndlv,ndl      
      
1     FORMAT(162(G15.8,1X))       
c -------- calcul de FNLu   
      call GRIGCO_8x1(83,GAUSS,POIDS,nbint)
      
      do IG=1,nbint
       TETA1=GAUSS(1,IG)
       TETA2=GAUSS(2,IG)
       TETA3=GAUSS(3,IG)

       call FORMS1(TETA1,N,NK) 
       
        do iordre=1,ORDRE-1
c  --------
       
        l=1
        m=1
        ll=1 
        do j=1,nel
         xn(j)=xnoe(MAIL(IM,j))
         yn(j)=ynoe(MAIL(IM,j)) 
         
         do k=1,ndpnv
          VEG0(l)=X0v(ndpnv*(MAIL(IM,j)-1)+k)
          VEG1(l)=Xiv(ndpnv*(MAIL(IM,j)-1)+k,iordre)
          VEG2(l)=Xiv(ndpnv*(MAIL(IM,j)-1)+k,ORDRE-iordre)
          l=l+1
         enddo 
        
        enddo 
c  --------    

         le=dsqrt((xn(2)-xn(1))**2+(yn(2)-yn(1))**2)
         Jc=le/2.
         
         call FORMS3(TETA1,H,HK,HKK,le)
         VKG1=0.           
         DO k=1,nel       
          VKG1 = VKG1+HK(2*k-1)*VEG1(2*k-1)+HK(2*k)*VEG1(2*k)               
         ENDDO !k=1,nel  
         VKG1=VKG1/Jc   
      
         VKG2=0. 
         DO k=1,nel       
          VKG2 = VKG2+HK(2*k-1)*VEG2(2*k-1)+HK(2*k)*VEG2(2*k)             
         ENDDO !k=1,nel  
         VKG2=VKG2/Jc 
         
c          write(*,1) IM,HK(1),HK(2),HK(3),HK(4)
         coef1= Jc*POIDS(IG)*YOUNG*Aire 
         
         do i=1,nel*ndpnu
          Fenlu(i)=Fenlu(i)-nlgm*.5*VKG2*VKG1*NK(i)*coef1/Jc
         enddo
         
        enddo !iordre=1,NORDRE-1
      enddo !IG=1,nbint    


c -------- calcul de FNLv      
      call GRIGCO_16x1(83,GAUSS,POIDS,nbint)
      
      
      do IG=1,nbint
       TETA1=GAUSS(1,IG)
       TETA2=GAUSS(2,IG)
       TETA3=GAUSS(3,IG)
        do iordre=1,ORDRE-1
c  --------      
        l=1
        do j=1,nel
         xn(j)=xnoe(MAIL(IM,j))
         yn(j)=ynoe(MAIL(IM,j))         
         do k=1,ndpnv
          VEG0(l)=X0v(ndpnv*(MAIL(IM,j)-1)+k)
          VEG1(l)=Xiv(ndpnv*(MAIL(IM,j)-1)+k,iordre)
          VEG2(l)=Xiv(ndpnv*(MAIL(IM,j)-1)+k,ORDRE-iordre)
          l=l+1
         enddo 
        enddo 
     
c  -------- 
         le=dsqrt((xn(2)-xn(1))**2+(yn(2)-yn(1))**2)
         Jc=le/2.
         call FORMS3(TETA1,H,HK,HKK,le) 
         call FORMS1(TETA1,N,NK)         
         
         VKG0=0. 
         DO k=1,nel       
          VKG0 = VKG0+HK(2*k-1)*VEG0(2*k-1)+HK(2*k)*VEG0(2*k)           
         ENDDO !k=1,nel  
         VKG0=VKG0/Jc  
         
         VKG1=0. 
         DO k=1,nel       
          VKG1 = VKG1+HK(2*k-1)*VEG1(2*k-1)+HK(2*k)*VEG1(2*k)           
         ENDDO !k=1,nel  
         VKG1=VKG1/Jc   
      
         VKG2=0. 
         DO k=1,nel       
          VKG2 = VKG2+HK(2*k-1)*VEG2(2*k-1)+HK(2*k)*VEG2(2*k)           
         ENDDO !k=1,nel  
         VKG2=VKG2/Jc
         
         terme=YOUNG*Aire*.5*VKG0*VKG1
  
c -----------calcul de l'effort normale Ni
         
c  --------      
         l=1
         do j=1,nel
           do k=1,ndpnu
           UEG(l)=Xiu(ndpnu*(MAIL(IM,j)-1)+k,iordre)
           l=l+1
          enddo 
         enddo
         call FORMS1(TETA1,N,NK) 
         UKG=0.
         DO k=1,nel       
          UKG = UKG+NK(k)*UEG(k)       
         ENDDO !k=1,nel 
         UKG=UKG/Jc
         
         Ni=YOUNG*Aire*UKG
         Ni=Ni+YOUNG*Aire*VKG0*VKG1
         
         VKG0=VKG2
c  --------          
         do sord=1,iordre
          l=1
          do j=1,nel
           do k=1,ndpnv
            VEG1(l)=Xiv(ndpnv*(MAIL(IM,j)-1)+k,sord)
            VEG2(l)=Xiv(ndpnv*(MAIL(IM,j)-1)+k,iordre-sord)
            l=l+1
           enddo 
          enddo           
c  --------  
          VKG1=0. 
          DO k=1,nel       
           VKG1 = VKG1+HK(2*k-1)*VEG1(2*k-1)+HK(2*k)*VEG1(2*k)           
          ENDDO !k=1,nel  
          VKG1=VKG1/Jc  
          
          VKG2=0. 
          DO k=1,nel       
           VKG2 = VKG2+HK(2*k-1)*VEG2(2*k-1)+HK(2*k)*VEG2(2*k)           
          ENDDO !k=1,nel  
          VKG2=VKG2/Jc      
          Ni=Ni+.5*YOUNG*Aire*VKG1*VKG2
        enddo !sord=1,iordre 

        coef1= Jc*POIDS(IG)  
        do i=1,ndpnv*nel
         Fenlv(i)=Fenlv(i)-(Ni+terme)*VKG0*HK(i)/Jc*coef1*nlgm
        enddo         
       
       enddo !iordre=1,ORDRE-1
       
      enddo !IG=1,nbint      
      
      
      
      end
