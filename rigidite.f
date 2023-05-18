      subroutine RTG(MTG,X0u,X0v,nel,ndpnu,ndpnv,ndlu,ndlv,ndl)
      implicit none
      include 'maillage.inc' 
c       include 'varstok.inc'
      real*8  xn(nel),yn(nel),VEGU(nel*ndpnu),VEGV(nel*ndpnv),Kuu(2,2),
     & Kvv(4,4),Kuv(2,4),KvG(4,4) 
      real*8  KGu(1000*1000),KGv(1000*1000),KGuv(1000*1000),
     & KGeo(1000*1000),MTG(ndl,*),X0u(*),X0v(*)
      integer i,j,k,l,m,n,nel,ndpn,nbint,ndl,IM,ndpnu,ndpnv,ndlu,ndlv
      
      ndl=ndlu+ndlv
      
      KGu=0.
      KGv=0.
      KGuv=0.
      KGeo=0.
c       MTG=0.
      
1     FORMAT(162(G15.8,1X))    
2     FORMAT(500(G15.8,1X))
    
      do IM=1,nmai
       Kuu=0.
       Kuv=0.
       Kvv=0. 
       KvG=0.
       call RT_el(X0u,X0v,Kuu,Kvv,Kuv,KvG,IM,nel,ndpnu,ndpnv)
       
       call ASSMATNXN(Kuu,KGu ,IM,nel,ndpnu,ndlu)
       call ASSMATNXN(Kvv,KGv ,IM,nel,ndpnv,ndlv)
       call ASSMATNXN(KvG,KGeo,IM,nel,ndpnv,ndlv)
       call ASSMATNXM(Kuv,KGuv,IM,nel,ndpnu,ndpnv,ndlu,ndlv)
       
       call ASSKG(MTG,KGu,KGv,KGuv,ndlu,ndlv,ndl)
      enddo !IM=1,nmai

c       do i=1,n
c       write(*,2)i,B(i)
c       enddo
      

      close (20)
      end
      
      subroutine RT_el(X0u,X0v,Kuu,Kvv,Kuv,KvG,IM,nel,ndpnu,
     & ndpnv)
      implicit none
      include 'donnees.inc'
      include 'maillage.inc' 
c       include 'varstok.inc'      
      real*8 le,xn(nel),yn(nel),GAUSS(3,16),POIDS(16),TETA1,TETA2,TETA3,
     & Jc,N(2),NK(2),MNN(2,2),Kuu(2,*),H(4),HK(4),HKK(4),Kuv(2,*),coef1,
     & MHH(4,4),MNH(2,4),Kvv(4,*),KvG(4,*),UKG,VKG,VEGU(2),
     & VEGV(4),XG,YG,NG0,X0u(*),X0v(*)
      integer i,j,k,l,m,nel,ndpnu,ndpnv,ndpn,nbint,ndle,ndl,IM,IG
      
1     FORMAT(162(G15.8,1X))   

      l=1
      do j=1,nel
       xn(j)=xnoe(MAIL(IM,j))
       yn(j)=ynoe(MAIL(IM,j))         
       do k=1,ndpnu
        VEGU(l)=X0u(ndpnu*(MAIL(IM,j)-1)+k)
        l=l+1
       enddo 
      enddo  
      
      l=1
      do j=1,nel        
       do k=1,ndpnv
        VEGV(l)=X0v(ndpnv*(MAIL(IM,j)-1)+k)
        l=l+1
       enddo 
      enddo       
      
      le=dsqrt((xn(2)-xn(1))**2+(yn(2)-yn(1))**2)      

c -------- calcul de Kuu        
      call GRIGCO_4x1(83,GAUSS,POIDS,nbint)
      do IG=1,nbint
       TETA1=GAUSS(1,IG)
       TETA2=GAUSS(2,IG)
       TETA3=GAUSS(3,IG)

       Jc=le/2.
       coef1= Jc*POIDS(IG)*YOUNG*Aire       
       
       call FORMS1(TETA1,N,NK) 
       CALL MatNH(NK,NK,MNN,2,2)
      
       Do i=1,2
        do j=1,2
         Kuu(i,j)=Kuu(i,j)+MNN(i,j)*coef1/Jc**2
        enddo
       enddo

      enddo !IG=1,nbint
      
c -------- calcul de KEvv      
      call GRIGCO_4x1(83,GAUSS,POIDS,nbint)  
      
      do IG=1,nbint
       TETA1=GAUSS(1,IG)
       TETA2=GAUSS(2,IG)
       TETA3=GAUSS(3,IG)
     
       Jc=le/2.
       coef1= Jc*POIDS(IG)*YOUNG*Iqua       
       call FORMS3(TETA1,H,HK,HKK,le)
       CALL MatNH(HKK,HKK,MHH,4,4)
       Do i=1,4
        do j=1,4
         Kvv(i,j)=Kvv(i,j)+MHH(i,j)*coef1/Jc**4
        enddo
       enddo
       

      enddo !IG=1,nbint      
      
c -------calcul de KGvv
      call GRIGCO_16x1(83,GAUSS,POIDS,nbint)  
      
      do IG=1,nbint
       TETA1=GAUSS(1,IG)
       TETA2=GAUSS(2,IG)
       TETA3=GAUSS(3,IG)
       
       call FORMS1(TETA1,N,NK)
       UKG=0.
       DO k=1,nel       
        UKG = UKG+NK(k)*VEGU(k)       
       ENDDO !k=1,nel 
       
       Jc=le/2.
       UKG=UKG/Jc

       call FORMS3(TETA1,H,HK,HKK,le) 
       
       XG=0.
       DO k=1,nel       
        XG = XG+H(2*k-1)*xn(k) 
       ENDDO !k=1,nel       
       
       VKG=0.
       DO k=1,nel       
        VKG = VKG+HK(2*k-1)*VEGV(2*k-1)+HK(2*k)*VEGV(2*k)  
       ENDDO !k=1,nel  
       VKG=VKG/Jc      
       NG0=YOUNG*Aire *(UKG+.5*VKG**2)
c        write(20,1)XG,NG0
       coef1= Jc*POIDS(IG) 
       
       CALL MatNH(HK,HK,MHH,4,4)
       
       Do i=1,4
        do j=1,4
         Kvv(i,j)=Kvv(i,j)+nlgm*
     &            MHH(i,j)*(NG0+YOUNG*Aire*VKG**2)*coef1/Jc**2
         KvG(i,j)=KvG(i,j)+nlgm*
     &            MHH(i,j)*(NG0+YOUNG*Aire*VKG**2)*coef1/Jc**2     
        enddo
       enddo
      enddo !IG=1,nbint   
      
    
c -------calcul de Kuv     
      call GRIGCO_4x1(83,GAUSS,POIDS,nbint)
      do IG=1,nbint
       TETA1=GAUSS(1,IG)
       TETA2=GAUSS(2,IG)
       TETA3=GAUSS(3,IG)
       
       Jc=le/2.
       call FORMS1(TETA1,N,NK)
       call FORMS3(TETA1,H,HK,HKK,le)
       CALL MatNH(NK,HK,MNH,2,4)
       coef1= Jc*POIDS(IG)*YOUNG*Aire 

      
       VKG=0.
       DO k=1,nel       
        VKG = VKG+HK(2*k-1)*VEGV(2*k-1)+HK(2*k)*VEGV(2*k)  
       ENDDO !k=1,nel  

       VKG=VKG/Jc
       
       
       Do i=1,2
        do j=1,4
         Kuv(i,j)=Kuv(i,j)+MNH(i,j)*nlgm*VKG*coef1/Jc**2
        enddo
       enddo       

      enddo !IG=1,nbint

      
 2     FORMAT(162(G15.8,1X))        
            
      
      end
      
      subroutine datloc(xn,yn,Uu,Uv,VEGU,VEGV,IM,nel,ndpnu,ndpnv)
      implicit none
      include 'maillage.inc'    
      real*8  xn(1),yn(1),Uu(1),Uv(1),VEGU(1),
     & VEGV(1)
      integer i,j,k,l,m,n,IM,nel,ndpn,ndl,ndpnu,ndpnv
      
       i=1
       
       do j=1,nel
        xn(j)=xnoe(MAIL(IM,j))
        yn(j)=ynoe(MAIL(IM,j))
        
       do k=1,ndpnu
        VEGU(i)=Uu(ndpnu*(MAIL(IM,j)-1)+k)
        i=i+1
       enddo       
        
       enddo 
       
       l=1
       do j=1,nel

       do k=1,ndpnv
        VEGV(l)=Uv(ndpnv*(MAIL(IM,j)-1)+k)
        l=l+1
       enddo         
        
       enddo        
      
      end
      
      
      
      SUBROUTINE MatNH(N,H,MNH,nN,nH) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION N(*),H(*),MNH(nN,*)
COMM      

      DO 1 I=1,nN
      DO 1 J=1,nH
1     MNH(I,J)=0.D0

      call TVECVEC(N,H,MNH,nN,nH)  
c       
c  2     FORMAT(162(G15.8,1X)) 
      END       
