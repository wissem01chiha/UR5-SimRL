      subroutine lcpar
      implicit none
      include 'donnees.inc'
      integer i,j,k,l,m,n,deb,fin,pas,nlist,serie(100),nserie
      real*8 val
      character*8 phra,mot1,mot2,mot3
      
      open(unit=1,file=fdat)
      read(1,*)
      read(1,*)  long,lrg,epai
                 Iqua=(lrg*epai**3)/12.
                 Aire=lrg*epai
      read(1,*)
      read(1,*)  nmaix
      read(1,*) !-------Materiau 
      read(1,*)  YOUNG, NU
                 Gcis=YOUNG/2/(1+NU)
      read(1,*) nlgm
      read(1,*) !-------LIAISON
      read(1,*)  nlist
      LIAI=0
      LIAIV=0
      k=1
      do i=1,nlist
       val=0.
       read(1,*)deb,fin,pas
       read(1,*)phra,val

       do j=deb,fin,pas
        LIAI(1,k)=j
        if(phra(1:2).eq.'U')then
         LIAI(2,k)=1
         LIAIV(k)=val       
        elseif(phra(1:2).eq.'V')then
         LIAI(2,k)=2
         LIAIV(k)=val    
        elseif(phra(1:2).eq.'RV')then
         LIAI(2,k)=3
         LIAIV(k)=val           
        endif              
        k=k+1     
       enddo
      enddo
      NLIAI=k-1

      
      read(1,*) !-------CHARGEMENT
      read(1,*) nlist
      k=1
      do i=1,nlist
       read(1,*)mot1
       TYPLOAD(k)=mot1
       if(TYPLOAD(k).ne.'FORCE'.and.TYPLOAD(k).ne.'DISP')then
        stop 'Type CHARGEMENT non implemente'
       endif
       read(1,*)nserie,(serie(j),j=1,nserie),val,phra
       do j=1,nserie
        LOAD(1,k)=serie(j)
        if(phra(1:1).eq.'1')then
         LOAD(2,k)=1
        elseif(phra(1:1).eq.'2')then
         LOAD(2,k)=2
        elseif(phra(1:1).eq.'3')then
         LOAD(2,k)=3
        elseif(phra(1:1).eq.'4')then
         LOAD(2,k)=4         
        endif  
        LOADV(k)=val
        k=k+1
       enddo
      enddo
      NLOAD=k-1
 

      read(1,*) !-------ANM
      read(1,*) CRIT,NORDRE,pilot,nbincr
      read(1,*) !-------VISUALISATION
      read(1,*) nctr,xctr,dincsto

      
      close(1)
      
      
      end
      
      
      
      subroutine mesh()
      implicit none
      include 'donnees.inc'
      include 'maillage.inc'
      real*8 tlx,tly,PROSCA
      integer noe,i,j,Icont,k
      open(unit=1,file='mail.ini' )
      tlx=long/nmaix
      noe=0
      do i=1,nmaix+1
        noe=noe+1
        xnoe(noe)=(i-1)*tlx
        ynoe(noe)=0.

c         write(*,1)noe,xnoe(noe),ynoe(noe)
      enddo
      
      
      nbnoe=noe
!       -----------Assemblage--------------
      nmai=nmaix
      Icont=0
      do i=1,nmaix
       Icont=Icont+1
       MAIL(Icont,1)=Icont
       MAIL(Icont,2)=Icont+1
      enddo
      do i=1,Icont
      write(1,3)xnoe(MAIL(i,1)),ynoe(MAIL(i,1)),MAIL(i,1)   
     
      write(1,3)xnoe(MAIL(i,2)),ynoe(MAIL(i,2)),MAIL(i,2) 
      write(1,3) 
      enddo
      
      close(1)
      
1     FORMAT(I4,4(1X,G15.8))      
2     FORMAT(I4,4(1X,G4.8)) 
3     FORMAT(1X,4(G15.8),6(G20.8),I5)

      end  
      
      
      
      
      

