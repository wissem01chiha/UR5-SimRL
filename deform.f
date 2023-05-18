      SUBROUTINE DEFORME(Xu,Xv,Xiu,Xiv,IP,ndpnu,ndpnv,ndlu,ndlv)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      include 'donnees.inc'
      include 'maillage.inc'
c       include 'varstok.inc'
      integer I,iordre,IP
      real*8 time,Xiu(ndlu,*),Xiv(ndlv,*),Xu(*),Xv(*),val
      CHARACTER*30 sweep_blanks
      CHARACTER*20 FICHIER1,FICHIER2,FICHIER3,tm,ord
      
      OPEN(2,FILE='DEF_SIM.plt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      
      Imilieme=IP/1000
	  Icentaine=MOD(IP,1000)/100
	  Idizaine=MOD(IP,100)/10
      Ireste=MOD(IP,10) 

          
       time=IP
       Write( tm, '(G15.4)' )time
       tm=sweep_blanks(tm)
       FICHIER1='Inc-'//trim(tm)//'.s' 


      OPEN(3,FILE=FICHIER1)
      
      WRITE(2,*)'plot ''Inc-'//trim(tm)//'.s'' u 1:2 w l ls -1 lw 2',
     & ' title ''' //trim(tm) //' sec'''
     

      do iordre=1,NORDRE
      Write( ord, '(I5)' )iordre+1
      ord=sweep_blanks(ord)      
      WRITE(2,*)'replot ''Inc-'//trim(tm)//'.s'' u 1:'//trim(ord)//
     & ' w lp lw 2',' title ''' //trim(ord) //' sec'''     
      enddo     
c       WRITE(2,*)"#pause 0.5"  

      fmax=0.
      do I=1,nbnoe
       do iordre=1,NORDRE
        val=max(dabs(Xiv(ndpnv*(I-1)+1,iordre)),val)
c         if(val.gt.fmax)fmax=val
       enddo
      
      enddo
      
c       DO I=1,NMAI
      do I=1,nbnoe
      
      write(3,40) xnoe(I),
     &   Xu(ndpnu*(I-1)+1),Xv(ndpnv*(I-1)+1),Xv(ndpnv*(I-1)+2),
     &   (Xiv(ndpnv*(I-1)+1,iordre)/val,iordre=1,NORDRE)     
c       write(3,40) xnoe(I2),ynoe(I2)      
c       WRITE(3,*)'    '
     
      enddo
      write(3,40)
      write(3,40) 'valmax',val
c       ENDDO
      
      close(3)
      close (2)
      
40    FORMAT(100(1X,G15.8))      
     
      END     
      
      
      character(30) function sweep_blanks(in_str)
      character(*), intent(in) :: in_str
      character(30) :: out_str
      character :: ch
      integer :: j

      out_str = " "
      do j=1, len_trim(in_str)
      ! get j-th char
      ch = in_str(j:j)
      if (ch .ne. " ") then
        out_str = trim(out_str) // ch
      endif
      sweep_blanks = out_str 
      end do
      end function sweep_blanks      
