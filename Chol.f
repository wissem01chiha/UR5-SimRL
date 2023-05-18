!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
      MODULE LU

      CONTAINS

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
       Subroutine LUDCMP(A,N,INDX,D,CODE)
       PARAMETER(NMAX=1d5,TINY=1.5D-30)
       REAL*8  AMAX,DUM, SUM, A(N,*),VV(NMAX)
       INTEGER CODE, D, INDX(*)

       D=1; CODE=0

       DO I=1,N
         AMAX=0.d0
         DO J=1,N
           IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
         END DO ! j loop
         IF(AMAX.LT.TINY) THEN
           CODE = 1
           RETURN
         END IF
         VV(I) = 1.d0 / AMAX
       END DO ! i loop

       DO J=1,N
         DO I=1,J-1
           SUM = A(I,J)
           DO K=1,I-1
             SUM = SUM - A(I,K)*A(K,J) 
           END DO ! k loop
           A(I,J) = SUM
         END DO ! i loop
         AMAX = 0.d0
         DO I=J,N
           SUM = A(I,J)
           DO K=1,J-1
             SUM = SUM - A(I,K)*A(K,J) 
           END DO ! k loop
           A(I,J) = SUM
           DUM = VV(I)*DABS(SUM)
           IF(DUM.GE.AMAX) THEN
             IMAX = I
             AMAX = DUM
           END IF
         END DO ! i loop  
   
         IF(J.NE.IMAX) THEN
           DO K=1,N
             DUM = A(IMAX,K)
             A(IMAX,K) = A(J,K)
             A(J,K) = DUM
           END DO ! k loop
           D = -D
           VV(IMAX) = VV(J)
         END IF

         INDX(J) = IMAX
         IF(DABS(A(J,J)) < TINY) A(J,J) = TINY
      
         IF(J.NE.N) THEN
           DUM = 1.d0 / A(J,J)
           DO I=J+1,N
             A(I,J) = A(I,J)*DUM
           END DO ! i loop
         END IF 
       END DO ! j loop

       RETURN
       END subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
       Subroutine LUBKSB(A,N,INDX,B)
       REAL*8  SUM, A(N,*),B(*)
       INTEGER INDX(*)

       II = 0

       DO I=1,N
         LL = INDX(I)
         SUM = B(LL)
         B(LL) = B(I)
         IF(II.NE.0) THEN
           DO J=II,I-1
             SUM = SUM - A(I,J)*B(J)
           END DO ! j loop
         ELSE IF(SUM.NE.0.d0) THEN
           II = I
         END IF
         B(I) = SUM
       END DO ! i loop

       DO I=N,1,-1
         SUM = B(I)
         IF(I < N) THEN
           DO J=I+1,N
             SUM = SUM - A(I,J)*B(J)
           END DO ! j loop
         END IF
         B(I) = SUM / A(I,I)
       END DO ! i loop
      
       RETURN
       END subroutine LUBKSB

      END MODULE LU

! end of file lu.f90

!*******************************************************
!* Solving a linear system AX = B by LU decomposition  *
!* with dynamic allocations                            *
!*                                                     *
!*                   F90 version by J-P Moreau, Paris  *
!*                           (www.jpmoreau.fr)         *
!* --------------------------------------------------- *
!* SAMPLE RUN:                                         *
!*                                                     *
!* Input file (test_lu.dat):                           *
!*                                                     *
!*  4                                                  *
!*  8  2    3  12     25.0                             *
!*  2  4    7   0.25  13.25                            *
!*  3  7    3   5     18.0                             *
!* 12  0.25 5   2     19.25                            *
!*                                                     *
!* Output file (test_lu.lst):                          *
!*                                                     *
!* --------------------------------------------------- *
!*  LINEAR SYSTEM TO BE SOLVED:                        *
!* --------------------------------------------------- *
!*  N=4                                                *
!*                                                     *
!*  8.000000  2.000000  3.000000  12.00000  25.00000   *
!*  2.000000  4.000000  7.000000  0.250000  13.25000   *
!*  3.000000  7.000000  3.000000  5.000000  18.00000   *
!*  12.00000  0.250000  5.000000  2.000000  19.25000   *
!*                                                     *
!*  System solution:                                   *
!*                                                     *
!*  X1=  1.000000                                      *
!*  X2=  1.000000                                      *
!*  X3=  1.000000                                      *
!*  X4=  1.000000                                      *
!* --------------------------------------------------- *
!*                                                     *
!* Uses: modules BASIS, LU.                            *
!*******************************************************
      subroutine RESOL_LU(A,B,n,deco)

      USE BASIS
      USE LU
      real*8 A(n,*)   !real matrix (n x n)
      real*8 B(*)     !real vector (n)
      
c       real*8  temp(n+1)  !real temporary vector (n+1)
      integer  INDX(n)  !integer vector (n)
      
c       real*8, pointer ::  A(:,:)   !real matrix (n x n)
c       real*8, pointer ::  B(:)     !real vector (n)
c       real*8, pointer ::  temp(:)  !real temporary vector (n+1)
c       integer,pointer ::  INDX(:)  !integer vector (n)

      integer d, rc, deco
      character*12 input, output
      character*8 s
  

2     FORMAT(500(G15.8,1X))

  !dynamic allocations
c       allocate(A(n,n),stat=ialloc)
c       allocate(B(n),stat=ialloc)
c       allocate(temp(n+1),stat=ialloc)
c       allocate(INDX(n),stat=ialloc)



c       do i=1,n
c       write(19,2)(A(i,j),j=1,n)
c       enddo
c       
c c       do i=1,n
c c       write(*,2)i,B(i)
c c       enddo
c       
c       stop 'nnnnnnnnnnnn'

!call LU decomposition routine
       
       call LUDCMP(A,n,INDX,D,rc)

!call appropriate solver if previous return code is ok
      if (rc.eq.0) then
        call LUBKSB(A,n,INDX,B)
      endif
        if (rc.eq.1) then
          stop ' The matrix is singular, no solution !'
        end if      


! section fin

50    format(' Input data file name (without .dat): ')
60    format(/' N = ',I2/)
200   format('    X',I1,' = ',F9.6)
201   format('    X',I2,'= ',F9.6)

      END

! end of file test_lu.f90






