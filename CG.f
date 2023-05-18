       subroutine RESOL_GC (MTG,SMG,VG ,ndl)
         
       
       !
        implicit none
       
!                integer ( kind = 4 ), parameter :: n = 79
       
        real ( kind = 8 ) a(ndl,ndl)
        real ( kind = 8 ) b(ndl)
        real ( kind = 8 ) bnrm2
        real ( kind = 8 ) err
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ii
        integer ( kind = 4 ) it
        integer ( kind = 4 ) it_max
        integer ( kind = 4 ) job
        integer ( kind = 4 ) nx
        integer ( kind = 4 ) ny
        real ( kind = 8 ) p(ndl)
        real ( kind = 8 ) q(ndl)
        real ( kind = 8 ) r(ndl)
        real ( kind = 8 ) rnrm2
        integer ( kind = 4 ) seed
        real ( kind = 8 ) tol
        real ( kind = 8 ) x(ndl)
        real ( kind = 8 ) x_exact(ndl)
        real ( kind = 8 ) z(ndl)
        integer ndl,n,j
        real*8 SMG(1000),VG(1000),MTG(1000,1000)
         n=ndl
!         write ( *, '(a)' ) ' '
!         write ( *, '(a)' ) 'TEST02'
!         write ( *, '(a)' ) '  Use CG_RC to solve a linear system'
!         write ( *, '(a)' ) '  involving the Wathen matrix.'
       
        nx = 5
        ny = 4
       
!         write ( *, '(a)' ) ' '
!         write ( *, '(a,i6)' ) '  NX = ', nx
!         write ( *, '(a,i6)' ) '  NY = ', ny
!         write ( *, '(a,i6)' ) '  N  = ', n
        
       
!                call wathen ( nx, ny, n, a )
         do i=1,ndl
          do j=1,ndl
           a(i,j)=MTG(i,j)
          enddo
         enddo
!       DO  i=1,ndl
!         write(*,2) i ,(MTG(i,j),j=1,ndl)
!       ENDDO         
        seed = 123456789
!                call    r8vec_uniform_01 ( n, seed, x_exact )
       
!                b(1:n) = matmul ( a(1:n,1:n), x_exact(1:n) )
        do i=1,ndl
         b(i)=SMG(i)
        enddo
       
        x(1:n) = 0.0D+00
!             Parameters for the stopping test.
       !
         it = 0
         it_max = 5000
         tol = 1.0D-08
         bnrm2 = sqrt ( sum ( b(1:n)**2 ) )
       !
!                Set parameters for the CG_RC code.
       !
         job = 1
       !
!                Repeatedly call the CG_RC code, and on return, do what JOB tells you.
       !
         do
       
           call cg_rc ( n, b, x, r, z, p, q, job )
       !
!                Compute q = A * p.
       !
           if ( job == 1 ) then
       
             q(1:n) = matmul ( a(1:n,1:n), p(1:n) )
       !
!                Solve M * z = r.
       !
           else if ( job == 2 ) then
       
             do i = 1, n
               z(i) = r(i) / a(i,i)
             end do
       !
!                Compute r = r - A * x.
       !
           else if ( job == 3 ) then
       
             r(1:n) = r(1:n) - matmul ( a(1:n,1:n), x(1:n) )
       !
!                Stopping test.
       !
           else if ( job == 4 ) then
       
             rnrm2 = sqrt ( sum ( r(1:n)**2 ) )
       
             if ( bnrm2 == 0.0D+00 ) then
               if ( rnrm2 <= tol ) then
                 exit
               end if
             else
               if ( rnrm2 <= tol * bnrm2 ) then
                 exit
               end if
             end if
       
             it = it + 1
       
             if ( it_max <= it ) then
               write ( *, '(a)' ) ' '
               write ( *, '(a)' ) '  Iteration limit exceeded.'
               write ( *, '(a)' ) '  Terminating early.'
               exit
             end if
       
           end if
       
           job = 2
       
         end do
         DO  i=1,n
          VG(i)=x(i)
         ENDDO          
         x_exact=0
c            write ( *, '(a)' ) ' '
c            write ( *, '(a,i5)' ) '  Number of iterations was ', it
c            write ( *, '(a,g14.6)' ) '  Estimated error is ', rnrm2
           
c            err = maxval ( abs ( x_exact(1:n) - x(1:n) ) )
c            write ( *, '(a,g14.6)' ) '  Loo error is ', err
       
!          write ( *, '(a)' ) ' '
!          write ( *, '(a)' ) '     I      X(I)         X_EXACT(I)    
!      &   B(I)'
!          write ( *, '(a)' ) ' '
!          do i = 1, n
!            write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
!      &        i, x(i), x_exact(i), b(i)
!          end do
       
         return
2     FORMAT(500(G15.8,1X))         
       end
       
       
       
       
       subroutine cg_rc ( n, b, x, r, z, p, q, job )
       
       
       !
        implicit none
       
        integer ( kind = 4 ) n
       
        real ( kind = 8 ) alpha
        real ( kind = 8 ) b(n)
        real ( kind = 8 ) beta
        integer ( kind = 4 ) iter
        integer ( kind = 4 ) job
        real ( kind = 8 ) p(n)
        real ( kind = 8 ) pdotq
        real ( kind = 8 ) q(n)
        real ( kind = 8 ) r(n)
        real ( kind = 8 ) rho
        real ( kind = 8 ) rho_old
        integer ( kind = 4 ) rlbl
        real ( kind = 8 ) x(n)
        real ( kind = 8 ) z(n)
             !
!                           Some local variables must be preserved between calls.
       !
        save iter
        save rho
        save rho_old
        save rlbl
       !
!                           Initialization.
!                    Ask the user to compute the initial residual.
       !
        if ( job == 1 ) then
       
          r(1:n) = b(1:n)
       
          job = 3
          rlbl = 2
       !
!                    Begin first conjugate gradient loop.
!                    Ask the user for a preconditioner solve.
       !
        else if ( rlbl == 2 ) then
       
          iter = 1
       
          job = 2
          rlbl = 3
       !
!                    Compute the direction.
!                    Ask the user to compute ALPHA.
!                    Save A*P to Q.
       !
        else if ( rlbl == 3 ) then
       
          rho = dot_product ( r, z )
       
          if ( 1 < iter ) then
            beta = rho / rho_old
            z(1:n) = z(1:n) + beta * p(1:n)
          end if
       
          p(1:n) = z(1:n)
       
          job = 1
          rlbl = 4
       !
!                    Compute current solution vector.
!                    Ask the user to check the stopping criterion.
       !
        else if ( rlbl == 4 ) then
       
          pdotq = dot_product ( p, q )
          alpha = rho / pdotq
          x(1:n) = x(1:n) + alpha * p(1:n)
          r(1:n) = r(1:n) - alpha * q(1:n)
       
          job = 4
          rlbl = 5
       !
!                    Begin the next step.
!                    Ask for a preconditioner solve.
       !
        else if ( rlbl == 5 ) then
       
          rho_old = rho
          iter = iter + 1
       
          job = 2
          rlbl = 3
       
        end if
       
        return
       end
       
             
      subroutine r8vec_uniform_01 ( n, seed, r )
      
      !
        implicit none
      
        integer ( kind = 4 ) n
      
        integer ( kind = 4 ) i
        integer ( kind = 4 ) k
        integer ( kind = 4 ) seed
        real ( kind = 8 ) r(n)
      
        if ( seed == 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
          write ( *, '(a)' ) '  Input value of SEED = 0.'
          stop
        end if
      
        do i = 1, n
      
          k = seed / 127773
      
          seed = 16807 * ( seed - k * 127773 ) - k * 2836
      
          if ( seed < 0 ) then
            seed = seed + 2147483647
          end if
      
          r(i) = real ( seed, kind = 8 ) * 4.656612875D-10
      
        end do
      
        return
      end       
