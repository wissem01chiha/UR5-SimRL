      program spinfem
      implicit none
      include 'donnees.inc'
      fdat='test.inp'
      fdat=trim(fdat)
      call lcpar
      call mesh()
      call increm()
      


      end
