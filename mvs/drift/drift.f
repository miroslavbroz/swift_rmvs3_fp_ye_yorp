c*************************************************************************
c                        DRIFT.F
c*************************************************************************
c This subroutine loops thorugh the particles and calls the danby routine
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xj,yj,zj      ==>  initial position in jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial position in jacobi coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xj,yj,zj      ==>  final position in jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final position in jacobi coord 
c                                       (real arrays)
c
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 9/5/94
c
c Modified: Miroslav Broz, miroslav.broz@email.cz
c Remarks: Added OpenMP parallelization directives and slightly
c          modified code for etajm1, etaj variables.
c Date: Oct 24rd 2001

      subroutine drift(nbod,mass,xj,yj,zj,vxj,vyj,vzj,dt)

      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod),dt

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Internals:
c      real*8 etajm1,etaj
      real*8 etaj(NTPMAX),mu
      integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth

c      etajm1 = mass(1)
      etaj(1) = mass(1)
      do j = 2, nbod
        etaj(j) = etaj(j-1) + mass(j)
      enddo

!$OMP PARALLEL DO
!$OMP& PRIVATE (j, mu, iflg)
!$OMP$ SHARED (nbod, mass, xj, yj, zj, vxj, vyj, vzj, dt, etaj)
      do j = 2, nbod
c        etaj = etajm1 + mass(j)
c        mu = mass(1)*etaj/etajm1
        mu = mass(1)*etaj(j)/etaj(j-1)
        call drift_one(mu,xj(j),yj(j),zj(j),vxj(j),vyj(j),vzj(j),
     &    dt,iflg)
        if (iflg.ne.0) then
          write(*,*) ' Planet ',j,' is lost !!!!!!!!!'
          write(*,*) mu,dt
          write(*,*) xj(j),yj(j),zj(j)
          write(*,*) vxj(j),vyj(j),vzj(j)
          write(*,*) ' STOPPING '
          call util_exit(1)
        endif
c        etajm1 = etaj
      enddo
!$OMP END PARALLEL DO

      return
      end
c--------------------------------------------------------------------------

