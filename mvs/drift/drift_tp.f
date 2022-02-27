c*************************************************************************
c                        DRIFT_TP.F
c*************************************************************************
c This subroutine loops thorugh the TEST particles and calls the danby routine
c
c             INPUT:
c                 ntp             ==>  number of test particles (int scalar)
c                 msun            ==>  mass of the sun (real scalar)
c                 xjt,yjt,zjt     ==>  initial position in jacobi coord 
c                                      (real arrays)
c                 vxjt,vyjt,vzjt  ==>  initial position in jacobi coord 
c                                      (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt              ==>  time step
c             OUTPUT:
c                 xjt,yjt,zjt     ==>  final position in jacobi coord 
c                                       (real arrays)
c                 vxjt,vyjt,vzjt  ==>  final position in jacobi coord 
c                                       (real arrays)
c
c Authors:  Hal Levison 
c Date:    2/18/93
c Last revision:
c
c Modified: Miroslav Broz, miroslav.broz@email.cz
c Remarks: Added OpenMP parallelization directives.
c Date: Oct 24rd 2001

      subroutine drift_tp(ntp,msun,xjt,yjt,zjt,vxjt,vyjt,vzjt,dt,istat)

      include '../../swift.inc'

c...  Inputs Only: 
      integer ntp
      real*8 msun,dt

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)

c...  Internals:
      integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth

!$OMP PARALLEL DO
!$OMP& PRIVATE (j, iflg)
!$OMP& SHARED (ntp, istat, msun, xjt, yjt, zjt, vxjt, vyjt, vzjt, dt)
      do j = 1, ntp
        if (istat(j,1).eq.0) then
          call drift_one(msun,xjt(j),yjt(j),zjt(j),
     &      vxjt(j),vyjt(j),vzjt(j),dt,iflg)
          if (iflg.ne.0) then
            istat(j,1) = 1
            istat(j,2) = -1
          endif
        endif
      enddo
!$OMP END PARALLEL DO

      return
      end
c--------------------------------------------------------------------------

