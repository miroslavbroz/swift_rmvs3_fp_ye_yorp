c*************************************************************************
c STEP_KDK2.F
c*************************************************************************
c Integrator of the 2nd order SBAB2 (see Laskar & Robutel, 2000).
c This subroutine takes a step in helio coordinates; both PL and TP.
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c
c Remarks: Adopted from martin's nbwh.f program
c Authors:  Hal Levison 
c Date:    2/19/93
c Last revision: 9/26/97
c
c Modified: Miroslav Broz, miroslav.broz@email.cz
c Date: Oct 10th 2001

      subroutine step_kdk2(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals
      integer i1sttp,i
      real*8 xint(NPLMAX,INTMAX),yint(NPLMAX,INTMAX),zint(NPLMAX,INTMAX)

c----
c...  Executable code 

      i1sttp = i1st

c...  first do the planets
      call step_kdk2_pl(i1st,nbod,mass,j2rp2,j4rp4,
     &     xint,yint,zint,
     &     xh,yh,zh,vxh,vyh,vzh,dt)

      if(ntp.ne.0) then

c...     next the test particles
         call step_kdk2_tp(i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &        xint,yint,zint,
     &        xht,yht,zht,vxht,vyht,vzht,istat,dt)

      endif

      return
      end   ! step_kdk2
c------------------------------------------------------------------------

