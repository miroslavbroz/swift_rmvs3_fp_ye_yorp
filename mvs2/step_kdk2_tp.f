c*************************************************************************
c STEP_KDK2_TP.F
c*************************************************************************
c Integrator of the 2nd order SBAB2 (see Laskar & Robutel, 2000).
c ONLY DOES TEST PARTICLES
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of test bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4    ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xint,yint,zint ==>  positions of PL's for step_kdk2_tp
c                                       (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c             Output:
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Adopted from martin's nbwhnew.f program
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 2/24/94
c
c Modified: Miroslav Broz, miroslav.broz@email.cz
c Remarks: adopted from step_kdk_tp, but 2nd order integrator
c          and common block /times/ for getacc_yarko subroutine
c Date: Oct 11th 2001

      subroutine step_kdk2_tp(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &  xint,yint,zint,
     &  xht,yht,zht,vxht,vyht,vzht,istat,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4  

      real*8 xint(NPLMAX,INTMAX),yint(NPLMAX,INTMAX),zint(NPLMAX,INTMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
      integer i
      real*8 dth
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)

      real*8 t
      common/times/ t

      save axht,ayht,azht     ! Note this !!

c----
c...  Executable code 

      if (i1st.eq.0) then 
c...    Get the accelerations in helio frame.
        call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,
     &    xint(1,1),yint(1,1),zint(1,1),
     &    xht,yht,zht,istat,axht,ayht,azht)
        i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a d1*dt
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,
     &  d_koef(1)*dt) 

      do i = 2, 3

c...    Take a drift forward by c1*dt
        dth = c_koef(i-1)*dt

        call drift_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,
     &    dth,istat)

c...    Get the accelerations in helio frame.
        t = t + dth

        call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,
     &    xint(1,i),yint(1,i),zint(1,i),
     &    xht,yht,zht,istat,axht,ayht,azht)

c...    Apply a heliocentric kick for a d2*dt 
        call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,
     &    d_koef(i)*dt) 

      enddo

      return
      end   ! step_kdk2_tp
c---------------------------------------------------------------------

