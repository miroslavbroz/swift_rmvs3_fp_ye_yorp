c*************************************************************************
c STEP_KDK2_PL.F
c*************************************************************************
c Integrator of the 2nd order SBAB2 (see Laskar & Robutel, 2000).
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 xint,yint,zint ==> positions of PL's for step_kdk2_tp
c
c Remarks: Adopted from martin's nbwhnew.f program
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 2/24/94
c
c Modified: Miroslav Broz, miroslav.broz@email.cz
c Date: Oct 10th 2001

      subroutine step_kdk2_pl(i1st,nbod,mass,j2rp2,j4rp4,
     &  xint,yint,zint,
     &  xh,yh,zh,vxh,vyh,vzh,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

      real*8 xint(NPLMAX,INTMAX),yint(NPLMAX,INTMAX),zint(NPLMAX,INTMAX)

c...  Internals:
      real*8 axh(NPLMAX),ayh(NPLMAX),azh(NPLMAX)
      real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
      real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)
      integer i, j

      save axh,ayh,azh,xj,yj,zj     ! Note this !!

c----
c...  Executable code 

      if (i1st.eq.0) then
c...    Convert to jacobi coords
        call coord_h2j(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &    xj,yj,zj,vxj,vyj,vzj)
c...    Get the accelerations in helio frame. if frist time step
        call getacch(nbod,mass,j2rp2,j4rp4,xj,yj,zj,
     &    xh,yh,zh,axh,ayh,azh)
        i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a d1*dt 
      call kickvh(nbod,vxh,vyh,vzh,axh,ayh,azh,d_koef(1)*dt)

      do i = 1, 2

c...    Remember these coordinates for step_kdk2_tp subroutine
        do j = 1, nbod
          xint(j,i) = xh(j)
          yint(j,i) = yh(j)
          zint(j,i) = zh(j)
        enddo

c...    Convert the helio. vels. to Jac. vels. (the Jac. positions are unchanged)
        call coord_vh2vj(nbod,mass,vxh,vyh,vzh,vxj,vyj,vzj)

c..     Drift in Jacobi coords for the step c1*dt
        call drift(nbod,mass,xj,yj,zj,vxj,vyj,vzj,c_koef(i)*dt)

c...    After drift, compute helio. xh and vh for acceleration calculations
        call coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &    xh,yh,zh,vxh,vyh,vzh)

c...    Get the accelerations in helio frame.
        call getacch(nbod,mass,j2rp2,j4rp4,xj,yj,zj,xh,yh,zh,
     &    axh,ayh,azh)

c...    Apply a heliocentric kick for a d2*dt
        call kickvh(nbod,vxh,vyh,vzh,axh,ayh,azh,d_koef(i+1)*dt)

      enddo

c...  Remember PL's coordinates after final drift
      do j = 1, nbod
        xint(j,3) = xh(j)
        yint(j,3) = yh(j)
        zint(j,3) = zh(j)
      enddo

      return
      end   ! step_kdk2_pl
c---------------------------------------------------------------------

