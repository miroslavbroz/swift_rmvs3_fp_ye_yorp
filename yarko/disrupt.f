c**********************************************************************
c DISRUPT.F
c**********************************************************************
c Decide on disruption of TP's and set corresponing flags.
c
c Input:
c  t			time
c  ntp			number of test particles
c  dt			disruption time step
c  tau(NTPMAX)		mean disruption times
c  idum			number to initialize ran1()
c
c Output:
c  istat(NTPMAX,NSTAT)	new stat-array for each TP with disruption flag
c  rstat(NTPMAX,NSTAT)	not yet used
c  idum			number to initialize ran1()
c  disrupted		.true. if any particle was disrupted, otherwise .false.
c
c Remarks: 
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Nov 8th 2000

      subroutine disrupt(t,ntp,dt,tau,idum,istat,rstat,disrupted)

      include '../swift.inc'

c  input
      integer ntp,idum,istat(NTPMAX,NSTAT)
      real*8 t,dt,tau(NTPMAX),rstat(NTPMAX,NSTAT)
      logical disrupted

c  internal
      integer i
      real*8 p,p_disr
      real ran1

c  ignore disruptions if tau < 0
      if (tau(1).lt.0.d0) return

c  main
      disrupted = .false.

      do i=1,ntp
        if (istat(i,1).eq.0) then
          p_disr = 1 - exp(-dt/tau(i))
          p  =1.0d0*ran1(idum)

c  discard the i-th TP due to disruption
          if (p.lt.p_disr) then
            disrupted = .true.
            istat(i,1) = 1
            istat(i,2) = -5
          endif

        endif ! istat
      enddo   ! ntp

      return
      end	! disrupt.f
c-----------------------------------------------------------------


