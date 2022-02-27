c**********************************************************************
c REORIENT.F
c**********************************************************************
c Collisional reorientation of the spin axes of test particles.
c
c Input:
c  t			time
c  ntp			number of test particles
c  s(3,NTPMAX)		array of spin axes orientations from /yarko/
c  dtreorient		reorientation time step
c  tau(NTPMAX)		mean reorientation times
c  beta_1		beta_1 exponent due to tau(omega) dependence
c  iseed			initial number for ran1
c
c Output:
c  s(3,NTPMAX)	new orientations of the spin axes in /yarko/
c
c Remarks: 
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Apr 18th 2001

      subroutine reorient(t,ntp,dt,tau,beta_1,omega_1,omega_2,
     :  iseed,istat,reoriented)

      include '../swift.inc'
      include 'spin.inc'

c  input
      integer ntp,iseed,istat(NTPMAX,NSTAT)
      real*8 t,dt,tau(NTPMAX),beta_1,omega_1,omega_2
      logical reoriented

c  internal
      integer i
      real*8 s1,coss2,sins2,p,p_reor,tau_reor
      real ran1

c  main
      reoriented = .false.

c  ignore reorientations if tau < 0
      if (tau(1).lt.0.d0) return

      do i=1,ntp
        if (istat(i,1).eq.0) then
          tau_reor = tau(i) * abs(omega(i))**beta_1
          p_reor = 1.d0 - exp(-dt/tau_reor)
          p = 1.0d0*ran1(iseed)

c          write(*,*) 'tau_reor: ', t/365.25, dt/365.25,
c     :      tau_reor/365.25, p_reor, p
         
          if (p.lt.p_reor) then     ! new random orientation of the spin axis
            reoriented = .true.
            s1 = TWOPI*ran1(iseed)
            sins2 = 2.0d0*ran1(iseed)-1.0d0
            coss2 = sqrt(1-sins2**2)
            s(1,i) = cos(s1)*coss2
            s(2,i) = sin(s1)*coss2
            s(3,i) = sins2

            omega(i) = omega_1 + (omega_2-omega_1)*ran1(iseed)    ! change rotational periods too...

            write(*,*) 'reorient: Particle id = ', i, ' was reoriented!'

          endif

        endif ! istat
      enddo   ! ntp

      return
      end       ! reorient.f
c-----------------------------------------------------------------


