c*************************************************************************
c GETACC_YARKO.F
c*************************************************************************
c This subroutine calculates and ADDS the THERMAL acceleration for the
c test particles in the HELIOCENTRIC frame. Yarkovsky diurnal and seasonal
c effects are included.
c
c See: Vokrouhlicky, D., Astron. Astrophys. 335, 1093-1100 (1998) and
c Vokrouhlicky, D., Farinella, P., Astron. J. 118, 3049-3060 (1999). 
c
c Input:
c   nbod		number of massive bodies
c   ntp			number of test particles
c   mass(nbod)		masses array
c   xh,yh,zh		massive bodies' heliocentric coordinates
c   xht,yht,zht		test particles' heliocentric coordinates
c   istat		status of the test paricles
c   irht		inverse distances of tp's from the Sun
c			(were already calculated earlier in getacch.f)
c Output:
c   axht,ayht,azht	new tp's accelerations in helio coord
c
c Created by:  Miroslav Broz, miroslav.broz@email.cz
c Date:  Aug 6th 2002

      subroutine getacc_yarko(nbod,ntp,mass,xh,yh,zh,
     &     xht,yht,zht,istat,irht,axht,ayht,azht)

      include '../swift.inc'
      include '../filter/filter.inc'
      include 'yarko.inc'
      include 'spin.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX)
      real*8 mass(NPLMAX),xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)

c  input orbital elements last computed in io_write_filter.f subroutine
c  in common block /elements/, we need also the epoch of mean anomaly t0
      include '../filter/cb_oscel.inc'

c  following common block /times/ serves another parameter for seasonal
c  yarkovsky effect computation
      real*8 t
      common/times/ t
c      real*8 dt
c      common/dttimes/ dt

c...  Outputs:
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)
                
c...  Internals:
      integer i,k
      real*8 irht(NTPMAX)
      real*8 sin_theta_0,cos_theta_0,Tstar4,lambda,rho,n0(3),f1(3),
     &  ex(3),ey(3),x1,Tstar3,ath,n,l,kl,gms,s0(3)
      real*8 a(3)
      complex*8 psi2
      real*8 pm

cc  OMP functions
c!$    logical omp_in_parallel
c!$    real *8 omp_get_wtime
cc  OMP timings, debugging
c!$    real *8 tbegin, tend
c!$    common /timing/ tbegin, tend

c----
c...  Executable code 

      if (.not.use_yarko) return

c!$    tend = omp_get_wtime()
c!$    write(*,*) 'getacc_yarko begin: ', tend - tbegin, ' sec'

      gms = sqrt(mass(1))

      if (is_forward) then
        pm = 1.d0
      else
        pm = -1.d0
      endif

!$OMP PARALLEL DO
!$OMP&  PRIVATE (i, n0, cos_theta_0, sin_theta_0, Tstar4, Tstar3,
!$OMP&    lambda, psi2, rho, f1, x1, ey, ath, n, l, s0)
!$OMP&  SHARED (ntp, istat, xht, yht, zht, irht, s, axht, ayht, azht,
!$OMP&    mass, osc_elmts, t, t0, koefsin, koefcos)
      do i = 1, ntp
        if (istat(i).eq.0) then		! only active particles
c
c  Yarkovsky diurnal force
c
          n0(1) = -xht(i)*irht(i)
          n0(2) = -yht(i)*irht(i)
          n0(3) = -zht(i)*irht(i)
          s0(1) = s(1,i)
          s0(2) = s(2,i)
          s0(3) = s(3,i)
          if (omega(i).lt.0.d0) then
            s0(1) = -s0(1)
            s0(2) = -s0(2)
            s0(3) = -s0(3)
          endif
          cos_theta_0 = n0(1)*s0(1)+n0(2)*s0(2)+n0(3)*s0(3)
          sin_theta_0 = sqrt(abs(1.d0-cos_theta_0**2))            ! NOTE the abs() due to round-off!
          Tstar4 = Tstar40(i)*irht(i)**2
          Tstar3 = sqrt(sqrt(Tstar4*Tstar4*Tstar4))               ! Tstar3=Tstar4**(3.d0/4.d0)
          lambda = lambda0(i)/Tstar3
          psi2 = 1.d0/(1.d0+(lambda/(1.d0+lambda))*psi(i))
          rho = rho0(i)*Tstar4/(1.d0+lambda)

          f1(1) =  rho*sin_theta_0*dble(psi2)
          f1(2) = -rho*sin_theta_0*imag(psi2)
          f1(3) =  rho*cos_theta_0
c
c  transformation to ecliptic coordinate frame
c
          if (sin_theta_0.ne.0.d0) then
            x1 = 1.d0/sin_theta_0
            do k = 1,3
              ex(k) = (n0(k)-cos_theta_0*s(k,i))*x1
            enddo
            ey(1) = (s0(2)*n0(3)-s0(3)*n0(2))*x1
            ey(2) = (s0(3)*n0(1)-s0(1)*n0(3))*x1
            ey(3) = (s0(1)*n0(2)-s0(2)*n0(1))*x1

            a(1) = pm*(ex(1)*f1(1)+ey(1)*f1(2)+s0(1)*f1(3))
            a(2) = pm*(ex(2)*f1(1)+ey(2)*f1(2)+s0(2)*f1(3))
            a(3) = pm*(ex(3)*f1(1)+ey(3)*f1(2)+s0(3)*f1(3))
          else
            a(1) = 0.d0
            a(2) = 0.d0
            a(3) = 0.d0
          endif
c
c  Yarkovsky seasonal force
c
          ath = 0.0d0
          n = gms * 1.0d0/sqrt(osc_elmts(1,i)**3)                  ! n =  sqrt(mass(1)) * osc_elmts(1,i)**(-1.5d0)
          l = osc_elmts(6,i) + n * (t-t0)

c  koefcos, koefsin, which depend on slow orbital elements are calculated
c  every dtfilter in yarko_seasonal.f

          do k = 1,k_seasonal
            kl = k*l
            ath = ath + koefcos(k,i)*cos(kl) + koefsin(k,i)*sin(kl)
          enddo

c  seasonal thermal force is aligned with spin axis

          ath = pm*ath
          a(1) = a(1) + ath * s0(1)
          a(2) = a(2) + ath * s0(2)
          a(3) = a(3) + ath * s0(3)

          if (.not.isnan(a(1))) then
            axht(i) = axht(i) + a(1)
            ayht(i) = ayht(i) + a(2)
            azht(i) = azht(i) + a(3)
          else
            write(*,*) 'Warning: NaN in getacc_yarko!',
     :        ' a close encounter with Vesta?'
          endif

        endif ! istat
      enddo   ! ntp
!$OMP END PARALLEL DO

c!$    tend = omp_get_wtime()
c!$    write(*,*) 'getacc_yarko end: ', tend - tbegin, ' sec'

      return
      end      ! getacc_yarko

c---------------------------------------------------------------------


