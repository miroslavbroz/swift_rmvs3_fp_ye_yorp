c**********************************************************************
c YARKO_SEASONAL.F
c**********************************************************************
c Calculate seasonal yarkovsky coefficients, which do NOT depend
c on mean anomaly.
c
c See Vokrouhlicky, D., Farinella, P., Astron. J. 118, 3049-3060 (1999). 
c
c Input:
c  ntp			number of TP's
c  mass(NPLMAX)		solar mass(1) is necessary only
c
c  thermal parameters in common block /yarko/
c  osculating orbital elements from io_write_filter.f in /elements/
c
c Output:
c  koefcos, koefsin in /yarko/
c
c Modified: Miroslav Broz, miroslav.broz@email.cz
c Date: Aug 6th 2002

      subroutine yarko_seasonal(ntp,mass,istat)

      include '../swift.inc'
      include '../filter/filter.inc'
      include 'yarko.inc'
      include 'spin.inc'

c...  Inputs: 
      integer ntp
      real*8 mass(NPLMAX)	! we need only mass(1) from this array
      integer istat(NTPMAX,NSTAT)

c  input orbital elements last computed in io_write_filter.f subroutine
c  in common block /elements/, we need also the epoch of mean anomaly t0
      include '../filter/cb_oscel.inc'

c...  Internals:
      integer i,k
      real*8 alpha(MAXk_seasonal), beta(MAXk_seasonal)
      real*8 e(MAXk_seasonal-1)		! powers of eccentricity
      real*8 sQ,sP,Xk,AXk,BXk,CXk,DXk,CXkDXk2,sinXk,cosXk,X1,
     :  expXk,sindeltak,cosdeltak,ath0,sininc,cosinc,sinperi,cosperi,
     :  sinnode,cosnode,n,lambdas,rdash,eta
c      real*8 Gk,AXkBXk,CXkDXk

c  debug code
c      real*8 radius(NTPMAX),da
c      common /debug/ radius

c----
c...  Executable code 

!$OMP PARALLEL DO
!$OMP& DEFAULT (PRIVATE)
!$OMP& SHARED (ntp, osc_elmts, s, mass, Rdash0, thetas, rhos0,
!$OMP&   koefcos, koefsin, k_seasonal)
      do i = 1, ntp
        if (istat(i,1).eq.0) then

c  precalculate powers of eccentricity
        e(1) = osc_elmts(2,i)
        if (e(1).gt.1.d0-TINY) goto 999	! exclude parabolic and hyperbolic orbits
        do k = 2,6
          e(k) = e(k-1)*e(1)
        enddo

c  calculate alpha and beta coefficients
        alpha(1) = 1
     :    -3.0 * (0.5**3)*e(2)
     :    +5.0/6.0 * (0.5**2)*e(4)
     :    -7.0/72.0 * (0.5**7)*e(6)
        alpha(2) = 4.0 * (
     :    +2.0 * (0.5**2)*e(1)
     :    -16.0/3.0 * (0.5**4)*e(3)
     :    +4.0 * (0.5**6)*e(5))
        alpha(3) = 9.0 * (
     :    +3.0 * (0.5**3)*e(2)
     :    -45.0/4.0 * (0.5**5)*e(4)
     :    +567.0/40.0 * (0.5**7)*e(6))
        alpha(4) = 16.0 * (
     :    +16.0/3.0 * (0.5**4)*e(3)
     :    -128.0/5.0 * (0.5**6)*e(5))
        alpha(5) = 25.0 * (
     :    +125.0/12.0 * (0.5**5)*e(4)
     :    -4375.0/72.0 * (0.5**7)*e(6))
        alpha(6) = 36.0 * (108.0/5.0 * (0.5**6)*e(5))
        alpha(7) = 49.0 * (16807.0/360.0 * (0.5**7)*e(6))

c  beta(k) / eta, see sQeta
        beta(1) = 1 - e(2)/8.0 + e(4)/192.0 - e(6)/9216.0
        beta(2) = 4.0/2.0*e(1) * (1 - e(2)/3.0 + e(4)/24.0)
        beta(3) = 27.0/8.0*e(2) * (1 - 9.0/16.0*e(2) + 81.0/640.0*e(4))
        beta(4) = 16.0/3.0*e(3) * (1 - 4.0/5.0*e(2))
        beta(5) = 25.0*125.0/384.0*e(4) * (1 - 25.0/24.0*e(2))
        beta(6) = 972.0/80.0*e(5)
        beta(7) = 823543.0/46080.0*e(6)

c  calculate P, Q vectors from inc, node, peri
        eta = sqrt(1-e(1)**2)
        sininc = sin(osc_elmts(3,i))
        cosinc = cos(osc_elmts(3,i))
        sinperi = sin(osc_elmts(4,i))
        cosperi = cos(osc_elmts(4,i))
        sinnode = sin(osc_elmts(5,i))
        cosnode = cos(osc_elmts(5,i))
        sP = (cosperi*cosnode-sinperi*sinnode*cosinc) * s(1,i)
     :    + (cosperi*sinnode+sinperi*cosnode*cosinc) * s(2,i)
     :    + (sinperi*sininc) * s(3,i)
        sQ = (-sinperi*cosnode-cosperi*sinnode*cosinc) * s(1,i)
     :    + (-sinperi*sinnode+cosperi*cosnode*cosinc) * s(2,i)
     :    + (cosperi*sininc) * s(3,i)

        n = sqrt(mass(1)) * osc_elmts(1,i)**(-1.5d0)
        Rdash = Rdash0(i) * sqrt(n/86400.0d0)
        X1 = sqrt(2.0d0) * Rdash
        lambdas = thetas(i)*(osc_elmts(1,i)**(3./4.))*
     :    eta**(3./4.)/X1
        ath0 = 4./9.*rhos0(i)/(osc_elmts(1,i)**2)/(1+lambdas)
     :    * 86400.**2/149597870e3

c        da = 0.0d0

        do k = 1,k_seasonal
          Xk = sqrt(2.0d0*k) * Rdash
          sinXk = sin(Xk)
          cosXk = cos(Xk)
c  canceled factor exp(Xk)
          expXk = exp(-Xk)
          AXk = -expXk*(Xk+2)-((Xk-2)*cosXk-Xk*sinXk)
          BXk = -expXk*Xk-(Xk*cosXk+(Xk-2)*sinXk)
          CXk = AXk + lambdas/(1+lambdas)*(expXk*3*(Xk+2)+
     :      (3*(Xk-2)*cosXk+Xk*(Xk-3)*sinXk))
          DXk = BXk + lambdas/(1+lambdas)*(expXk*Xk*(Xk+3)-
     :      (Xk*(Xk-3)*cosXk-3*(Xk-2)*sinXk))
c canceled factor AXkBXk
c          AXkBXk = sqrt(AXk**2+BXk**2)
c          CXkDXk = sqrt(CXk**2+DXk**2)
c          Gk = AXkBXk/CXkDXk
c          cosdeltak = (AXk*CXk+BXk*DXk)/(AXkBXk*CXkDXk)
c          sindeltak = (BXk*CXk-AXk*DXk)/(AXkBXk*CXkDXk)
          CXkDXk2 = CXk**2+DXk**2
          cosdeltak = (AXk*CXk+BXk*DXk)/CXkDXk2
          sindeltak = (BXk*CXk-AXk*DXk)/CXkDXk2
c          sindeltak = 0
c          cosdeltak = sqrt(1.d0-sindeltak**2)

          koefcos(k,i) = ath0 * (sP*alpha(k)*cosdeltak +
     :      sQ*eta*beta(k)*sindeltak)
          koefsin(k,i) = ath0 * (-sP*alpha(k)*sindeltak +
     :      sQ*eta*beta(k)*cosdeltak)

c          da = da + sindeltak / k *
c     :      (sP**2*alpha(k)**2+sQ**2*(eta*beta(k))**2)
        enddo

c        da = 4./(9.*n)*rhos0(i)/(osc_elmts(1,i)**2)/(1+lambdas)
c     :    * da *86400.**2/149597870e3
c        write(*,*) radius(i),abs(da*365.25d6)

999     continue

        endif  ! istat
      enddo    ! ntp
!$OMP END PARALLEL DO

c      stop

      return
      end      ! yarko_seasonal.f

c---------------------------------------------------------------------


