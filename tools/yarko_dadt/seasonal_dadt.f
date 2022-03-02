c**********************************************************************
c SEASONAL_DADT.F
c**********************************************************************
c Calculate seasonal yarkovsky drift of semimajor axis.
c
c Modified (from yarko_seasonal.f): Miroslav Broz, miroslav.broz@email.cz
c Date: Oct 10th 2001

      subroutine seasonal_dadt(ntp,mass,istat,elmts,dadt)

      include 'swift.inc'
      include 'yarko.inc'
      include 'spin.inc'

c...  Inputs: 
      integer ntp
      real*8 mass(NPLMAX)	! we need only mass(1) from this array
      integer istat(NTPMAX)
      real*8 elmts(6,1:NTPMAX)
      real*8 dadt(1:NTPMAX)

c...  Internals:
      integer i,k
      real*8 alpha(MAXk_seasonal), beta(MAXk_seasonal)
      real*8 e(MAXk_seasonal-1)		! powers of eccentricity
      real*8 sQ,sP,Xk,AXk,BXk,CXk,DXk,CXkDXk2,sinXk,cosXk,X1,
     :  expXk,sindeltak,cosdeltak,ath0,sininc,cosinc,sinperi,cosperi,
     :  sinnode,cosnode,n,lambdas,rdash,eta,dasum
c      real*8 Gk,AXkBXk,CXkDXk

c----
c...  Executable code 

      do i = 1, ntp
        if (istat(i).eq.0) then		! only active particles

c  precalculate powers of eccentricity
          e(1) = elmts(2,i)
          if (e(1).gt.1.d0-TINY) goto 999	! exclude parabolic and hyperbolic orbits
          do k = 2,6
            e(k) = e(k-1)*e(1)
          enddo

c  calculate alpha and beta coefficients
        alpha(1) = 1.d0
     :    -3.0d0 * (0.5d0**3)*e(2)
     :    +5.0d0/6.0d0 * (0.5d0**2)*e(4)
     :    -7.0d0/72.0d0 * (0.5d0**7)*e(6)
        alpha(2) = 4.0d0 * (
     :    +2.0d0 * (0.5d0**2)*e(1)
     :    -16.0d0/3.0d0 * (0.5d0**4)*e(3)
     :    +4.0d0 * (0.5d0**6)*e(5))
        alpha(3) = 9.0d0 * (
     :    +3.0d0 * (0.5d0**3)*e(2)
     :    -45.0d0/4.0d0 * (0.5d0**5)*e(4)
     :    +567.0d0/40.0d0 * (0.5d0**7)*e(6))
        alpha(4) = 16.0d0 * (
     :    +16.0d0/3.0d0 * (0.5d0**4)*e(3)
     :    -128.0d0/5.0d0 * (0.5d0**6)*e(5))
        alpha(5) = 25.0d0 * (
     :    +125.0d0/12.0d0 * (0.5d0**5)*e(4)
     :    -4375.0d0/72.0d0 * (0.5d0**7)*e(6))
        alpha(6) = 36.0d0 * (108.0d0/5.0d0 * (0.5d0**6)*e(5))
        alpha(7) = 49.0d0 * (16807.0d0/360.0d0 * (0.5d0**7)*e(6))

c  beta(k) / eta, see sQeta
        beta(1) = 1.d0 - e(2)/8.0d0 + e(4)/192.0d0 - e(6)/9216.0d0
        beta(2) = 4.0d0/2.0d0*e(1) * (1.d0 - e(2)/3.0d0 + e(4)/24.0d0)
        beta(3) = 27.0d0/8.0d0*e(2) * (1.d0 - 9.0d0/16.0d0*e(2)
     :    + 81.0d0/640.0d0*e(4))
        beta(4) = 16.0d0/3.0d0*e(3) * (1.d0 - 4.0d0/5.0d0*e(2))
        beta(5) = 25.0d0*125.0d0/384.0d0*e(4)
     :    * (1.d0 - 25.0d0/24.0d0*e(2))
        beta(6) = 972.0d0/80.0d0*e(5)
        beta(7) = 823543.0d0/46080.0d0*e(6)

c  calculate P, Q vectors from inc, node, peri
        eta = sqrt(1.d0-e(1)**2)
        sininc = sin(elmts(3,i))
        cosinc = cos(elmts(3,i))
        sinperi = sin(elmts(4,i))
        cosperi = cos(elmts(4,i))
        sinnode = sin(elmts(5,i))
        cosnode = cos(elmts(5,i))
        sP = (cosperi*cosnode-sinperi*sinnode*cosinc) * s(1,i)
     :    + (cosperi*sinnode+sinperi*cosnode*cosinc) * s(2,i)
     :    + (sinperi*sininc) * s(3,i)
        sQ = (-sinperi*cosnode-cosperi*sinnode*cosinc) * s(1,i)
     :    + (-sinperi*sinnode+cosperi*cosnode*cosinc) * s(2,i)
     :    + (cosperi*sininc) * s(3,i)

        n = sqrt(mass(1)) * elmts(1,i)**(-1.5d0)
        Rdash = Rdash0(i) * sqrt(n/86400.0d0)
        X1 = sqrt(2.0d0*1.0d0) * Rdash
        lambdas = thetas(i)*(elmts(1,i)**(3.d0/4.d0))*
     :    eta**(3.d0/4.d0)/X1
        ath0 = 4.d0/9.d0*rhos0(i)/(elmts(1,i)**2)/(1.d0+lambdas)
     :    * 86400.d0**2/149597870.d3

        dasum = 0.0d0

        do k = 1,k_seasonal
          Xk = sqrt(2.0d0*k) * Rdash
          sinXk = sin(Xk)
          cosXk = cos(Xk)
c  canceled factor exp(Xk)
          expXk = exp(-Xk)
          AXk = -expXk*(Xk+2.d0)-((Xk-2.d0)*cosXk-Xk*sinXk)
          BXk = -expXk*Xk-(Xk*cosXk+(Xk-2.d0)*sinXk)
          CXk = AXk + lambdas/(1.d0+lambdas)*(expXk*3.d0*(Xk+2.d0)+
     :      (3.d0*(Xk-2.d0)*cosXk+Xk*(Xk-3.d0)*sinXk))
          DXk = BXk + lambdas/(1.d0+lambdas)*(expXk*Xk*(Xk+3.d0)-
     :      (Xk*(Xk-3.d0)*cosXk-3.d0*(Xk-2.d0)*sinXk))
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
c          cosdeltak = sqrt(1.d0-sindeltak**2.d0)

c          koefcos(k,i) = ath0 * (sP*alpha(k)*cosdeltak +
c     :      sQ*eta*beta(k)*sindeltak)
c          koefsin(k,i) = ath0 * (-sP*alpha(k)*sindeltak +
c     :      sQ*eta*beta(k)*cosdeltak)

          dasum = dasum + sindeltak / k *
     :      (sP**2*alpha(k)**2+sQ**2*(eta*beta(k))**2)
        enddo

        dadt(i) = 4.d0/(9.d0*n)*rhos0(i)/(elmts(1,i)**2)/(1.d0+lambdas)
     :    * dasum *86400.d0**2/149597870.d3 * 365.25d6

999     continue

        endif
      enddo

      return
      end

c---------------------------------------------------------------------

