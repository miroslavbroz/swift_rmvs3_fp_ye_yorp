c*************************************************************************
c DIURNAL_DADT.F
c*************************************************************************
c This subroutine calculates diurnal drift of semimajor axes
c
c Modified (from getacc_yarko.f):  Miroslav Broz, miroslav.broz@email.cz
c Date:  Oct 10th 2001

      subroutine diurnal_dadt(nbod,ntp,mass,xh,yh,zh,
     &     xht,yht,zht,istat,irht,elmts,gamma,dadt)

      include 'swift.inc'
      include 'yarko.inc'
      include 'spin.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX)
      real*8 mass(NPLMAX),xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 elmts(6,NTPMAX)
      real*8 gamma(NTPMAX)

c...  Outputs:
      real*8 dadt(NTPMAX)
                
c...  Internals:
      integer i,k
      real*8 irht(NTPMAX)
      real*8 sin_theta_0,cos_theta_0,Tstar4,lambda,rho,n0(3),f1(3),
     &  ex(3),ey(3),x1,Tstar3,ath,n,l,s0(3)
      double complex psi2

c----
c...  Executable code 

      do i=1,ntp
        if (istat(i).eq.0) then		! only active particles
c
c  yarkovsky diurnal force
c
          n0(1) = -xht(i)*irht(i)
          n0(2) = -yht(i)*irht(i)
          n0(3) = -zht(i)*irht(i)
          s0(1) = sign(s(1,i), omega(i))
          s0(2) = sign(s(2,i), omega(i))
          s0(3) = sign(s(3,i), omega(i))
          cos_theta_0 = n0(1)*s0(1)+n0(2)*s0(2)+n0(3)*s0(3)
          sin_theta_0 = sqrt(abs(1.d0-cos_theta_0**2))            ! NOTE the abs() due to round-off!
          Tstar4 = Tstar40(i)*irht(i)**2
          Tstar3 = sqrt(sqrt(Tstar4*Tstar4*Tstar4))               ! Tstar3=Tstar4**(3.d0/4.d0)
          lambda = lambda0(i)/Tstar3
          psi2 = 1.d0/(1.d0+(lambda/(1.d0+lambda))*psi(i))
          rho = rho0(i)*Tstar4/(1.d0+lambda)

          n = sqrt(mass(1)) * elmts(1,i)**(-1.5d0)
          dadt(i) = 365.25d6 * rho*imag(psi2) * 2.0d0/n * cos(gamma(i))

        endif
      enddo

      return
      end

c---------------------------------------------------------------------

