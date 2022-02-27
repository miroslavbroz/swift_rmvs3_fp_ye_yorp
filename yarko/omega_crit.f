c omega_crit.f
c Check critical rotation rate.
c Miroslav Broz (miroslav.broz@email.cz), Jan 21st 2010

      subroutine omega_crit(ntp, omega_1, omega_2, iseed, istat)

      include '../swift.inc'
      include 'const.inc'
      include 'spin.inc'
      include 'yarko.inc'
      include 'yorp.inc'

c inputs
      integer ntp,iseed,istat(NTPMAX,NSTAT)
      real*8 omega_1, omega_2

c temporary variables
      integer i
      real*8 omega_g, omega_s, omega_c, P, omega_new
      real*8 s_f, alpha, beta, kappa, k, C

c functions
      real ran1

c Holsappple (2005); Eqs. (5.7), (5.9), (5.10)
      s_f = 0.25d0
      alpha = 0.7d0
      beta = 0.7d0
      kappa = 2.25d7 * 1.d-3*1.d-2**(-1.d0/2.d0)  ! g cm^-1/2 -> kg m^-1/2

      C = (alpha*beta)**(1.d0/3.d0) * sqrt(5.*(3.*s_f*(1.+beta**2)
     :  - sqrt(3.*(1.-beta**2+beta**2)))
     :  / (3.*s_f**2*(1.+beta**2)**2-1.+beta**2-beta**4))

c main
      do i = 1,ntp
        if (istat(i,1).eq.0) then

c gravity regime
          omega_g = sqrt(4.d0/3.d0*PI*capG*rho_bulk(i))

c strength regime
          k = kappa*(R(i)/1.d-2)**(-1.d0/2.d0)
          omega_s = C*sqrt(k/rho_bulk(i))*(1.d0/R(i))
          omega_c = omega_s+omega_g

c the sense of rotation stays the same after a mass shedding...
          if (abs(omega(i)).gt.omega_c) then
            omega_new = omega_1 + (omega_2-omega_1)*ran1(iseed)
            if (omega(i).lt.0.d0) then
              omega_new = -omega_new
            endif
            omega(i) = omega_new

c change a Gaussian sphere to another one (due to a change of shape)
            if (gauss_rnd) then
              fg_id(i) = int(1.0d0 + nGAUSS*ran1(iseed))
              fg_id(i) = min(max(fg_id(i),1),nGAUSS)
            endif

          endif

        endif ! istat
      enddo   ! ntp

      return
      end


