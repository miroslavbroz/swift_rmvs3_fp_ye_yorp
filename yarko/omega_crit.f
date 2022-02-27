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
      real*8 omega_c, P_1, P_2, P, omega_new

c functions
      real*8 ran1

c main
      P_1 = TWOPI/omega_1
      P_2 = TWOPI/omega_2

      do i = 1,ntp
        if (istat(i,1).eq.0) then

          omega_c = sqrt(8.d0/3.d0*PI*capG*rho_bulk(i))

c the sense of rotation stays the same after a mass shedding...
          if (abs(omega(i)).gt.omega_c) then
            P = P_1 + (P_2-P_1)*ran1(iseed)
            omega_new = TWOPI/P
            if (omega(i).lt.0.d0) then
              omega_new = -omega_new
            endif
            omega(i) = omega_new

c            omega(i) = sign(omega_1 + (omega_2-omega_1)*ran1(iseed),
c     :        omega(i))

c change a Gaussian sphere to another one (due to a change of shape)
            if (gauss_rnd) then
              fg_id(i) = int(0.999d0 + nGAUSS*ran1(iseed))
            endif

          endif

c 2DO: model a transition to a strength regime too (below 200 m)!

        endif ! istat
      enddo   ! ntp

      return
      end


