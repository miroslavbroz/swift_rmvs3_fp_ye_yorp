c yorp_evolve.f
c Evolve obliquities and spin rates due to YORP effect.
c Miroslav Broz (miroslav.broz@email.cz), Jan 19th 2010

      subroutine yorp_evolve(t,ntp,dt,mass,xht,yht,zht,vxht,vyht,vzht,
     :  istat)

      include '../swift.inc'
      include 'const.inc'
      include 'yorp.inc'
      include 'spin.inc'
      include 'yarko.inc'

c inputs (plus variable in /yarko/ common block!)
      integer ntp,iseed,istat(NTPMAX,NSTAT)
      real*8 t,dt
      real*8 mass(NPLMAX)
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

c locals
      integer i, ialpha, IER
      real*8 x(3), v(3), n(3), n_p(3), s_new(3)
      real*8 obliq, n0(3), f, g, domega_dt, depsil_dt, koef
      real*8 capa, e, inc, node, peri, capm
      real*8 tmp
      real*8 interpidx, dot_product3
      real*8 pm

c constants
      real*8 depsil_MAX
      parameter(depsil_MAX = 10.d0/DEGRAD)

      if (is_forward) then
        pm = 1.d0
      else
        pm = -1.d0
      endif

      do i = 1,ntp
        if (istat(i,1).eq.0) then

c calculate actual obliquity of a given body
          x(1) = xht(i)
          x(2) = yht(i)
          x(3) = zht(i)
          v(1) = vxht(i)
          v(2) = vyht(i)
          v(3) = vzht(i)
          
          call vector_product(x,v,n)
          call normalize(n)
          tmp = max(min(dot_product3(n, s(1,i)), 1.d0), -1.d0)
          obliq = acos(tmp)

c interpolation in f_ and g_func arrays
          call intep(obliq, f, obliq_func, f_func(0,fg_id(i)), nGDATA,
     :      IER)
          call intep(obliq, g, obliq_func, g_func(0,fg_id(i)), nGDATA,
     :      IER)
          
c get actual semimajor axis
          call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i),
     :      mass(1),ialpha,capa,e,inc,node,peri,capm)

c five scalings: with c_YORP, a, R, rho_bulk and also omega
c 2DO: scale by rho_surface too, due to thermal effects!
c (the day unit is here because of dt)

          koef = pm * c_YORP * (capa/a_0)**(-2) * (R(i)/R_0)**(-2)
     :      * (rho_bulk(i)/rho_0)**(-1) * day

          domega_dt = koef * f
          depsil_dt = koef * g/max(abs(omega(i)), 1.d-8)   ! prevent divisions by zeros!

c check depsil/dt acceleration (have to be smaller than ~10 deg per timestep)
          if (depsil_dt*dt.ge.depsil_MAX) then
c            write(*,*) 'yorp_evolve: Warning: depsil_dt = ',
c     :        depsil_dt*dt*DEGRAD, ' > ', depsil_MAX*DEGRAD,
c     :        ' deg per timestep for TPid = ', i
            depsil_dt = depsil_MAX/dt
          endif

c evolution equations (Capek & Vokrouhlicky 2004)

c change spin rate accordingly
          omega(i) = omega(i) + domega_dt*dt

c and rotate spin axis around (normal x spin)
          call vector_product(n, s(1,i), n_p)
          call normalize(n_p)
          call rotate_around_axis(n_p, depsil_dt*dt, s(1,i), s_new)
          call normalize(s_new)

          s(1,i) = s_new(1)
          s(2,i) = s_new(2)
          s(3,i) = s_new(3)

c          write(*,*) 'reorient: ', t/365.25d6, s(1,i), s(2,i), s(3,i),
c     :      omega(i), obliq, f, g, domega_dt, depsil_dt

        endif ! istat
      enddo   ! ntp

      return
      end


