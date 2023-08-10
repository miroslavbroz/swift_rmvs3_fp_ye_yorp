c  Rewrite from OrbFit by: Ondrej Chrenko (chrenko@sirrah.troja.mff.cuni.cz)
c  Cite: Yang et al. 2020 (Astronomy & Astrophysics, Volume 643, id.A38, 9 pp.)


      subroutine orbfit_forced(x,y,t,n,f,nf)

      include '../swift.inc'

c  input
      real*8 t(n),f(nf)
      integer n,nf
c  input & output
      real*8 x(n),y(n)

c  internals
      integer iremove_mode,i
      real*8 period,d0,d1,d2,d0k,d1k,d2k 

      iremove_mode = 1
      do i=1,nf
        period = 360.d0*3600.d0/f(i)
        call orbfit_peri2(t,y,n,period,iremove_mode,d0,d1,d2)
        call orbfit_peri2(t,x,n,period,iremove_mode,d0k,d1k,d2k)
      enddo

c      write(*,*) 'Ds: ',d0,d1,d2,d0k,d1k,d2k
c      stop

      return
      end

c--------------------------------------------------

      subroutine orbfit_argum(x,y,t,n,arg,freq)

      include '../swift.inc'

c  input
      real*8 x(n),y(n),t(n)
      integer n
c  output
      real*8 arg(n),freq

c  internals
      real*8 xm,ym
      real*8 nrevol(n)
      integer i

c  functions
      real*8 zero2pi

      xm = 0.d0
      ym = 0.d0
      do i=1,n
        xm = xm + x(i)
        ym = ym + y(i)
      enddo
      xm = xm/n
      ym = ym/n
      nrevol(1) = 0.d0
      arg(1) = atan2(y(1),x(1))
      do i=2,n
        arg(i) = atan2(y(i)-ym,x(i)-xm)
        arg(i) = zero2pi(arg(i))
        if (arg(i).gt.arg(i-1) + PI) then
          nrevol(i) = nrevol(i-1) - 1
        elseif (arg(i).lt.arg(i-1) - PI) then
          nrevol(i) = nrevol(i-1) + 1
        else
          nrevol(i) = nrevol(i-1)
        endif
      enddo
      do i=1,n
        arg(i) = arg(i) + TWOPI*dble(nrevol(i))
      enddo
      call linfit(t,arg,n,freq)

c      write(*,*) freq*3600.d0*180.d0/PI
c      stop      
 
      return
      end

c--------------------------------------------------

      subroutine orbfit_prop(t,x,y,n,period,amp,ph)

      include '../swift.inc'

c  input
      real*8 t(n),period
      integer n

c  input & output
      real*8 x(n),y(n)

c  output
      real*8 amp,ph      

c  internals
      real*8 d0,d1,d2,d0k,d1k,d2k
      integer iremove_mode

c  functions
      real*8 zero2pi

      iremove_mode = 1
      call orbfit_peri2(t,y,n,period,iremove_mode,d0,d1,d2)
      call orbfit_peri2(t,x,n,period,iremove_mode,d0k,d1k,d2k)
      amp = sqrt(d1*d1 + d2*d2)
      ph = zero2pi(atan2(d1,d2))

      return
      end

c--------------------------------------------------

      subroutine orbfit_peri2(t,dat,n,period,iremove_mode,d0,d1,d2)

      include '../swift.inc'

c  input 
      real*8 t(n),period
      integer n,iremove_mode
c  input & output
      real*8 dat(n)
c  output
      real*8 d0,d1,d2

c  internals
      integer i
      real*8 xk,sw,swc,swc2,sws,sws2,swcs,swcy,swsy,swy2,swy
      real*8 ym,tcur,x,dt,tol,xd,cd,sd,c,s
      real*8 a0,a1,a2,c1,c2,sp
      real*8 y(n),vc(n),vs(n)

      xk = TWOPI/period
      sw = 0.d0
      swc = 0.d0
      swc2 = 0.d0
      sws = 0.d0
      sws2 = 0.d0
      swcs = 0.d0
      swcy = 0.d0
      swsy = 0.d0
      swy2 = 0.d0
      swy = 0.d0
      do i=1,n
        swy = swy + dat(i)
      enddo
      sw = dble(n)
      ym = swy/sw
c      write(*,*) 'ym: ',ym
      do i=1,n
        y(i) = dat(i) - ym
      enddo
      tcur = t(1)
c      write(*,*) 'tcur first: ',tcur
      x = xk*t(1)
      vc(1) = cos(x)
      vs(1) = sin(x)
      dt = t(2) - t(1)
      tol = abs(dt * 1.d-4)
      xd = xk*dt
      cd = cos(xd)
      sd = sin(xd)
      do i=1,n
        if (abs(t(i) - tcur).gt.tol) then
          tcur = t(i)
          x = tcur*xk
          vc(i) = cos(x)
          vs(i) = sin(x)
        endif
        c = vc(i)
        s = vs(i)
        swc = swc + c
        swc2 = swc2 + c**2
        sws = sws + s
        sws2 = sws2 + s**2
        swcs = swcs + c*s
        swcy = swcy + c*y(i)
        swsy = swsy + s*y(i)
        swy2 = swy2 + y(i)**2
        if (i.lt.n) then
          vc(i+1) = c*cd - s*sd
          vs(i+1) = s*cd + c*sd
          tcur = tcur + dt
        endif
c        write(*,*) 'after loop 1'
c        write(*,*) c,s,swc,swc2,sws,sws2,swcs
c        write(*,*) swcy,swsy,swy2,vc(i+1),vs(i+1),tcur
c        stop
      enddo
      a0 = 1.d0/sqrt(sw)
      a1 = 1.d0/sqrt(swc2 - (a0*swc)**2)
      a2 = sws2 - (a0*sws)**2 - (a1*swcs)**2 - (a0*a0*a1*swc*sws)**2
      a2 = a2 + 2.d0*swc*sws*swcs*(a0*a1)**2
      a2 = 1.d0/sqrt(a2)
      c1 = a1*swcy
      c2 = a2*swsy - a1*a2*c1*(swcs - a0*a0*swc*sws)
      sp = c1*c1 + c2*c2
      sp = sp/swy2
      d2 = a2*c2
      d1 = a1*c1 + a2*a1*a1*c2*(a0*a0*swc*sws - swcs)
      d0 = -a0*a0*(d1*swc + d2*sws) + ym
c      write(*,*) 'Stuff 1'
c      write(*,*) c,s,swc,swc2,sws
c      write(*,*) 'Stuff 2'
c      write(*,*) sws2,swcs,swcy,swsy,swy2
c      write(*,*) 'Stuff 3'
c      write(*,*) a0,a1,a2,c1,c2,sp,d2,d1,d0
      if (iremove_mode.eq.0) then
        return
      endif
      do i=1,n
        dat(i) = y(i) - (d0 + d1*vc(i) + d2*vs(i))
      enddo

      return
      end

c--------------------------------------------------

      subroutine linfit(x,y,n,slope)

      include '../swift.inc'

c  input
      real*8 x(n),y(n)
      integer n

c  output
      real*8 slope

c  internals
      real*8 tm,t,det,tmp,cost
      real*8 c(2),b(2,2)
      integer i,j

      tm = 0.d0
      do i=1,n
        tm = tm + x(i)
      enddo
      tm = tm/n
      do i=1,2
        do j=1,2
          b(i,j) = 0.d0
        enddo
        c(i) = 0.d0
      enddo
      do i=1,n
        t = x(i) - tm
        b(1,2) = b(1,2) + t
        b(2,2) = b(2,2) + t*t
        c(1) = c(1) + y(i)
        c(2) = c(2) + y(i)*t
      enddo
      b(1,1) = n
      det = b(1,1)*b(2,2) - b(1,2)*b(1,2)
      b(2,1) = -b(1,2)/det
      b(1,2) = b(2,1)
      tmp = b(1,1)/det
      b(1,1) = b(2,2)/det
      b(2,2) = tmp
      cost = b(1,1)*c(1) + b(1,2)*c(2)
      slope = b(2,1)*c(1) + b(2,2)*c(2)

      return
      end
