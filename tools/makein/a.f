c  a.f
c  astronomical functions and subroutines library
c  miroslav broz (miroslav.broz@usa.net)
c  12. XI. 1997
c  thanks to martin solc (solc@mbox.cesnet.cz),
c    zdenek moravec (klet@klet.cz)
c
c  references:
c  Wolf, M. et al.: Atronomicka prirucka. Academia, Praha, 1992
c  Andrle, P.: Zaklady nebeske mechaniky. Academia, Praha, 1971
c  Dubyago, A. D.: The Determination of Orbits, The Macmillan
c    Company, New York, 1961 
c  Prihoda, P.: Astronomicke algoritmy pro kalkulatory.
c

c**********************************************************************
      double precision function sgn(x)
c**********************************************************************
c
c  signum of x
c      
      if (x.gt.0d0) then 
        sgn=1d0
        return
      elseif (x.lt.0d0) then
        sgn=-1d0
        return
      else
        sgn=0d0
      endif    
      end

c**********************************************************************
      double precision function frac(x)
c**********************************************************************
c
c  fractional part of x
c
      double precision x
      frac=x-int(x)
      return
      end
 
c**********************************************************************
      double precision function nula2pi(x)
c**********************************************************************
c
c  <0,2pi) truncation
c
      double precision x,frac,y,pi2
      parameter(pi2=6.283185307d0)
      y=frac(x/pi2)*pi2
      if (y.lt.0d0) y=y+pi2
      nula2pi=y
      return
      end 

c**********************************************************************
      subroutine hhms(h2,h,m,s)
c**********************************************************************
c
c  conversion of decimal hours to hours, minutes and seconds
c
      double precision h2,h,m,s
      h=int(h2)
      m=int((h2-h)*60)
      s=(h2-h)*3600-m*60
      return
      end 

c**********************************************************************
      double precision function hmsh(h,m,s)
c**********************************************************************
c
c  conversion of hours, minutes and seconds to decimal hours
c
      double precision h,m,s
      hmsh=h+m/60+s/3600
      return
      end

c**********************************************************************
      double precision function jd(y,m,d)
c**********************************************************************
c
c  julian date
c  y  year
c  m  month
c  d  day
c      
      double precision y,m,d
      double precision yy,mm,a,b,jd1
       
      if (m.gt.2) then
        yy=y
        mm=m
      else
        yy=y-1
        mm=m+12
      endif
      jd1=int(365.25*yy)+int(30.6001*(mm+1))+d+1720994.5d0
      if ((int(y).gt.1582).or.((int(y).eq.1582).and.(int(m).gt.10))
     &.or.((int(y).eq.1582).and.(int(m).eq.10).and.(int(d).gt.15)))
     &then
        a=int(yy/100)
        b=2-a+int(a/4)
        jd1=jd1+b
      endif
      jd=jd1
      return
      end

c**********************************************************************
      subroutine lst(jd,lambda,s0,ss,s)
c**********************************************************************
c
c  derivation of the local sidereal time from julian date an longitude
c  jd      julian date [d]
c  lambda  geocentric longitude [deg]
c  s0      mean sidereal time at 0 ut [h]
c  ss      mean sidereal time [h]
c  s       local sidereal time [h]   
c
      double precision jd,lambda,s0,ss,s
      double precision jd0,t,t3,frac
      jd0=aint(jd+0.5)-0.5
      t=24*(jd-jd0)
      t3=(jd0-2451545.0)/36525
      s0=6.697374558+2400.05133691*t3+0.0000258622*t3**2-1.7e-9*t3**3
      s0=24.*frac(s0/24.)
      if (s0.lt.0) s0=s0+24
      ss=s0+1.002737909350795*t
      s=ss+lambda/15.
      return
      end

c**********************************************************************
      integer function max(x,y,z)
c**********************************************************************
c
c  integer index of maximal double precision number
c  1 corresponds to x, ...  
c
      double precision x,y,z
      if (x.gt.y) then
        if (x.gt.z) then
          max=1
        else
          max=3
        endif
      else
        if (y.gt.z) then
          max=2
        else
          max=3
        endif
      endif
      return
      end

c**********************************************************************
      subroutine e3(t,alpha,delta,x0,y0,z0,elmts)
c**********************************************************************
c
c  determination of an orbit by three observations 
c  (Gauss-Lagrange method)
c
c  t(1:3)             julian date of observation [d]
c  alpha(), delta()   topocentric equatoreal coordinates [rad]
c  x0(1:3),y0(),z0()  rectangular equatoreal coordinates of sun [au]
c  elmts(1:6)         elements of orbit
c
c  parameter in common block /equinox/:
c  t0      equinox for elements, julian date [d]
c
c  parameters in common block /observer/:
c  lambda  geocentric longitude of the observatory
c  fi      geocentric latitude
c  height  height above sea level
c
      implicit double precision (a-z)
      parameter(eps=1d-13)

      integer i
      double precision pi,deg,rad,k_gauss,aber,eps_earth,r_earth,epsilon
      double precision t(3),alpha(3),delta(3),x(3),y(3),z(3),elmts(6)
      double precision h(3),a(3),b(3),c(3),d(3),r(3),rho(3)
      double precision rh(3),xh(3),yh(3),zh(3),tc(3),s(3),eta(3)
      double precision x0(3),y0(3),z0(3)
      double precision t0,lambda,fi,height
      common/equinox/t0
      common/observer/lambda,fi,height

      data pi,deg,rad/3.1415926535d0,0.017453292d0,57.29577951d0/
      data k_gauss,aber/0.01720209895d0,0.005706d0/
      data eps_earth,r_earth/0.003352329869d0,6378.140d3/
      data epsilon/0.409099609d0/
c
c  parallax correction
c
      do i=1,3
        call lst(t(i),lambda*rad,s0,ss,h(i))
        h(i)=nula2pi(h(i)/12*pi)
      enddo 

      e=sqrt(2*eps_earth-eps_earth**2)
      ss=(1-e**2)/sqrt(1-e**2*sin(fi)**2)
      cc=1/sqrt(1-e**2*sin(fi)**2)
      phi=atan((ss+height/r_earth)/(cc+height/r_earth)*tan(fi))
      rho0=(cc+height/r_earth)*cos(phi)/cos(fi)

      deltaxy=-sin(8.8/3600.*deg)*rho0*cos(fi)
      deltaz=-sin(8.8/3600.*deg)*rho0*sin(fi)
      deltaxy=0
      deltaz=0

      do i=1,3
        x(i)=x0(i)+deltaxy*cos(h(i))
        y(i)=y0(i)+deltaxy*sin(h(i))
        z(i)=z0(i)+deltaz
      enddo
c
c  auxiliary quantities
c
      do i=1,3
        a(i)=cos(delta(i))*cos(alpha(i))
        b(i)=cos(delta(i))*sin(alpha(i))
        c(i)=sin(delta(i))
      enddo
      a31=b(3)*c(1)-b(1)*c(3)
      b31=a(1)*c(3)-a(3)*c(1)
      c31=a(3)*b(1)-a(1)*b(3)
      d0=a(2)*a31+b(2)*b31+c(2)*c31
      if (abs(d0).lt.eps) then
        write(*,10)
10      format(' The orbit must be determined by four observations.')
      endif

      do i=1,3
        d(i)=x(i)*a31+y(i)*b31+z(i)*c31
        r(i)=sqrt(x(i)**2+y(i)**2+z(i)**2)
      enddo
      c0=-(a(2)*x(2)+b(2)*y(2)+c(2)*z(2))
c
c  zero approximation of geocentric distance
c
      tau1=k_gauss*(t(3)-t(2))
      tau2=k_gauss*(t(3)-t(1))
      tau3=k_gauss*(t(2)-t(1))

      n10=tau1/tau2
      n30=tau3/tau2

      nu1=1./6.*tau1*tau3*(1.+n10)
      nu3=1./6.*tau1*tau3*(1.+n30)
      p=-(n10*d(1)-d(2)+n30*d(3))/d0
      rho(2)=p

      n1=0.
      n3=0.

20    n11=n1
      n33=n3
c
c  solution of Lagranges equotions
c
      q=(nu1*d(1)+nu3*d(3))/d0
30    rh(2)=sqrt(r(2)**2+2*c0*rho(2)+rho(2)**2)
      f=rho(2)-p+q/rh(2)**3.
      f2=1.d0-3.d0*q*(rho(2)+c0)/rh(2)**5.
      rho2=rho(2)
      rho(2)=rho(2)-f/f2
      if (abs(rho(2)-rho2).lt.eps) goto 35
      goto 30
35    continue

      if (rho(2).lt.0) then
        write(*,40)
40      format(' Observations does not correspond with circular orbit.')
      endif
      rh(2)=sqrt(r(2)**2+2*c0*rho(2)+rho(2)**2)
      n1=n10+nu1/rh(2)**3
      n3=n30+nu3/rh(2)**3
c
c  derivation of rho1, rho3
c
      i=max(a31,b31,c31)
      if (i.eq.1) then
        rho(1)=-((n1*y(1)-y(2)+n3*y(3))*c(3)-(n1*z(1)-z(2)+
     &  n3*z(3))*b(3)+(b(2)*c(3)-c(2)*b(3))*rho(2))/(a31*n1)
      elseif (i.eq.2) then
        rho(1)=((n1*x(1)-x(2)+n3*x(3))*c(3)-(n1*z(1)-z(2)+
     &  n3*z(3))*a(3)+(a(2)*c(3)-c(2)*a(3))*rho(2))/(b31*n1)
      else
        rho(1)=-((n1*x(1)-x(2)+n3*x(3))*b(3)-(n1*y(1)-y(2)+
     &  n3*y(3))*a(3)+(a(2)*b(3)-b(2)*a(3))*rho(2))/(c31*n1)
      endif

      i=max(a(3),b(3),c(3))
      if (i.eq.1) then
        rho(3)=(a(2)*rho(2)-a(1)*n1*rho(1)+n1*x(1)-x(2)+
     &  n3*x(3))/(n3*a(3))
      elseif (i.eq.2) then
        rho(3)=(b(2)*rho(2)-b(1)*n1*rho(1)+n1*y(1)-y(2)+
     &  n3*y(3))/(n3*b(3))
      else
        rho(3)=(c(2)*rho(2)-c(1)*n1*rho(1)+n1*z(1)-z(2)+
     &  n3*z(3))/(n3*c(3))
      endif
c
c  heliocentric coordinates
c
      do i=1,3
        xh(i)=a(i)*rho(i)-x(i)
        yh(i)=b(i)*rho(i)-y(i)
        zh(i)=c(i)*rho(i)-z(i)
        rh(i)=sqrt(xh(i)**2+yh(i)**2+zh(i)**2)
      enddo
c
c  correction for aberration
c
      do i=1,3 
        tc(i)=t(i)-aber*rho(i)
      enddo
      tau1=k_gauss*(tc(3)-tc(2))
      tau2=k_gauss*(tc(3)-tc(1))
      tau3=k_gauss*(tc(2)-tc(1))
      n10=tau1/tau2
      n30=tau3/tau2
c
c  the ratio of the area of the sector to the area of the triangle
c
      k1=sqrt(2*(rh(2)*rh(3)+xh(2)*xh(3)+yh(2)*yh(3)+zh(2)*zh(3)))
      k2=sqrt(2*(rh(1)*rh(3)+xh(1)*xh(3)+yh(1)*yh(3)+zh(1)*zh(3)))
      k3=sqrt(2*(rh(1)*rh(2)+xh(1)*xh(2)+yh(1)*yh(2)+zh(1)*zh(2)))

      h(1)=tau1**2/(k1**2*(k1/3+(rh(2)+rh(3))/2))
      h(2)=tau2**2/(k2**2*(k2/3+(rh(1)+rh(3))/2))
      h(3)=tau3**2/(k3**2*(k3/3+(rh(1)+rh(2))/2))

      s(1)=1
      s(2)=1
      s(3)=1

50    s1=s(1)
      s2=s(2)
      s3=s(3)
      do i=1,3
        s(i)=1.d0+11.d0/9.d0*h(i)/s(i)
      enddo
      if ((abs(s(1)-s1).lt.eps).and.(abs(s(2)-s2).lt.eps).and.
     &(abs(s(3)-s3).lt.eps)) goto 55
      goto 50
55    continue

      do i=1,3
        eta(i)=1.+10./9.*h(i)/s(i)
      enddo
c
c  determination of new values of n1, n3 and nu1, nu3
c
      n1=n10*eta(2)/eta(1)
      n3=n30*eta(2)/eta(3)
      nu1=(n1-n10)*rh(2)**3
      nu3=(n3-n30)*rh(2)**3

      if ((abs(n1-n11).lt.eps).and.(abs(n3-n33).lt.eps)) goto 25
      goto 20
25    continue
c
c  derivation of elements
c
      sigma=(xh(1)*xh(3)+yh(1)*yh(3)+zh(1)*zh(3))/rh(1)**2
      x00=xh(3)-sigma*xh(1)
      y00=yh(3)-sigma*yh(1)
      z00=zh(3)-sigma*zh(1)
      r0=sqrt(x00**2+y00**2+z00**2)

      sinf=r0/rh(3)
      cosf=sigma*rh(1)/rh(3)
      f=atan2(r0,(sigma*rh(1)))
      p=(rh(1)*rh(3)*sinf*eta(2)/tau2)**2
      q1=p/rh(1)-1
      q3=p/rh(3)-1
      v1=atan2(cosf-q3/q1,sinf)
      v3=v1+f
      e0=q1/cos(v1)

      if (abs(e0).gt.1) then
        write(*,100)
100     format(1x,'Impossible to derivate eliptical elements.')
      endif

      a0=p/(1-e0**2)
      ee1=2*atan(sqrt((1-e0)/(1+e0))*tan(v1/2))
      ee3=2*atan(sqrt((1-e0)/(1+e0))*tan(v3/2))
      m1=ee1-e0*sin(ee1)
      m3=ee3-e0*sin(ee3)
      n=(m3-m1)/(t(3)-t(1))
      m0=m1+n*(t0-tc(1))
      p=a0**1.5

      px=xh(1)/rh(1)*cos(v1)-x00/r0*sin(v1)
      py=yh(1)/rh(1)*cos(v1)-y00/r0*sin(v1)
      pz=zh(1)/rh(1)*cos(v1)-z00/r0*sin(v1)
      qx=xh(1)/rh(1)*sin(v1)+x00/r0*cos(v1)
      qy=yh(1)/rh(1)*sin(v1)+y00/r0*cos(v1)
      qz=zh(1)/rh(1)*sin(v1)+z00/r0*cos(v1)

      omega0=nula2pi(atan2(pz*cos(epsilon)-py*sin(epsilon),
     &qz*cos(epsilon)-qy*sin(epsilon)))
      omega=nula2pi(atan2(py*cos(omega0)-qy*sin(omega0),
     &(px*cos(omega0)-qx*sin(omega0))*cos(epsilon)))
      i0=atan(-((pz*cos(epsilon)-py*sin(epsilon))*sin(omega))/
     &((px*sin(omega0)+qx*cos(omega0))*sin(omega0)))

      if (e0.lt.0d0) then
        omega0=nula2pi(omega0+pi)
        m0=nula2pi(m0+pi)
      endif

      elmts(1)=a0
      elmts(2)=abs(e0)
      elmts(3)=i0
      elmts(4)=omega0
      elmts(5)=omega
      elmts(6)=m0
      
      end
      
c**********************************************************************
      double precision function kepler(m,e)
c**********************************************************************
c
c  keplers equotion, eccentric anomaly, simple iteration method 
c  m  mean anomaly
c  e  eccentricity
c    
      double precision m,e
      double precision e1,e2,eps
      parameter(eps=1d-13) 
      e2=m+e*sin(m)
10    e1=e2
      e2=m+e*sin(e1)
      if (abs(e2-e1).lt.eps) goto 15
      goto 10
15    kepler=e2
      return 
      end

c**********************************************************************
      double precision function wegstein(m,e)
c**********************************************************************
c
c  keplers equotion, eccentric anomaly, wegstein's method as in
c  Taff, L. G., Brennan, T. A.: On solving Kepler's equotion,
c  STSCI preprint series No. 371, September 1989
c  m  mean anomaly
c  e  eccentricity
c
      double precision e,m,e1,e2,e3,y1,y2,y3,eps
      parameter(eps=1d-13) 
      e1=m+e
      y1=m+e*sin(e1)
      e2=y1
      y2=m+e*sin(e2)
1     e3=e2+(e2-e1)/((e1-y1)/(e2-y2)-1)
      y3=m+e*sin(e3)
      e1=e2
      e2=e3
      y1=y2
      y2=y3
      if (abs(e2-e1).ge.eps) goto 1
      wegstein=e2
      return
      end

c**********************************************************************
      double precision function ekepl1(m,e)
c**********************************************************************
c       
c  keplers equotion, eccentric anomaly, legendre-based starter and
c  halley iterator, method described in
c  Odell, A. W., Gooding, R. H.: Procedures for solving Kepler's
c  equotion. Tech. Memo Space 356, p. 35, January 1986
c  m  mean anomaly
c  e  eccentricity
c 
      implicit none
      double precision m,e,eps,c,s,psi,psi0,xi,eta,fd,fdd,f
      parameter(eps=1d-13) 
      c=e*cos(m)
      s=e*sin(m)
      psi=s/sqrt(1d0-c-c+e*e)
1     xi=cos(psi)
      eta=sin(psi)
      fd=(1d0-c*xi)+s*eta
      fdd=c*eta+s*xi
      f=psi-fdd
      psi0=psi
      psi=psi-f*fd/(fd*fd-5d-1*f*fdd)
      if (abs(psi-psi0).ge.eps) goto 1
      ekepl1=m+psi
      return
      end

c**********************************************************************
      double precision function emkepl(m,e)
c**********************************************************************
c       
c  keplers equotion, eccentric anomaly, method described in
c  Odell, A. W., Gooding, R. H.: Procedures for solving Kepler's
c  equotion. Tech. Memo Space 356, p. 37, January 1986
c  m  mean anomaly
c  e  eccentricity
c  
      double precision m,e,x,ee2,term,d,x0,eps
      parameter(eps=1d-13) 
      x=(1-e)*sin(m)
      ee2=-m*m
      term=m
      d=0d0
1     d=d+2d0
      term=term*ee2/(d*(d+1d0))
      x0=x
      x=x-term
      if (abs(x-x0).ge.eps) goto 1
      emkepl=x
      return
      end

c
c  rotation matrices
c
c**********************************************************************
      subroutine rotat1(x0,x,a)
c**********************************************************************
      double precision x0(3),x(3),a,c,s
      c=cos(a)
      s=sin(a)
      x(1)=x0(1)
      x(2)=x0(2)*c-x0(3)*s
      x(3)=x0(2)*s+x0(3)*c
      end

c**********************************************************************
      subroutine rotat2(x0,x,a)
c**********************************************************************
      double precision x0(3),x(3),a,c,s
      c=cos(a)
      s=sin(a)
      x(1)=x0(1)*c-x0(3)*s
      x(2)=x0(2)
      x(3)=x0(1)*s+x0(3)*c
      end

c**********************************************************************
      subroutine rotat3(x0,x,a)
c**********************************************************************
      double precision x0(3),x(3),a,c,s
      c=cos(a)
      s=sin(a)
      x(1)=x0(1)*c-x0(2)*s
      x(2)=x0(1)*s+x0(2)*c
      x(3)=x0(3)
      end

c
c  mean longitudes of planets
c
      double precision function lm1(t)
      double precision t
      lm1=4.4026d0+2608.7903d0*t
      end
      double precision function lm2(t)
      double precision t
      lm2=3.1761d0+1021.3286d0*t
      end
      double precision function lm3(t)
      double precision t
      lm3=1.7535d0+628.3076d0*t
      end
      double precision function lm4(t)
      double precision t
      lm4=6.2035d0+334.0612d0*t
      end
      double precision function lm5(t)
      double precision t
      lm5=0.5995d0+52.9691d0*t
      end
      double precision function lm6(t)
      double precision t
      lm6=0.8740d0+21.3299d0*t
      end
      double precision function lm7(t)
      double precision t
      lm7=5.4813d0+7.4782d0*t
      end
      double precision function lm8(t)
      double precision t
      lm8=5.3119d0+3.8133d0*t
      end

c**********************************************************************
      subroutine elmts_earth(t,a,l,k,h,q,p)
c**********************************************************************
c
c  elements of earth
c
      parameter(maxi=17)
      double precision t,a,l,k,h,q,p
      double precision l2,l3,l4,l5,arg(maxi),s(maxi),c(maxi)
      double precision lm2,lm3,lm4,lm5
 
      l2=lm2(t)
      l3=lm3(t)
      l4=lm4(t)
      l5=lm5(t)

      arg(1)=4*l3-8*l4+3*l5
      arg(2)=2*l2-3*l3
      arg(3)=l3-2*l5
      arg(4)=l5
      arg(5)=2*l3-2*l5 
      arg(6)=l2-l3 
      arg(7)=l2-2*l3 
      arg(8)=2*l2-2*l3
      arg(9)=l2 
      arg(10)=3*l2-4*l3
      arg(11)=l3-2*l4 
      arg(12)=8*l2-13*l3 
      arg(13)=l3-3*l5 
      arg(14)=3*l2-3*l3 
      arg(15)=l3-l5 
      arg(16)=3*l2-5*l3 
      arg(17)=2*l5

      do i=1,maxi
        s(i)=sin(arg(i))
        c(i)=cos(arg(i))
      enddo
     
      a=1.0000010+
     &1d-7*(15*s(2)+6*s(3)-3*s(7)+5*s(10)-3*s(11)-s(16)
     &+3*c(3)+112*c(5)+76*c(6)+c(7)-41*c(8)-3*c(11)-25*c(14)+11*c(15)
     &-4*c(16))

      l=1.7534703+628.3075849*t-0.0000001*t**2+
     &1d-7*(-97*s(1)-2*s(2)-13*s(3)-125*s(4)-206*s(5)+166*s(6)+3*s(7)
     &+127*s(8)-65*s(11)+49*s(12)-3*s(13)+60*s(14)-51*s(15)-46*s(16)
     &-4*s(17)
     &+322*c(1)+109*c(2)+20*c(3)+17*c(4)-c(7)+14*c(10)+56*c(11)+76*c(12)
     &+4*c(13)-c(15)+12*c(16)+c(17))

      k=-0.0037408-0.0000823*t+0.0000003*t**2+
     &1d-7*(-3*s(4)-2*s(6)+s(11)-s(12)-6*s(13)-s(15)-29*s(16)+22*s(17)
     &-199*c(2)+186*c(3)-150*c(4)+82*c(7)+48*c(9)-48*c(10)-44*c(11)
     &-2*c(12)+35*c(13)-7*c(15)-c(16)-16*c(17))
    
      h=0.0162845-0.0000620*t-0.0000003*t**2+
     &1d-7*(199*s(2)-186*s(3)-151*s(4)-82*s(7)+48*s(9)+48*s(10)
     &+45*s(11)+2*s(12)-35*s(13)+s(15)+s(16)-17*s(17)
     &+c(4)-c(12)-6*c(13)-2*c(15)-29*c(16)-22*c(17))

      q=-0.0001135d0*t+0.0000001d0*t**2
      p=0.0000102d0*t+0.0000005d0*t**2
      end
  
c**********************************************************************
      double precision function eps_earth(t)
c**********************************************************************
c
c  inclination of ecliptic to equator
c      
      double precision pi,deg
      parameter(pi = 3.14159265358979d0, deg=pi/180.d0)
      double precision t,hmsh
      eps_earth=hmsh(23d0,26d0,
     &21.448d0-46.8150d0*t-0.00059d0*t**2+0.001813d0*t**3)*deg
      return
      end
      
c**********************************************************************
      subroutine preces_angle(t0,t,dzeta,zz,fi)
c**********************************************************************
c
c  precession angles
c      
      double precision k
      parameter(k=4.848136811d-6)
      double precision t0,t,dzeta,zz,fi
      dzeta=(2306.2181d0+1.39656d0*t0-0.00139d0*t0**2)*t
     &+(0.30188d0-0.000344d0*t0)*t**2+0.017998d0*t**3
      zz=(2306.2181d0+1.39656d0*t0-0.00139d0*t0**2)*t
     &+(1.09468d0+0.000066d0*t0)*t**2+0.018203d0*t**3
      fi=(2004.3109d0-0.85330d0*t0-0.000217d0*t0**2)*t
     &-(0.42665d0+0.000217d0*t0)*t**2-0.041833d0*t**3
      dzeta=k*dzeta
      zz=k*zz
      fi=k*fi
      return
      end

c**********************************************************************
      subroutine preces(x,dzeta,zz,fi)
c**********************************************************************
c
c  correction for precession
c      
      double precision x(3),dzeta,zz,fi
      double precision x1(3),x2(3)
      call rotat3(x,x1,-dzeta)
      call rotat2(x1,x2,fi)
      call rotat3(x2,x,-zz)
      return
      end
      
c**********************************************************************
      subroutine nutate_angle(t,deltapsi,deltaepsilon)
c**********************************************************************
c
c  nutation angles
c      
      double precision k
      parameter(k=4.848136811d-6)
      double precision t,deltapsi,deltaepsilon
      double precision l,l2,f,d,omega
      common/nutate_anglepar/d
      l=485866.733+1717915922.633*t+31.310*t**2+0.064*t**3
      l2=1287099.804+129596518.224*t-0.577*t**2-0.012*t**3
      f=335778.667+1739527263.137*t-13.257*t**2+0.011*t**3
      d=1072261.3077+1602961601.328*t-6.891*t**2+0.019*t**3
      omega=450160.280-6962890.539*t+0.455*t**2+0.008*t**3
      l=k*l
      l2=k*l2
      f=k*f
      d=k*d
      omega=k*omega
      deltapsi=-17.2*sin(omega)-1.32*sin(2*f-2*d+2*omega)
     1-0.23*sin(2*f+2*omega)+0.21*sin(2*omega)+0.14*sin(l2)+0.07*sin(l)
     2-0.05*sin(l2+2*f-2*d+2*omega)-0.04*sin(2*f+omega)
     3-0.03*sin(l+2*f+2*omega)
      deltaepsilon=9.2*cos(omega)+0.57*cos(2*f-2*d+2*omega)
     1+0.10*cos(2*f+2*omega)-0.09*cos(2*omega)+0.01*cos(l2)
     2+0.02*cos(l2+2*f-2*d+2*omega)+0.02*cos(2*f+omega)
     3+0.01*cos(l+2*f+2*omega)
      deltapsi=k*deltapsi
      deltaepsilon=k*deltaepsilon
      return
      end

c**********************************************************************
      subroutine nutate(x,epsilon,deltapsi,deltaepsilon)
c**********************************************************************
c
c  correction for nutation
c      
      double precision x(3),epsilon,deltapsi,deltaepsilon
      double precision x1(3),x2(3)
      call rotat1(x,x1,epsilon)
      call rotat3(x1,x2,-deltapsi)
      call rotat1(x2,x,-epsilon-deltaepsilon)
      return
      end

c**********************************************************************
      subroutine elmts_planet(planet,t,elmts)
c**********************************************************************
c
c  elements of planets, derivation of keplers elements
c
      integer planet
      double precision t
      double precision elmts(6)
      double precision a,l,k,h,q,p
      if (planet.eq.1) then
      elseif (planet.eq.2) then
      elseif (planet.eq.3) then
        call elmts_earth(t,a,l,k,h,q,p)
      elseif (planet.eq.4) then
      elseif (planet.eq.5) then
      elseif (planet.eq.6) then
      elseif (planet.eq.7) then
      elseif (planet.eq.8) then
      endif

      elmts(1)=a 
      elmts(2)=sqrt(k*k+h*h)
      elmts(3)=2*asin(sqrt(q*q+p*p))
      elmts(4)=atan2(h,k)
      elmts(5)=atan2(p,q)
      elmts(6)=l-elmts(4)
      elmts(4)=elmts(4)-elmts(5)
      return
      end
 
c**********************************************************************
      subroutine eph(elmts,x)
c**********************************************************************
c
c  calculation of an ephemeris
c  elmts(1:6)  elements of orbit
c  x(1:3)      rectangular coordinates
c
      double precision elmts(6),x(3)
      double precision e,kepler
      double precision x0(3)
      e=kepler(elmts(6),elmts(2))
      x0(1)=elmts(1)*(cos(e)-elmts(2))
      x0(2)=elmts(1)*sqrt(1-elmts(2)**2)*sin(e)
      x0(3)=0d0
      call rotat3(x0,x,elmts(4))
      call rotat1(x,x0,elmts(3))
      call rotat3(x0,x,elmts(5))
c      
c  calculation with right anomaly      
c
c      v=2*atan(sqrt((1+elmts(2))/(1-elmts(2)))*tan(e/2))
c      u=v+elmts(4)
c      l=elmts(5)+atan2(cos(elmts(3))*sin(u),cos(u))
c      b=asin(sin(elmts(3))*sin(u))
c      r=elmts(1)*(1-elmts(2)*cos(e))
c      x(1)=r*cos(b)*cos(l)
c      x(2)=r*cos(b)*sin(l)
c      x(3)=r*sin(b)
      return
      end
      
c**********************************************************************
      subroutine e2(t,rho,alpha,delta,x,y,z,elmts)
c**********************************************************************
c
c  determination of an orbit by two observations
c
c  t(1:2)                 julian date of observation [d]            
c  rho(),alpha(),delta()  geocentric equatoreal coordinates [rad]
c  x(),y(),z()            rectangular equatoreal coords of sun [rad]
c  elmts(1:6)             elements of orbit
c
c  parameter in common block /equinox/:
c  t0  equinox for elements, julian date [d]
c
c  parameters in common block /observer/:
c  lambda  geocentric longitude of the observatory
c  fi      geocentric latitude
c  height  height above sea level
c
      implicit double precision (a-z) 
      integer i
      double precision t(2),rho(2),alpha(2),delta(2),x(2),y(2),z(2)
      double precision elmts(6)
      double precision x00,y00,z00,x0(2),y0(2),z0(2),xx(3)
      double precision r1,r2,a,b,c,sinf,sini,tau,k,h,s,s0,eta,q1,q2
      double precision tc1,tc2,p,e,ee1,ee2,m1,m2,n,v1,v2,f
      double precision eps,k_gauss,aber,epsilon
      double precision t0,lambda,fi,height
      parameter(eps=1d-13)
      parameter(epsilon=0.409099609d0)
      parameter(k_gauss=0.01720209895d0,aber=0.005706d0)
      parameter(pi=3.1415926535d0)
      common/equinox/t0
      common/observer/lambda,fi,height
c
c  transformation for rectagular ecliptical coords
c
      do i=1,2
c        call parcor(t(i),xx)
        xx(1)=0
        xx(2)=0
        xx(3)=0
        x00=rho(i)*cos(alpha(i))*cos(delta(1))-x(i)-xx(1)
        y00=rho(i)*sin(alpha(i))*cos(delta(1))-y(i)-xx(2)
        z00=rho(i)*sin(delta(i))-z(i)-xx(3)
        x0(i)=x00
        y0(i)=y00*cos(epsilon)+z00*sin(epsilon)
        z0(i)=-y00*sin(epsilon)+z00*cos(epsilon)
      enddo  
      r1=sqrt(x0(1)**2+y0(1)**2+z0(1)**2)
      r2=sqrt(x0(2)**2+y0(2)**2+z0(2)**2)
c
c  correction for aberration
c
      tc1=t(1)-aber*rho(1)
      tc2=t(2)-aber*rho(2)
      tau=k_gauss*(tc2-tc1)
c
c  derivation of f, i, omega
c
      a=(y0(1)*z0(2)-y0(2)*z0(1))/(r1*r2)
      b=(z0(2)*x0(1)-z0(1)*x0(2))/(r1*r2)
      c=(x0(1)*y0(2)-x0(2)*y0(1))/(r1*r2)
      sinf=sqrt(a**2+b**2+c**2)
      f=asin(sinf)
      sini=sqrt(a**2+b**2)/sinf
      elmts(3)=asin(sini)
      elmts(5)=nula2pi(atan2(a,b))
c
c  the ratio of the area of the sector to the area of the triangle
c
      k=sqrt(2*(r1*r2+x0(1)*x0(2)+y0(1)*y0(2)+z0(1)*z0(2)))
      h=tau**2/(k**2*(k/3+(r1+r2)/2))
      s=1
10    s0=s
      s=1.d0+11.d0/9.d0*h/s
      if (abs(s-s0).lt.eps) goto 15
      goto 10
15    continue
      eta=1.d0+10.d0/9.d0*h/s
c
c  elements a, e, omega, m 
c     
      p=(r1*r2*sinf*eta/tau)**2
      q1=p/r1-1
      q2=p/r2-1
      v1=atan2(cos(f)-q2/q1,sinf)
      v2=v1+f
      e=q1/cos(v1)
      elmts(2)=abs(e)
      elmts(1)=p/(1-e**2)
      ee1=2*atan(sqrt((1-e)/(1+e))*tan(v1/2))
      ee2=2*atan(sqrt((1-e)/(1+e))*tan(v2/2))
      m1=ee1-e*sin(ee1)
      m2=ee2-e*sin(ee2)
      n=(m2-m1)/(t(2)-t(1))
      elmts(6)=m1+n*(t0-tc1)
      elmts(4)=nula2pi(atan2(z0(1),(x0(1)*b+y0(1)*a)/sinf)-v1)
      if (e.lt.0) then
        elmts(4)=nula2pi(elmts(4)+pi)
        elmts(6)=nula2pi(elmts(6)+pi)
      endif
      end

c**********************************************************************
      subroutine parcor(t,dx)
c**********************************************************************
c
c  correction for parallax
c
c  t      julian date [d]
c  dx(3)  corrections of equatoreal rectangular coordinates
c
c  parameters in common block /observer/:
c  lambda  geocentric longitude of the observatory
c  fi      geocentric latitude
c  height  height above sea level
c
      double precision t,dx(3)
      double precision h,nula2pi,s0,ss,e,cc,phi,rho0,deltaxy,deltaz
      double precision pi,deg,rad,eps_earth,r_earth
      double precision lambda,fi,height
      common/observer/lambda,fi,height
      parameter(pi=3.1415926535d0,deg=0.0174532920,rad=57.29577951d0)
      parameter(eps_earth=0.003352329869d0,r_earth=6378.140d3)
      
      call lst(t,lambda*rad,s0,ss,h)
      h=nula2pi(h/12*pi)

      e=sqrt(2*eps_earth-eps_earth**2)
      ss=(1-e**2)/sqrt(1-e**2*sin(fi)**2)
      cc=1/sqrt(1-e**2*sin(fi)**2)
      phi=atan((ss+height/r_earth)/(cc+height/r_earth)*tan(fi))
      rho0=(cc+height/r_earth)*cos(phi)/cos(fi)

      deltaxy=-sin(8.8/3600.*deg)*rho0*cos(fi)
      deltaz=-sin(8.8/3600.*deg)*rho0*sin(fi)

      dx(1)=deltaxy*cos(h)
      dx(2)=deltaxy*sin(h)
      dx(3)=deltaz
      return
      end

c**********************************************************************
      subroutine pv(elmts,x,v)
c**********************************************************************
c
c  calculation of position and velocity from elements
c  elmts(6)   elements of orbit
c  x(3),v(3)  position and velocity in rectangular coords
c
      double precision elmts(6),x(3),v(3)
      double precision a,e,i0,omega,omega0,m0,e0,
     &px,py,pz,qx,qy,qz,ax,ay,az,bx,by,bz,r,kg,ekepl1
      parameter(kg=0.017202098902128322400d0)

      a=elmts(1)
      e=elmts(2)
      i0=elmts(3)
      omega0=elmts(4)
      omega=elmts(5)
      m0=elmts(6)

      px=cos(omega0)*cos(omega)-sin(omega0)*sin(omega)*cos(i0)
      py=cos(omega0)*sin(omega)+sin(omega0)*cos(omega)*cos(i0)
      pz=sin(omega0)*sin(i0)
      qx=-sin(omega0)*cos(omega)-cos(omega0)*sin(omega)*cos(i0)
      qy=-sin(omega0)*sin(omega)+cos(omega0)*cos(omega)*cos(i0)
      qz=cos(omega0)*sin(i0)

      ax=a*px
      ay=a*py
      az=a*pz
      bx=a*sqrt(1-e**2)*qx
      by=a*sqrt(1-e**2)*qy
      bz=a*sqrt(1-e**2)*qz

      e0=ekepl1(m0,e)

      x(1)=ax*(cos(e0)-e)+bx*sin(e0)
      x(2)=ay*(cos(e0)-e)+by*sin(e0)
      x(3)=az*(cos(e0)-e)+bz*sin(e0)
      r=sqrt(x(1)**2+x(2)**2+x(3)**2)
      v(1)=kg/(r*sqrt(a))*(-ax*sin(e0)+bx*cos(e0))
      v(2)=kg/(r*sqrt(a))*(-ay*sin(e0)+by*cos(e0))
      v(3)=kg/(r*sqrt(a))*(-az*sin(e0)+bz*cos(e0))

      return
      end

c**********************************************************************
      subroutine e1(x,v,elmts,gm1)
c**********************************************************************
c
c  determination of elemnts from the position and velocity
c  x(3)      rectagular coordinates of the body [au]
c  v(3)      velocity [au/day] 
c  elmts(6)  elements of orbit
c
      double precision x(3),v(3),elmts(6)
      double precision v2,r,a,xx,yy,zz,p,e,i0,e0,n0,m0,v0,
     &omega,omega0,sinv,cosv,nula2pi,gms,kg
      real*8 gm1
      parameter(gms=2.9591220663860403d-4,kg=0.017202098902128322400d0)

      gms2=gms+gm1
      v2=v(1)**2+v(2)**2+v(3)**2
      r=sqrt(x(1)**2+x(2)**2+x(3)**2)
      a=1/(2/r-v2/gms2)
      xx=x(2)*v(3)-x(3)*v(2)
      yy=x(1)*v(3)-x(3)*v(1)
      zz=x(1)*v(2)-x(2)*v(1)
      p=(xx**2+yy**2+zz**2)/gms2
      e=sqrt(1-p/a)
      i0=atan(sqrt(xx**2+yy**2)/zz)
      omega=nula2pi(atan2(xx,yy))
      sinv=(x(1)*v(1)+x(2)*v(2)+x(3)*v(3))/r*sqrt(p/gms2)
      cosv=(p/r-1)
      v0=nula2pi(atan2(sinv,cosv))
      e0=2*atan(sqrt((1-e)/(1+e))*tan(v0/2))
      n0=kg*a**(-1.5)
      m0=nula2pi(e0-e*sin(e0))
      omega0=nula2pi(atan2(x(3)/sin(i0),x(1)*cos(omega)+
     &x(2)*sin(omega))-v0)

      elmts(1)=a
      elmts(2)=e
      elmts(3)=i0
      elmts(4)=omega0
      elmts(5)=omega
      elmts(6)=m0
      return
      end

c**********************************************************************
      integer function length(s)
c**********************************************************************
c
c  length of string s without spaces
c
      character s*(*)
      do 10 i=len(s),1,-1
        if (s(i:i).ne.' ') goto 20
10    continue
20    length=i
      return
      end

c**********************************************************************
      subroutine read_mpc(fname,name,t,alpha,delta,max)
c**********************************************************************
c
c  read mpc file
c  fname            filename
c  name             substring to find on line
c  t(1:max)         julian date
c  alpha(),delta()  equatoreal coordinates
c  max              maximal index, return the number of lines
c
      character fname*30,name*7
      integer max
      double precision t(max),alpha(max),delta(max)
      integer i,length
      double precision x1,x2,x3,jd,pi,deg
      character fline*80
      parameter(pi=3.1415926535d0,deg=0.017453292d0)
      
      open(1,file=fname)
      i=0
10    if (i.ge.max) goto 100 
        read(1,20,end=100,err=110) fline
20      format(a)
        if ((index(fline(1:12),name(1:length(name))).ne.0).and.
     &  (fline(15:15).eq.'C')) then
          i=i+1
          read(fline(16:20),'(d4.0)') x1
          read(fline(21:23),'(d2.0)') x2
          read(fline(24:32),'(d8.5)') x3
          t(i)=jd(x1,x2,x3)
          read(fline(33:35),'(d2.0)') x1
          read(fline(36:38),'(d2.0)') x2
          read(fline(39:44),'(d5.2)') x3
          alpha(i)=(x1+x2/60+x3/3600)/12*pi
          read(fline(46:48),'(d2.0)') x1
          read(fline(49:51),'(d2.0)') x2
          read(fline(52:56),'(d5.2)') x3
          delta(i)=(x1+x2/60+x3/3600)*deg
          if (fline(45:45).eq.'-') delta(i)=-delta(i)
        endif
      goto 10
100   continue
      close(1)
      max=i
      return
110   write(*,115)
115   format(1x,'Error reading input.')     
      return
      end

c**********************************************************************
      subroutine decode_elmts(fline,elmts)
c**********************************************************************
c
c  decode elements from a given string
c  fline*131     string to decode
c  elmts(6)      julian date
c
      
      character fline*131
      double precision elmts(6)
      double precision x1,x2,x3,jd,deg,nula2pi,k_gauss
      parameter(deg=0.017453292d0,k_gauss=0.01720209895d0)
      
      read(fline(36:45),'(d10.7)') elmts(1)
      read(fline(46:55),'(d9.7)') elmts(2)
      read(fline(56:64),'(d8.5)') elmts(3)
      read(fline(65:74),'(d9.5)') elmts(5)
      read(fline(75:84),'(d9.5)') elmts(4)
      elmts(3)=elmts(3)*deg
      elmts(4)=elmts(4)*deg
      elmts(5)=elmts(5)*deg
      if (fline(85:85).eq.'M') then
        read(fline(86:95),'(d9.5)') x1
        read(fline(96:105),'(d9.1)') x2
        elmts(6)=nula2pi(k_gauss*elmts(1)**(-1.5)*
     &    (2451545-x2)+x1*deg)
      else
        read(fline(85:86),'(d2.0)') x1
        read(fline(87:88),'(d2.0)') x2
        read(fline(89:96),'(d8.5)') x3
        elmts(6)=nula2pi(k_gauss*elmts(1)**(-1.5)*
     &    (2451545-jd(1900+x1,x2,x3)))
      endif
      return
      end

c**********************************************************************
      subroutine read_elmts(fname,name,elmts,max)
c**********************************************************************
c
c  read elements from dat file, this subroutine calls decode_elmts
c  after finding the proper line
c  fname         filename
c  name          substring to find on line
c  elmts(1:6)    julian date
c  max           maximal number of elements
c
      character fname*30,name*9
      integer max
      double precision elmts(6)
      integer i,length
      character fline*131
      
      open(1,file=fname)
      i=0
10    read(1,20,end=100,err=110) fline
20    format(a)
      if ((index(fline,'BEGIN')).ne.0) goto 10
      
30    if (i.ge.max) goto 100
        read(1,20,end=100,err=110) fline
        if (index(fline(1:29),name(1:length(name)+1)).ne.0) then
          i=i+1
          call decode_elmts(fline,elmts)
        endif
      goto 30
100   continue
      close(1)
      return
110   write(*,115)
115   format(1x,'Error reading input.')     
      return
      end

c**********************************************************************
      subroutine eclequ(elmts)
c**********************************************************************
c
c  tranformation of elements of orbit from ecliptical coordinate frame
c  to equatoreal
c
      double precision elmts(6)
      double precision i,omega,omega0,nula2pi

      parameter(pi=3.1415926535d0,epsilon=0.409099609d0)

      i=elmts(3)
      omega0=elmts(4)
      omega=elmts(5)

      elmts(5)=atan2(sin(i)*sin(omega),
     &cos(i)*sin(epsilon)+sin(i)*cos(epsilon)*cos(omega))

      elmts(3)=atan2(sin(i)*sin(omega)/sin(elmts(5)),
     &cos(i)*cos(epsilon)-sin(i)*sin(epsilon)*cos(omega))

      elmts(4)=nula2pi(elmts(4)+atan2(sin(epsilon)*sin(omega),
     &sin(i)*cos(epsilon)+cos(i)*sin(epsilon)*cos(omega)))

      return
      end

c**********************************************************************
      subroutine equecl(elmts)
c**********************************************************************
c
c  tranformation of elements of orbit from equatoreal coordinate frame
c  to eliptical
c
      double precision elmts(6)
      double precision i,omega,omega0,nula2pi

      parameter(pi=3.1415926535d0,epsilon=0.409099609d0)

      i=elmts(3)
      omega0=elmts(4)
      omega=elmts(5)

      elmts(5)=atan2(sin(i)*sin(omega),
     &-cos(i)*sin(epsilon)+sin(i)*cos(epsilon)*cos(omega))

      elmts(3)=atan2(sin(i)*sin(omega)/sin(elmts(5)),
     &cos(i)*cos(epsilon)+sin(i)*sin(epsilon)*cos(omega))

      elmts(4)=nula2pi(elmts(4)-atan2(sin(epsilon)*sin(omega),
     &sin(i)*cos(epsilon)-cos(i)*sin(epsilon)*cos(omega)))

      return
      end

c**********************************************************************
      subroutine de_rk4(x,y0,y1,h,x0,de_y2)
c**********************************************************************
c
c  integration of the 2nd order differential equotions, 
c  method runge-kutta, 4th order
c  2nd derivative is calculated by calling subroutine de_y2(x,y0,y1,y2)
c  x      independent variable
c  y0(n)  positions
c  y1(n)  velocities
c  h      integrating step
c  x0     time to stop integration
c
c  parameters in common block /numeq/:
c  n  number of equotions
c
      integer i,n
      common/numeq/n
c
c  max number of equotions
c
      parameter(m=36)
      double precision x,y0(n),y1(n),h,x0
      double precision k1(m),k2(m),k3(m),k4(m)
      double precision xp,y0p(m),y1p(m)
      external de_y2

1     continue

      call de_y2(x,y0,y1,k1)
      xp=x+h*0.5
      do i=1,n
        y1p(i)=y1(i)+k1(i)*h*0.5
        y0p(i)=y0(i)+y1(i)*h*0.5+k1(i)*h**2*0.125
      enddo
      call de_y2(xp,y0p,y1p,k2)
      do i=1,n
        y1p(i)=y1(i)+k2(i)*h*0.5
        y0p(i)=y0(i)+y1(i)*h*0.5+k2(i)*h**2*0.125
      enddo
      call de_y2(xp,y0p,y1p,k3)
      xp=x+h
      do i=1,n
        y1p(i)=y1(i)+k3(i)*h
        y0p(i)=y0(i)+y1(i)*h+k3(i)*h**2*0.5
      enddo
      call de_y2(xp,y0p,y1p,k4)
      do i=1,n
        y0(i)=y0(i)+y1(i)*h+(k1(i)+k2(i)+k3(i))*h**2/6.
        y1(i)=y1(i)+(k1(i)+2*(k2(i)+k3(i))+k4(i))*h/6.
      enddo
      x=x+h

      if (x.lt.x0) goto 1
      return
      end

c**********************************************************************
      subroutine de_pcm(t,x,v,h,t0,de_a)
c**********************************************************************
c
c  integration of equotions of motion, predictor-corrector method
c  of the 2nd order with modification, acceleration is calculated
c  by calling subroutine de_a(t,x,v,a)
c  t       time
c  x(n)  positions of n bodies
c  v(n)  velocities
c  h       time step
c  t0      time to stop integration
c
c  parameters in common block /numeq/:
c  n  number of equotions
c
      integer i,n,m
      common/numeq/n
c
c  max number of bodies
c
      parameter(m=39)
      double precision t,x(n),v(n),h,t0
      double precision t1,x1(m),v1(m),a1(m),t2,a2(m)
      double precision px1(m),pv1(m),kx1(m),kv1(m)
      double precision px2(m),pv2(m),kx2(m),kv2(m)
      double precision mx2(m),mv2(m)
      external de_a
c
c  one step with rk4 method
c
      t1=t
      do i=1,n
        x1(i)=x(i)
        v1(i)=v(i)
      enddo
      call de_rk4(t1,x1,v1,h,t1+h,de_a)
      t2=t1+h
c
c  another one without modificator
c
      call de_a(t1,x1,v1,a1)
      do i=1,n
        pv1(i)=v(i)+2*h*a1(i)
        px1(i)=x(i)+2*h*v1(i)
      enddo

      call de_a(t2,px1,pv1,a2)
      do i=1,n
        kv1(i)=v1(i)+h/2.*(a1(i)+a2(i))
        kx1(i)=x1(i)+h/2.*(v1(i)+pv1(i))

        x(i)=x1(i)
        v(i)=v1(i)

        v1(i)=kv1(i)+1./5.*(pv1(i)-kv1(i))
        x1(i)=kx1(i)+1./5.*(px1(i)-kx1(i))
      enddo
c
c  continue pcm method
c
1     continue

      call de_a(t1,x1,v1,a1)
      do i=1,n
        pv2(i)=v(i)+2*h*a1(i)
        px2(i)=x(i)+2*h*v1(i)

        mv2(i)=pv2(i)-4./5.*(pv1(i)-kv1(i))
        mx2(i)=px2(i)-4./5.*(px1(i)-kx1(i))
      enddo

      call de_a(t2,mx2,mv2,a2)
      do i=1,n
        kv2(i)=v1(i)+h/2.*(a1(i)+a2(i))
        kx2(i)=x1(i)+h/2.*(v1(i)+mv2(i))
c
c  it is not necessery to declare arrays x2(m),v2(m)
c
        x(i)=x1(i)
        v(i)=v1(i)

        v1(i)=kv2(i)+1./5.*(mv2(i)-kv2(i))
        x1(i)=kx2(i)+1./5.*(mx2(i)-kx2(i))
c
c  save predictors and correctors for modification in the next step
c
        px1(i)=px2(i)
        pv1(i)=pv2(i)
        kx1(i)=kx2(i)
        kv1(i)=kv2(i)
      enddo

      t=t1
      t1=t2
      t2=t2+h

      if (t2.lt.t0) goto 1

      t=t2
      do i=1,n
        x(i)=x1(i)
        v(i)=v1(i)
      enddo

      return
      end

c**********************************************************************
      subroutine de_ham(t,x,v,h,t0,de_a)
c**********************************************************************
c
c  integration of equotions of motion, hamming method of the 4th order
c  with modification, acceleration is calculated by calling subroutine
c  de_a(t,x,v,a)
c  t       time
c  x(n)  positions of n bodies
c  v(n)  velocities
c  h       time step
c  t0      time to stop integration
c
c  parameters in common block /numeq/:
c  n  number of equotions
c
      integer i,n,m
      common/numeq/n
c
c  max number of bodies
c
      parameter(m=39)
      double precision t,x(n),v(n),h,t0
      double precision t1,x1(m),v1(m),a1(m)
      double precision t2,x2(m),v2(m),a2(m)
      double precision t3,x3(m),v3(m),a3(m),t4,a4(m)
      double precision px3(m),pv3(m),kx3(m),kv3(m)
      double precision px4(m),pv4(m),kx4(m),kv4(m)
      double precision mx4(m),mv4(m)
      external de_a
c
c  three steps with rk4 method
c
      t1=t
      do i=1,n
        x1(i)=x(i)
        v1(i)=v(i)
      enddo
      call de_rk4(t1,x1,v1,h,t1+h,de_a)
      t2=t1
      do i=1,n
        x2(i)=x1(i)
        v2(i)=v1(i)
      enddo
      call de_rk4(t2,x2,v2,h,t2+h,de_a)
      t3=t2
      do i=1,n
        x3(i)=x2(i)
        v3(i)=v2(i)
      enddo
      call de_rk4(t3,x3,v3,h,t3+h,de_a)
      t4=t3+h
c
c  another one without modificator
c
      call de_a(t1,x1,v1,a1)
      call de_a(t2,x2,v2,a2)
      call de_a(t3,x3,v3,a3)
      do i=1,n
        pv3(i)=v(i)+4./3.*h*(2*(a1(i)+a3(i))-a2(i))
        px3(i)=x(i)+4./3.*h*(2*(v1(i)+v3(i))-v2(i))
      enddo

      call de_a(t4,px3,pv3,a4)
      do i=1,n
        kv3(i)=1./8.*(9*v3(i)-v1(i))
     &    +3./8.*h*(a4(i)+2*a3(i)-a2(i))
        kx3(i)=1./8.*(9*x3(i)-x1(i))
     &    +3./8.*h*(pv3(i)+2*v3(i)-v2(i))

        x(i)=x1(i)
        v(i)=v1(i)
        x1(i)=x2(i)
        v1(i)=v2(i)
        a1(i)=a2(i)
        x2(i)=x3(i)
        v2(i)=v3(i)
        a2(i)=a3(i)

        v3(i)=kv3(i)+9./121.*(pv3(i)-kv3(i))
        x3(i)=kx3(i)+9./121.*(px3(i)-kx3(i))
      enddo
c
c  continue hamming method
c
1     continue

      call de_a(t3,x3,v3,a3)
      do i=1,n
        pv4(i)=v(i)+4./3.*h*(2*(a1(i)+a3(i))-a2(i))
        px4(i)=x(i)+4./3.*h*(2*(v1(i)+v3(i))-v2(i))

        mv4(i)=pv4(i)-112./121.*(pv3(i)-kv3(i))
        mx4(i)=px4(i)-112./121.*(px3(i)-kx3(i))
      enddo

      call de_a(t4,mx4,mv4,a4)
      do i=1,n
        kv4(i)=1./8.*(9*v3(i)-v1(i))
     &    +3./8.*h*(a4(i)+2*a3(i)-a2(i))
        kx4(i)=1./8.*(9*x3(i)-x1(i))
     &    +3./8.*h*(mv4(i)+2*v3(i)-v2(i))

        x(i)=x1(i)
        v(i)=v1(i)
        x1(i)=x2(i)
        v1(i)=v2(i)
        a1(i)=a2(i)
        x2(i)=x3(i)
        v2(i)=v3(i)
        a2(i)=a3(i)

        v3(i)=kv4(i)+9./121.*(mv4(i)-kv4(i))
        x3(i)=kx4(i)+9./121.*(mx4(i)-kx4(i))
c
c  save predictors and correctors for modification in the next step
c
        px3(i)=px4(i)
        pv3(i)=pv4(i)
        kx3(i)=kx4(i)
        kv3(i)=kv4(i)
      enddo

      t=t1
      t1=t2
      t2=t3
      t3=t4
      t4=t4+h

      if (t4.lt.t0) goto 1

      t=t4
      do i=1,n
        x(i)=x3(i)
        v(i)=v3(i)
      enddo

      return
      end

c**********************************************************************
      subroutine de_ada(t,x,v,h,t0,de_a)
c**********************************************************************
c
c  integration of equotions of motion, adams-moulton method 
c  of the 4th order with modification, acceleration is calculated 
c  by calling subroutine de_a(t,x,v,a)
c  t       time
c  x(n)  positions of n bodies
c  v(n)  velocities
c  h       time step
c  t0      time to stop integration
c
c  parameters in common block /numeq/:
c  n  number of equotions
c
      integer i,n,m
      common/numeq/n
c
c  max number of bodies
c
      parameter(m=39)
      double precision t,x(n),v(n),h,t0
      double precision t1,x1(m),v1(m),a1(m),a(m)
      double precision t2,x2(m),v2(m),a2(m)
      double precision t3,x3(m),v3(m),a3(m),t4,a4(m)
      double precision px3(m),pv3(m),kx3(m),kv3(m)
      double precision px4(m),pv4(m),kx4(m),kv4(m)
      double precision mx4(m),mv4(m)
      external de_a
c
c  three steps with rk4 method
c
      t1=t
      do i=1,n
        x1(i)=x(i)
        v1(i)=v(i)
      enddo
      call de_rk4(t1,x1,v1,h,t1+h,de_a)
      t2=t1
      do i=1,n
        x2(i)=x1(i)
        v2(i)=v1(i)
      enddo
      call de_rk4(t2,x2,v2,h,t2+h,de_a)
      t3=t2
      do i=1,n
        x3(i)=x2(i)
        v3(i)=v2(i)
      enddo
      call de_rk4(t3,x3,v3,h,t3+h,de_a)
      t4=t3+h
c
c  another one without modificator
c
      call de_a(t,x,v,a)
      call de_a(t1,x1,v1,a1)
      call de_a(t2,x2,v2,a2)
      call de_a(t3,x3,v3,a3)
      do i=1,n
        pv3(i)=v3(i)
     &    +h/24.*(55.*a3(i)-59.*a2(i)+37.*a1(i)-9.*a(i))
        px3(i)=x3(i)
     &    +h/24.*(55.*v3(i)-59.*v2(i)+37.*v1(i)-9.*v(i))
      enddo

      call de_a(t4,px3,pv3,a4)
      do i=1,n
        kv3(i)=v3(i)
     &    +h/24.*(9.*a4(i)+19.*a3(i)-5.*a2(i)+a1(i))
        kx3(i)=x3(i)
     &    +h/24.*(9.*pv3(i)+19.*v3(i)-5.*v2(i)+v1(i))
        x(i)=x1(i)
        v(i)=v1(i)
        a(i)=a1(i)
        x1(i)=x2(i)
        v1(i)=v2(i)
        a1(i)=a2(i)
        x2(i)=x3(i)
        v2(i)=v3(i)
        a2(i)=a3(i)

        v3(i)=kv3(i)+19./270.*(pv3(i)-kv3(i))
        x3(i)=kx3(i)+19./270.*(px3(i)-kx3(i))
      enddo
c
c  continue adams-moulton method
c
1     continue

      call de_a(t3,x3,v3,a3)
      do i=1,n
        pv4(i)=v3(i)
     &    +h/24.*(55.*a3(i)-59.*a2(i)+37.*a1(i)-9.*a(i))
        px4(i)=x3(i)
     &    +h/24.*(55.*v3(i)-59.*v2(i)+37.*v1(i)-9.*v(i))

        mv4(i)=pv4(i)-251./270.*(pv3(i)-kv3(i))
        mx4(i)=px4(i)-251./270.*(px3(i)-kx3(i))
      enddo

      call de_a(t4,mx4,mv4,a4)
      do i=1,n
        kv4(i)=v3(i)
     &    +h/24.*(9.*a4(i)+19.*a3(i)-5.*a2(i)+a1(i))
        kx4(i)=x3(i)
     &    +h/24.*(9.*mv4(i)+19.*v3(i)-5.*v2(i)+v1(i))

        x(i)=x1(i)
        v(i)=v1(i)
        a(i)=a1(i)
        x1(i)=x2(i)
        v1(i)=v2(i)
        a1(i)=a2(i)
        x2(i)=x3(i)
        v2(i)=v3(i)
        a2(i)=a3(i)

        v3(i)=kv4(i)+19./270.*(mv4(i)-kv4(i))
        x3(i)=kx4(i)+19./270.*(mx4(i)-kx4(i))
c
c  save predictors and correctors for modification in the next step
c
        px3(i)=px4(i)
        pv3(i)=pv4(i)
        kx3(i)=kx4(i)
        kv3(i)=kv4(i)
      enddo

      t=t1
      t1=t2
      t2=t3
      t3=t4
      t4=t4+h
         
      if (t4.lt.t0) goto 1

      t=t4
      do i=1,n
        x(i)=x3(i)
        v(i)=v3(i)
      enddo

      return
      end

c**********************************************************************
      subroutine de_adc(t,x,v,h,t0,de_a)
c**********************************************************************
c
c  integration of equotions of motion, adams-moulton method 
c  of the 4th order with modification and time step control, 
c  acceleration is calculated by calling subroutine de_a(t,x,v,a)
c  t     time
c  x(n)  positions of n bodies
c  v(n)  velocities
c  h     time step
c  t0     time to stop integration
c
c  parameters in common block /numeq/:
c  n  number of equotions
c
c  parameters in common block /local_err/:
c  max_err, min_err  maximal and minimal local relative error
c
      integer i,j,n,m
      common/numeq/n
      common/local_err/max_err,min_err
      double precision min_err,max_err
c
c  max number of bodies
c
      parameter(m=39)
      double precision t,x(n),v(n),h,t0
      double precision t1,x1(m),v1(m),a1(m),a(m)
      double precision t2,x2(m),v2(m),a2(m)
      double precision t3,x3(m),v3(m),a3(m),t4,a4(m)
      double precision px3(m),pv3(m),kx3(m),kv3(m)
      double precision px4(m),pv4(m),kx4(m),kv4(m)
      double precision mx4(m),mv4(m),temp
      logical incr,decr
      external de_a
c
c  three steps with rk4 method
c
100   continue
      incr=.true.
      decr=.false.
      t1=t
      do i=1,n
        x1(i)=x(i)
        v1(i)=v(i)
      enddo
      call de_rk4(t1,x1,v1,h,t1+h,de_a)
      t2=t1
      do i=1,n
        x2(i)=x1(i)
        v2(i)=v1(i)
      enddo
      call de_rk4(t2,x2,v2,h,t2+h,de_a)
      t3=t2
      do i=1,n
        x3(i)=x2(i)
        v3(i)=v2(i)
      enddo
      call de_rk4(t3,x3,v3,h,t3+h,de_a)
      t4=t3+h
c
c  another one without modificator
c
      call de_a(t,x,v,a)
      call de_a(t1,x1,v1,a1)
      call de_a(t2,x2,v2,a2)
      call de_a(t3,x3,v3,a3)
      do i=1,n
        pv3(i)=v3(i)
     &    +h/24.*(55.*a3(i)-59.*a2(i)+37.*a1(i)-9.*a(i))
        px3(i)=x3(i)
     &    +h/24.*(55.*v3(i)-59.*v2(i)+37.*v1(i)-9.*v(i))
      enddo

      call de_a(t4,px3,pv3,a4)
      do i=1,n
        kv3(i)=v3(i)
     &    +h/24.*(9.*a4(i)+19.*a3(i)-5.*a2(i)+a1(i))
        kx3(i)=x3(i)
     &    +h/24.*(9.*pv3(i)+19.*v3(i)-5.*v2(i)+v1(i))
        x(i)=x1(i)
        v(i)=v1(i)
        a(i)=a1(i)
        x1(i)=x2(i)
        v1(i)=v2(i)
        a1(i)=a2(i)
        x2(i)=x3(i)
        v2(i)=v3(i)
        a2(i)=a3(i)

        v3(i)=kv3(i)+19./270.*(pv3(i)-kv3(i))
        x3(i)=kx3(i)+19./270.*(px3(i)-kx3(i))
      enddo
c
c  continue adams-moulton method
c
1     continue

      call de_a(t3,x3,v3,a3)
      do i=1,n
        pv4(i)=v3(i)
     &    +h/24.*(55.*a3(i)-59.*a2(i)+37.*a1(i)-9.*a(i))
        px4(i)=x3(i)
     &    +h/24.*(55.*v3(i)-59.*v2(i)+37.*v1(i)-9.*v(i))

        mv4(i)=pv4(i)-251./270.*(pv3(i)-kv3(i))
        mx4(i)=px4(i)-251./270.*(px3(i)-kx3(i))
      enddo

      call de_a(t4,mx4,mv4,a4)
      do i=1,n
        kv4(i)=v3(i)
     &    +h/24.*(9.*a4(i)+19.*a3(i)-5.*a2(i)+a1(i))
        kx4(i)=x3(i)
     &    +h/24.*(9.*mv4(i)+19.*v3(i)-5.*v2(i)+v1(i))

        x(i)=x1(i)
        v(i)=v1(i)
        a(i)=a1(i)
        x1(i)=x2(i)
        v1(i)=v2(i)
        a1(i)=a2(i)
        x2(i)=x3(i)
        v2(i)=v3(i)
        a2(i)=a3(i)

        v3(i)=kv4(i)+19./270.*(mv4(i)-kv4(i))
        x3(i)=kx4(i)+19./270.*(mx4(i)-kx4(i))
c
c  save predictors and correctors for modification in the next step
c
        px3(i)=px4(i)
        pv3(i)=pv4(i)
        kx3(i)=kx4(i)
        kv3(i)=kv4(i)
        temp=abs((kx4(i)-px4(i))/kx4(i))
        if (temp.gt.max_err) decr=.true.
        if (temp.gt.min_err) incr=.false.
      enddo
      if (decr.or.incr) then
        t=t2
        do j=1,n
          x(j)=x2(j)
          v(j)=v2(j)
        enddo
        if (decr) then
          h=h/2
        else
          h=2*h
        endif
        goto 100
      endif
      incr=.true.
      t=t1
      t1=t2
      t2=t3
      t3=t4
      t4=t4+h
         
      if (t4.lt.t0) goto 1

      t=t4
      do i=1,n
        x(i)=x3(i)
        v(i)=v3(i)
      enddo

      return
      end

