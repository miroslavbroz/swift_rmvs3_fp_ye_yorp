
c**********************************************************************
      subroutine pv1(elmts,x,v,gm1)
c**********************************************************************
c
c  calculation of position and velocity from elements
c  elmts(6)   elements of orbit
c  x(3),v(3)  position and velocity in rectangular coords
c  gm1        mass of the planet (or zero)
c
      implicit none
      double precision elmts(6),x(3),v(3),gm1
      double precision a,e,i0,omega,omega0,m0,e0,
     &px,py,pz,qx,qy,qz,ax,ay,az,bx,by,bz,r,gms,kg,gms2,kg2,ekepl1
c      parameter(gms=2.9591220663860403d-4,kg=0.017202098902128322400d0)
      parameter(gms=0.2959122082855911d-03)     ! a new value from JPLEPH

      gms2=gms+gm1
      kg2=sqrt(gms2)

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
      v(1)=kg2/(r*sqrt(a))*(-ax*sin(e0)+bx*cos(e0))
      v(2)=kg2/(r*sqrt(a))*(-ay*sin(e0)+by*cos(e0))
      v(3)=kg2/(r*sqrt(a))*(-az*sin(e0)+bz*cos(e0))

      return
      end


