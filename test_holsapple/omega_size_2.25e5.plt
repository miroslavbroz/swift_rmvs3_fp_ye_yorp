#!/usr/bin/gnuplot

G = 6.67403e-11  # kg^-1 m^3 s^-2
gram = 1.e-3  # kg
cm = 1.e-2  # m
km = 1.e3  # m
hour = 3600.  # s
day = 24.*hour
rho = 2000.  # kg m^-3
D = 1.  # m

omega(P) = 2.*pi/(P*hour)
omega_g(rho) = sqrt(4./3.*pi*G*rho)

print "rho = ", rho, " kg m^-3"
print "omega_g = ", omega_g(rho), " rad/s = ", omega_g(rho)/(2.*pi/day), " rev/day"
print "P_crit = ", 2.*pi/omega_g(rho)/hour, " h"

# Holsapple (2004), Eqs. (5.7), (5.9)
s = 0.25
alpha = 0.7
beta = 0.7
kappa = 2.25e5 * gram*cm**(-1./2.)
k(D) = kappa*(D/cm/2.)**(-1./2.)
C = (alpha*beta)**(1./3.) * sqrt(  5.*(3.*s*(1.+beta**2) - sqrt(3.*(1.-beta**2+beta**2))) / (3.*s**2*(1.+beta**2)**2-1.+beta**2-beta**4)  )

omega_s(D,rho) = C*sqrt(k(D)/rho)*(1./(D/2.))

print "D = ", D, " m"
print "omega_s = ", omega_s(D,rho), " rad/s"

########################################################################

set xl "D [km]"
set yl "omega [rad/s]"
set cbl "A [mag]"

set xr [0.001:1e4]
set yr [1e-6:1]
set cbr [0:1.0]
set logscale x
set logscale y
set palette rgbformulae 33,13,10
 
set arrow from 10.,graph 0 rto 0,graph 1 nohead lt 0 front
set arrow from 1.,graph 0 rto 0,graph 1 nohead lt 0 front

p "lc_summary.out" u 1:(omega($2)):3 lc palette z,\
  omega_g(1000.) w l lt 0 t "rho = 1000 kg m^{-3}",\
  omega_g(2000.) w l lt 0 t "2000",\
  omega_g(3000.) w l lt 0 t "3000",\
  omega_s(x*km,rho) lc 'black' lw 1,\
  omega_g(rho)+omega_s(x*km,rho) lc 'red' lw 2 dt 2,\
  "omega_s.dat" u 1:2 w l lt 0 t "Holsapple (2005)"

pa -1

set term png small
set out "omega_size_2.25e5.png"
rep


