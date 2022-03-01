#!/usr/bin/gnuplot

G = 6.67403e-11  # kg^-1 m^3 s^-2
gram = 1.e-3  # kg
cm = 1.e-2  # m
km = 1.e3  # m
hour = 3600.  # s
day = 24.*hour
rho = 2500.  # kg m^-3
D = 10.  # m

omega(P) = 2.*pi/(P*hour)
omega_g(rho) = sqrt(4./3.*pi*G*rho)

print "rho = ", rho, " kg m^-3"
print "omega_g = ", omega_g(rho), " rad/s = ", omega_g(rho)/(2.*pi/day), " rev/day"
print "P_crit = ", 2.*pi/omega_g(rho)/hour, " h"

# Holsapple (2004), Eqs. (5.7), (5.9)
s = 0.25
alpha = 0.7
beta = 0.7
kappa = 2.25e7 * gram*cm**(-1./2.)
k(D) = kappa*(D/cm/2.)**(-1./2.)
C = (alpha*beta)**(1./3.) * sqrt(  5.*(3.*s*(1.+beta**2) - sqrt(3.*(1.-beta**2+beta**2))) / (3.*s**2*(1.+beta**2)**2-1.+beta**2-beta**4)  )

omega_s(D,rho) = C*sqrt(k(D)/rho)*(1./(D/2.))

print "D = ", D, " m"
print "omega_s = ", omega_s(D,rho), " rad/s = ", omega_s(D,rho)/(2.*pi/day), " rev/day"

########################################################################

set xl "t [My]"
set yl "omega [rad/s]"
set cbl "TP"

set zeroaxis
set key left
set palette rgbformulae 33,13,10

p "yorp.out" u 2:6:1 w p lc palette z,\
  omega_g(rho) w l lt 0,\
  0.0 w l lt 0 not

pa -1

set term png small
set out "omegat.png"
rep

q


  omega_s(D, rho) w l lt 0,\
  -omega_s(D, rho) w l lt 0 not,\

