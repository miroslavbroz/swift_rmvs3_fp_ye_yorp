#!/usr/bin/gnuplot

hour = 3600.  # s

P = 4.*hour
omega = 2.*pi/P

print "omega = ", omega, " rad s^-1"

GM = 2.9591221047429078E-04  # au^3 d^-2
a = 3.42  # au
e = 0.15  # 1
i = 1.e-3  # rad

vkepl = sqrt(GM/a)
rp = a*(1.-e)
vp = vkepl*sqrt((1.+e)/(1.-e))
x = rp*cos(i)
y = 0.0
z = rp*sin(i)
vx = 0.0
vy = vp*cos(i)
vz = vp*sin(i)

print "vkepl = ", vkepl, " au d^-1"
print "vp = ", vp, " au d^-1"
print "rp = ", rp, " au"
print ""
print x, y, z
print vx, vy, vz
print 0
print 0.0
print ""

