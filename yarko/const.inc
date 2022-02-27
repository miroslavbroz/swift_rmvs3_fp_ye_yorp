c*************************************************************************
c CONST.INC
c*************************************************************************
c Fundamental constatns.
c
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Feb 12th 2002

      real*8 sigma,cv,S0,k_gauss,AU,day
c  Boltzmann constant [J.K^(-4).m^(-2).s^(-1)]
      parameter(sigma=5.66962d-8)
c  speed of light [m.s^(-1)]
      parameter(cv=299792458.0d0)
c  solar constant [W.m^(-2)]
      parameter(S0=1371)
c  Gauss constant
      parameter(k_gauss=0.01720209895d0)
c  astronomical unit [m]
      parameter(AU=149597870.d3)
c  day [sec]
      parameter(day=86400.d0)
c  gravitational constant [SI units]
      real*8 capG
      parameter(capG = 6.67428d-11)


c  older values, above corrected according to DV's code
c      parameter(sigma=5.6703d-8)
c      parameter(S0=1370)


