#!/usr/bin/gnuplot

set xl "t [yr]"
set yl "phi360 [deg]"

set zeroaxis

p \
  "phi360.dat" u ($3/365.25):4 w lp,\
  "phi360.dat" u ($3/365.25):5 w lp,\
  "phi360.dat" u ($3/365.25):6 w lp

pa -1


