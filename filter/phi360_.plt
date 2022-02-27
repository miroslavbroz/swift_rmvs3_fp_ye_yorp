#!/usr/bin/gnuplot

set xl "t [yr]"
set yl "a-a' [au]"
set y2l "lambda-lambda'-sigma [deg]"

set zeroaxis
set ytics nomirror
set y2tics

p \
  "phi360.dat" u ($3/365.25):($7/0.2783) w lp,\
  "phi360.dat" u ($3/365.25):8 ax x1y2 w lp,\
  "phi360.dat" u ($3/365.25):($4/360.) ax x1y2 w lp,\
  "phi360.dat" u ($3/365.25):(atan2($7/0.2783,$8/180.*pi)/pi*180./360.) ax x1y2 w lp,\
  0.0 ax x1y2 w l lt 0

pa -1


