#!/usr/bin/gnuplot

load "config.plt"

set xl "t [My]"
set yl "a [au]"
set cbl "TP"

#set yr [3.3:3.6]
set zeroaxis
set key left width -1

set palette rgbformulae 33,13,10

p "<awk '($1>0)' follow.out" u 2:3:1 lc 'gray' t "osculating",\
  "<awk '($1>0)' follow.filter.out" u 2:3:1 lc palette z t "mean",\
  "<awk '($1>0)' follow.proper.out" u 2:3:1 lc palette z t "proper",\
  a0+dadt*x w l lc 'cyan' t sprintf("%.4e au/My", dadt),\
  a0 w l lt 0

pa -1

set term png small
set out "at.png"
rep


