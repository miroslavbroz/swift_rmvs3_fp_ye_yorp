#!/usr/bin/gnuplot

set xl "t [My]"
set yl "sz []"
set cbl "TP"

tmp=1.2
set yr [-tmp:tmp]

set zeroaxis
set palette rgbformulae 33,13,10

p "yorp.out" u 2:5:1 w p lc palette z,\
  -1.0 w l lt 0 not,\
  1.0 w l lt 0 not

pa -1

set term png small
set out "sz.png"
rep

q


