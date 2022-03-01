#!/usr/bin/gnuplot

set xl "t [My]"
set yl "e []"
set cbl "TP"

set zeroaxis
set key left width -1

set palette rgbformulae 33,13,10

p "<awk '($1>0)' follow.out" u 2:4:1 lc 'gray' t "osculating",\
  "<awk '($1>0)' follow.filter.out" u 2:4:1 lc palette z t "mean",\
  "<awk '($1>0)' follow.proper.out" u 2:4:1 lc palette z t "proper"

pa -1

set term png small
set out "et.png"
rep


