#!/usr/bin/gnuplot

set xl "obliquity gamma [deg]"
set yl "f-function [rad/s]"

p "<./intep_test 1.y 0 180 1" u 1:2 w lp,\
  "1.y" u 1:2 w lp
pa -1

set term png small
set out "f_functions1.png"
rep


