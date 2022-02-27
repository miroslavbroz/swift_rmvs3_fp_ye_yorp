#!/usr/bin/gnuplot

set xl "obliquity gamma [deg]"
set yl "f-function [rad/s]"

p "<./intep_test 1.y 0 180 1" u 1:2,\
  "1.y" u 1:2
pa -1

set term png small
set out "f_functions.png"
rep


