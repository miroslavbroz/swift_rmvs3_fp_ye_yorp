#!/usr/bin/gnuplot

set xl "obliquity gamma [deg]"
set yl "f-function [rad/s]"

set nokey
set zeroaxis

p "<./intep_test  1.y 0 180 1" u 1:2 w l,\
  "<./intep_test  2.y 0 180 1" u 1:2 w l,\
  "<./intep_test  3.y 0 180 1" u 1:2 w l,\
  "<./intep_test  4.y 0 180 1" u 1:2 w l,\
  "<./intep_test  5.y 0 180 1" u 1:2 w l,\
  "<./intep_test  6.y 0 180 1" u 1:2 w l,\
  "<./intep_test  7.y 0 180 1" u 1:2 w l,\
  "<./intep_test  8.y 0 180 1" u 1:2 w l,\
  "<./intep_test  9.y 0 180 1" u 1:2 w l,\
  "<./intep_test 10.y 0 180 1" u 1:2 w l
pa -1

set term png small
set out "f_functions.png"
rep


