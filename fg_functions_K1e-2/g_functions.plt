#!/usr/bin/gnuplot

set xl "obliquity gamma [deg]"
set yl "g-function [rad/s]"

set nokey
set zeroaxis

p "<./intep_test.sh   1.y 0 180 1  1 3" w l,\
  "<./intep_test.sh   2.y 0 180 1  1 3" w l,\
  "<./intep_test.sh   3.y 0 180 1  1 3" w l,\
  "<./intep_test.sh   4.y 0 180 1  1 3" w l,\
  "<./intep_test.sh   5.y 0 180 1  1 3" w l,\
  "<./intep_test.sh   6.y 0 180 1  1 3" w l,\
  "<./intep_test.sh   7.y 0 180 1  1 3" w l,\
  "<./intep_test.sh   8.y 0 180 1  1 3" w l,\
  "<./intep_test.sh   9.y 0 180 1  1 3" w l,\
  "<./intep_test.sh  10.y 0 180 1  1 3" w l
pa -1

set term png small
set out "g_functions.png"
rep


