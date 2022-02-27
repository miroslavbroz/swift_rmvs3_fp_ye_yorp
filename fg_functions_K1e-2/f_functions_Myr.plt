#!/usr/bin/gnuplot

set xl "obliquity gamma [deg]"
set yl "f-function [10^-4 s^-1/Myr]"

set yr [-2:2]
set nokey
set zeroaxis
set xtics 30
set mxtics 3

f(x) = x/1e-4 * (1.e6*365.25*86400.)

p "<./intep_test  1.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test  2.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test  3.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test  4.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test  5.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test  6.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test  7.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test  8.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test  9.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 10.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 11.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 12.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 13.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 14.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 15.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 16.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 17.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 18.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 19.y 0 180 1" u 1:(f($2)) w l,\
  "<./intep_test 20.y 0 180 1" u 1:(f($2)) w l
pa -1

set term png small
set out "f_functions_Myr.png"
rep


