#!/usr/bin/gnuplot

set xl "obliquity gamma [deg]"
set yl "g-function [deg/Myr]"

set nokey
set zeroaxis
set xtics 30
set mxtics 3

# navic musim zapocitat nasobeni omega <=> P = 6 h!
f(x) = x/pi*180. * (1.e6*365.25*86400.) / (2.*pi/(6.*3600.))

p "<./intep_test.sh   1.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh   2.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh   3.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh   4.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh   5.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh   6.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh   7.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh   8.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh   9.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  10.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  11.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  12.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  13.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  14.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  15.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  16.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  17.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  18.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  19.y 0 180 1  1 3" u 1:(f($2)) w l,\
  "<./intep_test.sh  20.y 0 180 1  1 3" u 1:(f($2)) w l
pa -1

set term png small
set out "g_functions_Myr.png"
rep


