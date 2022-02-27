#!/usr/bin/gnuplot

set xl "obliquity gamma [deg]"
set yl "g-function [rad/s]"

set nokey
set zeroaxis
set xtics 30
set mxtics 3

p "<./intep_test.sh   1.y 0 180 1  1 3",\
  "1.y" u 1:3
pa -1

set term png small
set out "g_functions1.png"
rep


