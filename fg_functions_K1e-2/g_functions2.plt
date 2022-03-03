#!/usr/bin/gnuplot

set xl "obliquity gamma [deg]"
set yl "g-function [rad/s]"
set zeroaxis
set grid noxtics noytics front

p "<awk '(FNR==1){ print s; }{ print; }' *.y" u 1:3 w lp

pa -1

set term png small
set out "g_functions2.png"
rep


