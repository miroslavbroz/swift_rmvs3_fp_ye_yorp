#!/usr/bin/gnuplot

set xl "obliquity gamma [deg]"
set yl "f-function [rad/s]"
set zeroaxis
set grid noxtics noytics front

p "<awk '(FNR==1){ print s; }{ print; }' *.y" u 1:2 w lp

pa -1

set term png small
set out "f_functions2.png"
rep


