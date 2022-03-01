#!/usr/bin/gnuplot

M_S = 2.9591397931969080E-04
M_J = 2.8253459095242259E-07

RHSCALE = 3.5  # see rmvs.inc

load "config.plt"

set xl "t [My]"
set yl "q, Q [au]"
set cbl "TP"

set yr [3:5.5]
set zeroaxis
set key left width -1

set palette rgbformulae 33,13,10

p \
  "<awk '($1>0)' ../test_cybele_12_jupiter_e0.15_OUTWRD/follow.out" u 2:($3*(1+$4)):1 w p lc 'gray',\
  "<awk '($1>0)' follow.out" u 2:($3*(1+$4)):1 w p lc 'black' t "Q",\
  "<awk '($1==-2)' follow.out" u 2:($3*(1-$4)):1 w d lc 'cyan' t "q_J",\
  "<awk '($1==-2)' follow.out" u 2:($3*(1-$4)*(1.-(M_J/(3.*M_S))**(1./3.))):1 w d lc 'green' t "Hill",\
  "<awk '($1==-2)' follow.out" u 2:($3*(1-$4)*(1.-RHSCALE*(M_J/(3.*M_S))**(1./3.))):1 w d lc 'orange' t "RHSCALE"

pa -1

set term png small
set out "qt.png"
rep


