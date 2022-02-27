#!/bin/sh

awk '($1==0){ print $2; }' *.y > f_hist.dat

