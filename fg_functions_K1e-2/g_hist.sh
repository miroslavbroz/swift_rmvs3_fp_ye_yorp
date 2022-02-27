#!/bin/sh

awk '($1==0){ print $3; }' *.y > g_hist.dat

