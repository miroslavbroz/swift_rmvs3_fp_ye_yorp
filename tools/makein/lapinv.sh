#!/bin/sh

if [ "$1" = "" ] ; then
  echo "Usage: lapinv.sh [FILE]"
  echo "FILE is formatted like follow2.out, laplac.out must be in the CWD"
fi

awk '(NR==4){Ls2=$3; Ls3=$4;} (NR==8){gm=$1; print gm,Ls2,Ls3;}' laplac.out > lapinv.tmp
cat lapinv.tmp $@ | lapinv > lapinv.out

