#!/bin/sh

make

./makein < makein.in > makein.out
awk 'BEGIN{ print 1; } ((NR>=3) && (NR<=5)) || (NR>=33)' < makein.out > makein.tmp
./pltp.awk < makein.tmp

#./baryc2 2 3 4 5 6 7 8 9 10 < makein.out > baryc.out
#./pltp.awk < baryc.out

./yarko_dadt < yarko_dadt.in


