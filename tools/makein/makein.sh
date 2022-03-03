#!/bin/sh

# makein.sh
# Create pl.in and tp.in input files (x's and v's) for SWIFT from orbital elements.
# Planetary ephemeris is taken from JPL's DE405.
# Miroslav Broz (miroslav.broz@email.cz), May 28th 2008

#./astorb_makein.awk librating.astorb > makein.in	# this is for AstOrb-like input file
./astdysequ_makein.awk *.equ > makein.in		# this is for *.eq1 files from AstDyS

DIR=.
PATH=$DIR:$PATH

ln -sf $DIR/JPLEPH		# this is for another directory...
makein < makein.in > makein.out

# apply barycentric correction for terestrial planets (and Pluto)
baryc2 10 < makein.out > baryc2.out

# rotate into invariant plane (suitable for studies of resonant asteroids)
laplac < baryc2.out > laplac.out

# split the pl.in and tp.in files and calculate *.elmts to check the results
pltp.awk < laplac.out
mv pl.in pl_inv.in
mv tp.in tp_inv.in
xvpl2el < pl_inv.in > pl_inv.elmts
xvtp2el < tp_inv.in > tp_inv.elmts

pltp.awk < baryc2.out
xvpl2el < pl.in > pl.elmts
xvtp2el < tp.in > tp.elmts


