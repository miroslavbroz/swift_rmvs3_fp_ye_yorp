#!/bin/sh

# filter.sh
# Calculate passbands/stopbands of the filters in filter.in.
# See the algorithm in Quinn et al. (1991).
# Miroslav Broz (miroslav.broz@email.cz), Jun 3rd 2008

FILTERIN=filter.in
if [ ! -e $FILTERIN ] ; then FILTERIN=""; fi

./filter.awk 2000 filter.in  >  filter.out
./filter.awk 150  filter2.in >> filter.out

cat filter.out


