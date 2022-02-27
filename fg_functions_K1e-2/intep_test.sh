#!/bin/sh

awk -vFIRST=$5 -vSECOND=$6 '!/^#/{ print $FIRST, $SECOND; }' < $1 > intep_test.tmp
./intep_test intep_test.tmp $2 $3 $4

