#!/bin/sh

FILES=`cat files.txt`

for FILE in $FILES ; do
  OUT=`echo $FILE | awk '{ gsub("\\\\.", "_2ND."); print; }'`
  echo "../filter/$FILE -> $OUT"

  ./2nd.awk ../filter/$FILE > $OUT

done

exit


