#!/usr/bin/gawk -f

BEGIN{
  FS=",";
  print "# D P PFlag";
  print "# km h 1";
}
(i==1){
  gsub("\"", "");
  D = $9;
  P = $19;
  A = $23;
  U = $24;
  if ((D != "-9.99") && (P !="-9.99") && (U+0==3)) {
    print D,P,A,U;
  }
}
/# Data/{
  i=1;
}




