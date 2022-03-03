#!/usr/bin/awk -f

# astorb_makein.awk
# Convert astorb.dat to makein input format
# (=> tp.in and pl.in for numerical integrators).
# Miroslav Broz (miroslav.broz@email.cz), Mar 12th 2003

BEGIN{
  n=0;
  TMP="astorb_makein.tmp";
}
!/^#/{
  if (n==0) {
    od=substr($0,107,8);
    year=substr(od,1,4);
    month=substr(od,5,2);
    day=substr(od,7,2);
    CMD = "jdate.pl jd " year " " month " " day;
    CMD | getline jd;
  }
  n++;
  a=substr($0,169,14);		# this is for NEW AstOrb format
  e=substr($0,158,12);
  i=substr($0,148,11);
  peri=substr($0,126,12);
  node=substr($0,137,12);
  m=substr($0,115,12);
  print a " " e " " i " " peri " " node " " m > TMP;
}
END{
  fflush("");

  print jd "\n11\nT\n(6(1x,e22.16))\n" n;
  system("cat " TMP);
  system("rm " TMP);
}


