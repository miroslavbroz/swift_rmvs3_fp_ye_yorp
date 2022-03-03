#!/usr/bin/awk -f

# astdysequ_makein.awk
# Convert *.eq1 files from Pisa (AstDyS) to makein input format
# (=> tp.in and pl.in for numerical integrators).
# Miroslav Broz (miroslav.broz@email.cz), May 29th 2008

BEGIN{
  pi=3.1415926535; pi2=2.*pi; deg=pi/180.; rad=180./pi;
  n=0;
  TMP="astdysequ_makein.tmp";
  OFMT=CONVFMT="%22.16f";
}
/^ MJD/{
  mjd = $2;
  jd = mjd + 2400000.5
}
/^ EQU/{
  n++;
  a = $2+0;
  esinlp = $3;		# varpi = peri + node
  ecoslp = $4;
  tani2sinln = $5;
  tani2cosln = $6;
  meanlong = $7*deg;	# lambda = M + varpi

  e     = sqrt(esinlp**2 + ecoslp**2);
  varpi = atan2(esinlp, ecoslp);
  i     = 2. * atan(sqrt(tani2sinln**2 + tani2cosln**2));
  node  = atan2(tani2sinln, tani2cosln);
  peri  = zero2pi(varpi - node);
  m     = zero2pi(meanlong - node - peri);

  i    = i * rad;
  peri = peri * rad;
  node = node * rad;
  m    = m * rad;

  print a " " e " " i " " peri " " node " " m > TMP;
}
END{
  fflush("");

  print jd "\n11\nT\n(6(1x,e22.16))\n" n;
  system("cat " TMP);
  system("rm " TMP);
}

func tan(x) {
  return sin(x)/cos(x);
}

func atan(x) {
  return atan2(x,1.);
}

func zero2pi(x, y) {
  y=x/pi2
  y=(y-int(y))*pi2
  if (y<0.) { y=y+pi2; }
  return y;
}


