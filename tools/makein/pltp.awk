#!/usr/bin/awk -f

# Split makein output to pl.in and tp.in files.
# Eventually drop Pluto (ie. the last planet) - "-pluto" switch.
# Miroslav Broz (miroslav.broz@email.cz), Mar 13th 2003

BEGIN{
  PLIN = "pl.in";
  TPIN = "tp.in";

  if (ARGV[1] == "-pluto") {
    ARGV[1] = "";
    ntp_minus = -1;
  } else {
    ntp_minus = 0;
  }

  n = 0;
  ntp = 0;
}
!/^#/{
  n++;
  if (n == 1) {
    ntp = $1;
    print (ntp+ntp_minus) > PLIN;
  } else if ((n > 1) && (n <= (ntp+ntp_minus)*3+1)) {
    print $0 > PLIN;
  } else if (n > ntp*3+1) {
    print $0 > TPIN;
  }
}

