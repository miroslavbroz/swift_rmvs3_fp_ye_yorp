#!/usr/bin/awk -f

BEGIN{

  f_slow = 1./ARGV[1]; ARGV[1]="";
#  f_slow = 1./150.;	# yr^-1; I generally do NOT want frequencies larger than this limit
			# (ALL secular frequencies are definitely lower!)

  passband_A = 0.005;	# frequencies f/f_sampling < than this are passed (with ripple_A)
  stopband_A = 0.05;	#             f/f_sampling >               stopped (with suppression_A)

  passband_B = 0.05;	# the same for filter B
  stopband_B = 0.15;

  ripple_A   = 3.5e-5;	# maximum distortion of the signal in the passband of filter A
  suppress_A = 1.e-5;	# minimum suppression of the signal in the stopband

  ripple_B   = 1.e-9;	# the same for filter B
  suppress_B = 1.e-9;

  print "========================================================================";

  print "f_slow = " f_slow " yr^-1 = 1/" 1./f_slow " yr";
  print "<- lower frequencies (larger periods) than this have to be PRESERVED!\n";

  print "passband_A <=> f/f_sampling < " passband_A;
  print "stopband_A <=>              > " stopband_A;
  print "passband_B <=> f/f_sampling < " passband_B;
  print "stopband_B <=>              > " stopband_B;

  print "ripple_A = " ripple_A "	, suppress_A = " suppress_A;
  print "ripple_B = " ripple_B "	, suppress_B = " suppress_B;

  print "========================================================================";

  total_ripple = 0.;
}
{
  gsub("d", "e");
}
(FNR==1){
  nfilters = $1;
}
(FNR==2){
  for (i=1; i<=length($1); i++) {
    filter[i] = substr($1,i,1);
  }
}
(FNR==3){
  gsub("!.*", "");
  for (i=1; i<=NF; i++) {
    K[i] = $i;
  }
}
(FNR==6){
  dt[0] = $1/365.25;
}

END{
  for (i=1; i<=nfilters; i++) {
    dt[i] = FILTER(i,dt[i-1],filter[i],K[i]);
  }

  print "\n========================================================================";
  print "final sampling timestep dt = " dt[nfilters] " yr";
  print "\nfinal passband <=> f < " f_min0 " yr^-1 <=> P > " 1./f_min0 " yr";
  print "total ripple in the passband = " total_ripple;
  print "\nfinal stopband <=> f > " f_max0 " yr^-1 <=> P < " 1./f_max0 " yr";
  print "minimum suppression in the stopband = " suppress "\n";
}

########################################################################

func FILTER(i,dt0,F,K){

  f_sampling0 = 1./dt0;
  print "dt" i-1 " = " dt0 " yr, f_sampling" i-1 " = " f_sampling0 " yr^-1";

  print "filter no. " i " = " F;
  print "decimation factor K" i " = " K;

  if (F == "A") { passband = passband_A; stopband = stopband_A; ripple = ripple_A; suppress = suppress_A; }
  else          { passband = passband_B; stopband = stopband_B; ripple = ripple_B; suppress = suppress_B; }

## passband check
  print "f_slow = " f_slow " yr^-1 <= passband_" F " * f_sampling" i-1 " = " passband * f_sampling0 " yr^-1";

  if (f_slow <= passband * f_sampling0) {
    print "Passband check passed OK."
  } else {
    print "ERROR: PASSBAND CHECK NOT PASSED!!!\n"
    exit 1;
  }

  f_min0 = f_sampling0 * passband;
  f_max0 = f_sampling0 * stopband;
  print "\nf_max" i-1 " = " f_max0 " yr^-1 = 1/" 1/f_max0 " yr";

  dt1 = K * dt0;
  print "dt" i " = " dt1 " yr";

## aliasing check
  f_sampling1 = 1./dt1;
  print "f_sampling" i " = " f_sampling1 " yr^-1 >= f_slow + f_max" i-1 " = " f_slow + f_max0;

  if (f_sampling1 >= f_slow + f_max0) {
    print "Aliasing check passed OK."
  } else {
    print "ERROR: ALIASING CHECK NOT PASSED!!!\n"
    exit 1;
  }

  total_ripple = total_ripple + ripple  ;

  print "\n------------------------------------------------------------------------";

  return dt1;
}


