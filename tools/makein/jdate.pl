#!/usr/bin/perl

if ($#ARGV < 0) {
  print "Usage: jdate.pl [jd YYYY MM DD|od JD]\n";
  exit;
}

if ($ARGV[0] eq "jd") {
  print &jdate($ARGV[1], $ARGV[2], $ARGV[3]) . "\n";
} else {
  ($y, $m, $d) = &od($ARGV[1]);
  print "$y $m $d\n";
}

exit;

sub jdate {
  local($y,$m,$d) = @_;
  local($yy,$mm,$jd,$a,$b);

  if ($m > 2) {
    $yy = $y;
    $mm = $m;
  } else {
    $yy = $y - 1.0;
    $mm = $m + 12.0;
  }
  $jd = int(365.25*$yy) + int(30.6001*($mm+1)) + $d + 1720994.5;
  if ($y > 1582 || $y == 1582 && $m > 10 || $y == 1582 && $m == 10 && int($d) > 15) {
    $a = int($yy/100.0);
    $b = 2.0-$a+int($a/4.0);
    $jd = $jd + $b;
  }
  return $jd;
}

sub od {
  local($jd) = @_;
  local($y,$m,$d,$z,$f,$alfa,$a,$b,$c,$x,$e);

  $z = int($jd+0.5);
  $f = &frac($jd+0.5);
  if ($z<2299163.0) {
    $a = $z;
  } else {
    $alfa = int(($z-1867216.25)/36524.25);
    $a = $z+1.0+$alfa-int($alfa/4.0);
  }
  $b = $a + 1524.0;
  $c = int(($b-122.1)/365.25);
  $x = int(365.25*$c);
  $e = int(($b-$x)/30.6001);
  $d = $b-$x-int(30.6001*$e)+$f;
  if ($e<14.0) {
    $m = $e - 1.0;
  } else {
    $m = $e - 13.0;
  }
  if ($m>2.0) {
    $y = $c - 4716.0;
  } else {
    $y = $c - 4715.0;
  }
  return ($y,$m,$d);
}

sub frac {
  local($x) = @_;
  return ($x-int($x));
}

