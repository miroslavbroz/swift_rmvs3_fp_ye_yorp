c spin.inc
c Include file for spin-related routines.
c Miroslav Broz (miroslav.broz@email.cz), Jan 19th 2010

c  spin axes data
      real*8 s(3,NTPMAX), omega(NTPMAX)

c  input/output thermal parameters in common block /spin/
      common /spin/
     :  s, omega


