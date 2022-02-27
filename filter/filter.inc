c*************************************************************************
c FILTER.INC
c*************************************************************************
C Include file for filtering routines
c
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Oct 23rd 2001

        integer FMAXM,FMAXC,FMAXN,FMAXE,FMAXNPL,FMAXNTP
        parameter(FMAXM=80)		! max size of LF filter -FMAXM:FMAXM
	parameter(FMAXC=4)		! max # of filter/decimation cycles
	parameter(FMAXN=2*FMAXM+1)	! max # of samples in buffer
	parameter(FMAXE=9)		! max # of elements to be filtered
c
c  max # of planets and test particles, which orbital elements will be
c  filtered (we need to allocate array elmts(FMAXN,7,-FMAXNPL:FMAXNTP,FMAXC))
c
	parameter(FMAXNPL=NPLMAX)
	parameter(FMAXNTP=NTPMAX)

c 2DO: CALL A PREPROCESSOR!!!
c  there are no TPs in SyMBA:
c	parameter(FMAXNPL=NTPMAX)
c	parameter(FMAXNTP=1)

 
