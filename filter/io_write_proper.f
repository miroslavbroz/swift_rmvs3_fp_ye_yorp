c*************************************************************************
c IO_WRITE_PROPER
c*************************************************************************
c accumulate array of orbital elements for test particles and planets,
c filter it (by the given method) and write output to a real*4 binary file.
c
c Input:
c  t         ==> current time (real scalar)
c  nbod      ==> number of massive bodies (int scalar)
c  ntp       ==> number of test particles (int scalar)
c  istat     ==>  status of the test paricles
c  oname     ==> output file name (character string) 
c  iu        ==> unit number to write to
c  fopenstat ==> the status flag for the open statements
c
c other parameters or data in common blocks:
c  /proper/
c  /elements/
c  /elements_mean/
c
c Output:
c  stored in internal arrays and then written to binary file `oname'
c
c Authors: Miroslav Broz, miroslav.broz@email.cz
c Date:    Apr 23rd 2003
c
c Revisions:
c Nov 18th 2003: additonal parameters in proper-filter input file,
c   namely required range of g and s frequencies and possible
c   (vectorial) addition of output frequencies
c Aug 25th 2006: adapted to SyMBA (which handles PLs only and no TPs)

      subroutine io_write_proper(t,nbod,ntp,istat,oname,iu,fopenstat)

      include '../swift.inc'
      include 'filter.inc'
      include 'proper.inc'
      include 'cb_oscel.inc'
      include 'cb_meanel.inc'

      real*8 t
      integer nbod, ntp, iu
      integer istat(NTPMAX,NSTAT)
      character*(*) oname, fopenstat

c  temporal variables
      integer i1st,ierr
      real*8 tstart,tstop,eps

      save i1st,tstart,tstop

c  data
      data i1st,eps/0,1.d-3/

c=======================================================================
c  ... executable code

c  return immediately, if the filter is switched off or
c  mean or osculating elements were NOT updated!
      if ((prop_ftype.eq.0).or.
     :  (.not.((upflg_mean).or.
     :  ((upflg_osc).and.(prop_ftype.eq.4)))))
     :  then
        return
      endif
      upflg_mean = .false.
      upflg_osc = .false.

c  1st time through (only a single test every call)
      if (i1st.eq.0) then

c  try to open outfile (right at the beginning of the integration)
        call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
        if (ierr.ne.0) then
          write(*,*) ' SWIFT ERROR: in io_write_proper: '
          write(*,*) '     Could not open binary output file:'
          call util_exit(1)
        endif
c modified by Miroslav Broz (miroslav.broz@email.cz), Feb 12th 2008
        close(iu)

c  initial values
        tstart = t
        tstop = tstart + prop_win
        i1st = 1
      endif
c
c  possibly skip some amount of data (ie. if twin.lt.dt)
c
      if (t.lt.tstart-eps) then
        return
      endif
c
c  decide on filter type
c
      if (prop_ftype.eq.1) then

        call proper_runavg(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.2) then

        call proper_minmax(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.3) then

        call proper_minmax2(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.4) then

        call proper_sigma(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.6) then

        call proper_sigma2(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.5) then

        call proper_fmft(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.7) then

        call proper_orbfit(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      endif     ! prop_ftype.eq.6

      return
      end


