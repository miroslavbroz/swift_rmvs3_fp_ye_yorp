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
c  /proper_2ND/
c  /elements_2ND/
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
c Aug 14th 2008: added new (7th) filter for Trojans

      subroutine io_write_proper_2ND(t,nbod,ntp,istat,oname,iu,
     :  fopenstat)

      include '../swift.inc'
      include 'filter.inc'
      include 'proper_2ND.inc'
      include 'cb_oscel_2ND.inc'
      include 'cb_meanel_2ND.inc'

      real*8 t
      integer nbod, ntp, iu
      integer istat(NTPMAX,NSTAT)
      character*(*) oname, fopenstat

c  temporal variables
      integer i1st,ierr
      real*8 tstart,tstop,eps

      save i1st,tstart,tstop

c  data
      data i1st,eps/0,1.d-8/

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
          write(*,*) ' SWIFT ERROR: in io_write_proper_2ND: '
          write(*,*) '     Could not open binary output file:'
          call util_exit(1)
        endif
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

        call proper_runavg_2ND(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.2) then

        call proper_minmax_2ND(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.3) then

        call proper_minmax2_2ND(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.4) then

        call proper_sigma_2ND(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.6) then

        call proper_sigma2_2ND(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.5) then

        call proper_fmft_2ND(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      elseif (prop_ftype.eq.7) then

        call proper_trojan_2ND(t,tstart,tstop,nbod,ntp,istat,oname,iu)

      endif     ! prop_ftype.eq.6

      return
      end


