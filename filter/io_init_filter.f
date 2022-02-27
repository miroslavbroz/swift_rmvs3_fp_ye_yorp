c**********************************************************************
c IO_INIT_FILTER.F
c**********************************************************************
c Read in filter parameters, precompute filters.
c
c Filter/decimation process was published in Quinn, T. R., et al.:
c Integration of Earth's orbit. AJ 102, p. 2287, 1991
c
c These subroutines and functions are used to obtain long term
c variations of orbital elements.
c
c Input:
c  infile	name of the filter parameters input file,
c		for structure of this file refer to example filter.in
c  ntp		number of TP's (due to check the buffer size)
c
c Output:
c  dtfilter			sampling period
c  outfile			name of binary output file
c  iflgchk			real*8 and real*4 output flags
c
c  in common block /filter/
c  see filter_cb.inc
c
c Remarks: 
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Apr 8th 2003
c
c Revisions:
c Aug 23rd 2006: an error message when open fails was added,
c   both filters A and B can have different dimensions now

      subroutine io_init_filter(infile,nbod,ntp,dtfilter_param,
     :  outfile,iflgchk)

      include '../swift.inc'
      include 'filter.inc'
      include 'cb_flt.inc'

c  input
      character*(*) infile
      integer nbod,ntp

c  output
      real*8 dtfilter_param
      character*(*) outfile
      integer iflgchk

c  internals
      integer i,ierr,AM,BM
      real*8 Ax0,Bx0,Abeta,Bbeta
      logical write_r8

c  functions
      integer length

c  main ----------------------------------------------------------

c Check number of test particles, filtration subroutines needs to allocate
c buffer for orbital elements elmts(FMAXN,6,FMAXNTP)
c Filter parameters are given to io_write_filter.f in /filter/ common block

      if (ntp.gt.FMAXNTP) then
        write(*,*) 'Number of bodies nbod =',ntp,' exceeds',
     &    ' max number FMAXNTP =',FMAXNTP, '.'
        call util_exit(1)
      endif

c  read input parameters
      write(*,*) 'Filter parameters file called ',infile
      call io_open(7,infile,'old','formatted',ierr)
      if (ierr.ne.0) then
        write(*,*) 'io_init_filter: Error opening file ',infile,'.'
        call util_exit(1)
      endif

      read(7,*) filter_N
      read(7,10) filter_stat
10    format(a)
      read(7,*) (filter_D(i),i=1,filter_N)
      read(7,*) filter_MA,Ax0,Abeta
      read(7,*) filter_MB,Bx0,Bbeta
c
c  sampling period, output to outfile will be produced every
c  filter_D(1) * filter_D(2) * ... * filter_D(N) * dtfilter, see FBUFN
c
      read(7,*) dtfilter
      read(7,10) outfile
      read(7,*) filter_direct
      read(7,*) write_r8
      read(7,*) filter_write
      read(7,*) filter_write_r8
      read(7,*) filter_plid
      read(7,*) filter_p
      read(7,*) filter_q

      close(unit = 7)
      write(*,*) ' '

c  check the number of filters
      if (filter_stat(length(filter_stat):length(filter_stat))
     :  .eq.char(13)) then
        filter_stat = filter_stat(1:length(filter_stat)-1)
        write(*,20) infile(1:length(infile))
20      format('io_init_filter: Warning: DOS-like EOLs in ', a)
      endif
      if (filter_N.ne.length(filter_stat)) then
        write(*,*) 'Number of filters ',filter_N,' differs',
     &    ' from length of filter sequence ',
     &    filter_stat(1:length(filter_stat)), '.'
        call util_exit(1)
      endif

c  chelck PL id
      if ((filter_plid.gt.-2).or.(filter_plid.lt.-FMAXNPL)) then
        write(*,*) 'io_init_filter: Error: filter_plid is outside',
     &    ' interval [-FMAXNPL, -2].'
        call util_exit(1)
      endif

      call filter_create(filter_A,filter_MA,Ax0,Abeta)
      call filter_create(filter_B,filter_MB,Bx0,Bbeta)

      FBUFN = 2*max(filter_MA,filter_MB)+1      ! width of buffers

c  give dtfilter both as a parameter and in common block 
      dtfilter_param = dtfilter
c  give outfile also in flt_cb common block
      outfilterfile = outfile

c  setup iflgchk bits, according to options
c  the summary of iflgchk bits follows:
c     6 ... reorientation
c     7 ... disruption
c     8 ... write real*8 bin.dat

      if (write_r8) iflgchk=ibset(iflgchk,8)

c  FNdtout * dtfilter ... output (filtered) timestep
c  (ie. multiplied decimation factors) 
      FNdtout = 1
      do i = 1, filter_N
        FNdtout = FNdtout * filter_D(i)
      enddo

c  FNinibuf ... number of points needed to fill the buffer INITIALLY
c  (we will have to shift output time by FNinibuf * dtfilter / 2)
      FNinibuf = FBUFN
      do i = filter_N-1, 1, -1
        FNinibuf = filter_D(i)*(FNinibuf-1) + FBUFN
      enddo 
      FNinibuf = FNinibuf - 2   ! empirically "measured" factor?!
c  2DO: compute FNinibuf correctly even for a single filter!

      return
      end    ! io_init_filter
c-----------------------------------------------------------------

