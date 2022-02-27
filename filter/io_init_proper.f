c**********************************************************************
c IO_INIT_PROPER.F
c**********************************************************************
c Read parameters of proper-elements filter. There are 5 variants
c of the filter:
c
c 1) Running averages of mean elements
c 2) Minimum/maximum filter for `resonant elements'
c 3) Running minimum/maximum of individual elements
c 4) Representative surface filter (another `res. el.')
c 5) Frequency Modified Fourier Transform
c    ref.: Sidlichovsky & Nesvorny (1997), Cel. Mech. 65, 137
c
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: May 9th 2003
c
c Revisions:
c Nov 18th 2003: additonal parameters in proper-filter input file,
c   namely required range of g and s frequencies and possible
c   (vectorial) addition of output frequencies
c Jul 25th 2006: added line numbers to error messages for easier
c   debugging of proper.in
c Aug 23rd 2006: an error message when open fails was added

      subroutine io_init_proper(infile,dtproper,outfile,iflgchk)

      include '../swift.inc'
      include 'filter.inc'
      include 'cb_flt.inc'
      include 'proper.inc'

c  input
      character*(*) infile

c  output
      integer iflgchk
      real*8 dtproper
      character*(*) outfile

c temporal variables
      integer i, j, k, iu, ierr, l
      real*8 dtoutfilter

c  functions
      integer read_i
      logical read_l
      real*8 read_f

c  executable code...
      iu = 7

c  read input parameters
      write(*,*) 'Proper parameters file called ',infile
      call io_open(iu,infile,'old','formatted',ierr)
      if (ierr.ne.0) then
        write(*,*) 'io_init_proper: Error opening file ',infile,'.'
        call util_exit(1)
      endif

      l=0
      call read_s(iu,outfile,l)
      prop_ftype = read_i(iu,l)
      prop_writer8 = read_l(iu,l)

      dtproper = read_f(iu,l)

      prop_win = read_f(iu,l)
      prop_win = prop_win * 365.25d6
      prop_dt = read_f(iu,l)
      prop_dt = prop_dt * 365.25d6

      prop_ielmt = read_i(iu,l)
      prop_minmax = read_l(iu,l)      

      prop_mmel(1) = read_l(iu,l)
      prop_mmel(2) = read_l(iu,l)
      prop_mmel(3) = read_l(iu,l)

      prop_plid = read_i(iu,l)
      prop_p = read_i(iu,l)
      prop_q = read_i(iu,l)
      prop_sigma = read_f(iu,l)
      prop_sigma = prop_sigma / degrad
      prop_dsigma = read_f(iu,l)
      prop_dsigma = prop_dsigma / degrad
      prop_varpi = read_f(iu,l)
      prop_varpi = prop_varpi / degrad
      prop_dvarpi = read_f(iu,l)
      prop_dvarpi = prop_dvarpi / degrad
      prop_Omega = read_f(iu,l)
      prop_Omega = prop_Omega / degrad
      prop_dOmega = read_f(iu,l)
      prop_dOmega = prop_dOmega / degrad
      prop_sigmadot = read_i(iu,l)

      fmft_flag = read_i(iu,l)
      fmft_nfreq = read_i(iu,l)
      fmft_ndata = read_i(iu,l)
      fmft_writegs = read_l(iu,l)
      fmft_minf = read_f(iu,l)
      fmft_maxf = read_f(iu,l)
      fmft_freqsmall = read_f(iu,l)
      fmft_freqtoler = read_f(iu,l)
      fmft_nfreqknown = read_i(iu,l)

      if (fmft_nfreqknown.gt.MAX_FMFT_NFREQKNOWN) then
        write(*,40) fmft_nfreqknown,MAX_FMFT_NFREQKNOWN
40      format('Number of known frequencies ',l6,' have to be .lt. ',l6,
     :    ' (see proper.inc file)')
        call util_exit(1)
      endif
      do i = 1, fmft_nfreqknown
        fmft_freqknown(i) = read_f(iu,l)
      enddo

      fmft_ming = read_f(iu,l)
      fmft_maxg = read_f(iu,l)
      fmft_mins = read_f(iu,l)
      fmft_maxs = read_f(iu,l)
      fmft_vectadd = read_l(iu,l)
      fmft_chkrst = read_l(iu,l)
      fmft_alsopl = read_l(iu,l)
      fmft_propgs = read_l(iu,l)

      close(unit = iu)
      write(*,*) ' '

c  check input parameters and available array dimensions

c  filter type
      if ((prop_ftype.lt.0).or.(prop_ftype.gt.7)) then
        write(*,50) prop_ftype
50      format('Unknown proper filter type no. ',i2,'.')
        call util_exit(1)
      endif

c  sampling time-step
      if (dtproper.lt.0.d0) then
        dtproper = dtfilter
      endif

c  give output filename also in proper common block
      outproperfile = outfile

c  mean elements output time-step
      k = 1
      do i = 1, filter_N
        k = k * filter_D(i)
      enddo
      dtoutfilter = k * dtfilter
      j = int(prop_win/dtoutfilter)
      if ((j.gt.PMAXN).and.(prop_ftype.le.3)) then
        write(*,20) prop_win/365.25d6,dtoutfilter/365.25d3,j,PMAXN
20      format('Running window width over mean element output (',f5.2,
     :    ' Myr/',f6.3,' kyr = ',i6,') have to be .lt. ',i6,
     :    ' (see proper.inc file)')
        call util_exit(1)
      endif

      if (fmft_ndata.gt.PMAXN) then
        write(*,30) fmft_ndata,PMAXN
30      format('Width of FMFT buffer ',i6,' have to be .lt. ',i6,
     :    ' (see proper.inc file)')
        call util_exit(1)
      endif

      return
      end       ! io_init_proper.f

c-----------------------------------------------------------------


