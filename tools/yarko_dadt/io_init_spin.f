c**********************************************************************
c IO_INIT_SPIN.F
c**********************************************************************
c Read in spin axes data.
c
c Input:
c  infile	name of the spin axes data input file
c  ntp		number of test particles
c
c Output:
c  iseed		number to initialize ran1()
c
c  in common block /yarko/
c  s(3,NTPMAX)		orientation of the spin axis	
c  omega(NTPMAX)	spin rates
c
c Remarks: 
c Author:  Miroslav Broz, miroslav.broz@email.cz
c Date:  Jan 19th 2009

      subroutine io_init_spin(infile,ntp,iseed)

      include 'swift.inc'
      include 'spin.inc'

c  input
      character*(*) infile
      integer ntp

c  output
      integer iseed

c  internal
      integer i,ierr,ntpchk,iform
      real*8 s1,s2

c  main     
      write(*,*) 'Spin axes data file called ',infile
      call io_open(7,infile,'old','formatted',ierr)

      read(7,*,err=99,end=99) ntpchk
      if(ntp.eq.0) ntp = ntpchk
      if(ntpchk.ne.ntp) then
        write(*,*) ' SWIFT ERROR: in io_init_spin'
        write(*,*) '     Number of test particles differs in ',infile
        call util_exit(1)
      endif

      read(7,*,err=99,end=99) iseed

      read(7,*,err=99,end=99) iform
      if (iform.eq.0) then
        do i=1,ntpchk
          read(7,*,err=99,end=99) s1,s2,omega(i)
          s1=s1/DEGRAD
          s2=s2/DEGRAD
          s(1,i)=cos(s1)*cos(s2)
          s(2,i)=sin(s1)*cos(s2)
          s(3,i)=sin(s2)
        enddo
      elseif (iform.eq.1) then
        do i=1,ntpchk
          read(7,*,err=99,end=99) s(1,i),s(2,i),s(3,i),omega(i)
        enddo
      else
        write(*,*) ' SWIFT ERROR: in io_init_spin'
        write(*,*) '     Spin axes format ',iform,' not implemented.'
        call util_exit(1)
      endif

      close(unit = 7)
      write(*,*) ' '
      return

99    continue
      write(*,*) 'Error reading file',infile
      stop

      end    ! io_init_spin.f
c-----------------------------------------------------------------


