c**********************************************************************
c IO_INIT_YARKO.F
c**********************************************************************
c Read in thermal data, precompute thermal parameters for getacch subroutine
c
c Input:
c  infile	name of the thermal data input file
c  ntp		number of test particles
c
c Output:
c  in common block /yarko/
c
c  in common block /prdrag/
c  pr(NTPMA_BondX)		precomputed parameter of P-R drag
c
c Remarks: 
c A_Bonduthor: Miroslav Broz, miroslav.broz@email.cz
c Date: Feb 12th 2002

      subroutine io_init_yarko(infile,ntp)

      include 'swift.inc'
      include 'const.inc'
      include 'yarko.inc'
      include 'pr.inc'
      include 'spin.inc'

c  input
      character*(*) infile
      integer ntp

c  internal
      integer i,ierr,ntpchk

c  main
      write(*,*) 'Yarkovsky data file called ',infile
      call io_open(7,infile,'old','formatted',ierr)

      read(7,*,err=99,end=99) ntpchk
      if(ntp.eq.0) ntp = ntpchk
      if(ntpchk.ne.ntp) then
        write(*,*) ' SWIFT ERROR: in io_init_yarko'
        write(*,*) '     Number of test particles differs in ',infile
        call util_exit(1)
      endif

      do i = 1,ntp
        read(7,*,err=99,end=99) R(i),rho_bulk(i),rho_surf(i),K_th(i),
     :    C_th(i),A_Bond(i),eps_IR(i)

c  Precompute parameters also for getacch_prdrag subroutine.
c  (This is useful for combined *_mvs2fypr integrator; see also io_init_pr
c  - there is a similar line, except units of R and rho are different!)

        pr(i) = 3.d0/4.d0* S0/(cv**2*rho_bulk(i)*R(i)) *day	! note unit conversion!

c  in case omega(i).eq.-1, scale omega as 1/R
        if (omega(i).lt.-0.99d0) then
          omega(i) = TWOPI/(5.d0*3600.d0)/(2.d0*R(i)/1000.d0)
        endif

      enddo  ! ntp

      close(unit = 7)
      write(*,*) ' '
      return

99    continue
      write(*,*) 'Error reading file',infile
      stop

      end    ! io_init_yarko.f
c-----------------------------------------------------------------


