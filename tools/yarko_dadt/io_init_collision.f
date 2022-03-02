c**********************************************************************
c IO_INIT_COLLISION.F
c**********************************************************************
c Read collision-related parameters
c
c Input:
c  infile	name of the yarko-data input file
c  nbod		number of massive bodies (including the Sun)
c  ntp		number of test particles
c
c  in common block /yarko/
c  R(NTPMAX)		radii of TP's from io_init_yarko.f subroutine
c
c Output:
c  dtreorient		reorientation time step
c  tau_reor(NTPMAX)	mean reorientation times
c  tau_disr(NTPMAX)	mean disruption times
c  iflgchk		if bits 6 or 7 are set, then write *.out file
c
c  in common block /hill/
c  hill_scale		individual scaling of Hill spheres
c
c Remarks: 
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Apr 18th 2001

      subroutine io_init_collision(infile,nbod,ntp,
     &  dtreorient,dtdisrupt,tau_reor,tau_disr,beta_1,
     &  omega_1,omega_2,a_mean1,a_mean2,iflgchk)

      include 'swift.inc'
      include 'yarko.inc'

c  input
      character*(*) infile
      integer nbod,ntp,iflgchk

c  output
      real*8 dtreorient,dtdisrupt,tau_reor(NTPMAX),tau_disr(NTPMAX),
     :  beta_1,omega_1,omega_2,a_mean1,a_mean2

c  Hill spheres scaling factors (see modified rmvs3_chk.f subroutine)
      real*8 hill_scale(NPLMAX)
      common /hill/ hill_scale

c  internal
      integer i,ierr
      real*8 tau_reor0,tau_disr0,B,beta_2,D_0,omega_0,P_1,P_2
      logical write_reorient,write_disrupt

c  main ----------------------------------------------------------
      write(*,*) 'Collisional data file called ',infile
      call io_open(7,infile,'old','formatted',ierr)

      read(7,*,err=99,end=99) dtreorient
      read(7,*,err=99,end=99) B
      read(7,*,err=99,end=99) beta_1
      read(7,*,err=99,end=99) beta_2
      read(7,*,err=99,end=99) D_0
      read(7,*,err=99,end=99) omega_0
      read(7,*,err=99,end=99) dtdisrupt
      read(7,*,err=99,end=99) tau_disr0
      read(7,*,err=99,end=99) (hill_scale(i),i=2,nbod)
      read(7,*,err=99,end=99) P_1, P_2
      read(7,*,err=99,end=99) a_mean1, a_mean2
      read(7,*,err=99,end=99) write_reorient
      read(7,*,err=99,end=99) write_disrupt

c  dt's are given in yr
      dtreorient = dtreorient*365.25d0
      dtdisrupt  = dtdisrupt *365.25d0

c  periods are in hours
      omega_1 = 2.d0*PI/(P_1*3600.d0)
      omega_2 = 2.d0*PI/(P_2*3600.d0)

      do i=1,ntp

c  Expecting dependence tau_reor(omega, D) = B (omega/omega_0)^beta_1 (D/D_0)^beta_2
c  Beware! We have to multiply by omega^beta_1 later on due to YORP evolution!

        tau_reor(i) = 365.25d0*B * (1.0d0/omega_0)**beta_1
     :    * (2.d0*R(i)/D_0)**beta_2

c  expecting dependence of disruption time scale of the form: tau*sqrt(R [m])
        tau_disr(i) = 365.25d0*tau_disr0 * sqrt(R(i))
      enddo

      if (write_reorient) iflgchk=ibset(iflgchk,6)
      if (write_disrupt) iflgchk=ibset(iflgchk,7)

      close(unit = 7)
      write(*,*) ' '
      return

99    continue
      write(*,*) 'Error reading file',infile
      stop

      end    ! io_init_collision.f
c-----------------------------------------------------------------


