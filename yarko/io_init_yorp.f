c io_init_yorp.f
c Read YORP parameters from file and initialize f- and g-functions.
c Miroslav Broz (miroslav.broz@email.cz), Jan 19th 2010

      subroutine io_init_yorp(infile,ntp,dt,dtout)

      include '../swift.inc'
      include 'yorp.inc'

c  input (plus variables in /yorp/ common block!)
      character*(*) infile
      integer ntp
      real*8 dt, dtout

c  internal
      integer i,j,ierr,ntpchk,id
      real*8 dobliq,dummy,tmp
      character*80 fgfile,fgdir,str

c  functions
      integer length

c  main
      write(*,*) 'YORP data file called ',infile
      call io_open(7,infile,'old','formatted',ierr)

      read(7,*) ntpchk
      if (ntp.eq.0) ntp = ntpchk
      if (ntpchk.ne.ntp) then
        write(*,*) 'io_init_yorp: Error: Number of test particles ',
     :    'differs in ', infile
        call util_exit(1)
      endif

      read(7,*) dt
      read(7,*) dtout
      read(7,*) nGAUSS
      read(7,*) nGDATA
      read(7,*) dobliq
      read(7,*) a_0
      read(7,*) R_0
      read(7,*) rho_0
      read(7,*) omega_0
      read(7,*) c_YORP
      read(7,*) gauss_rnd
      read(7,10) fgdir
10    format(a)

      do i = 1,ntp
        read(7,*) id, tmp
        if (id.gt.ntp) then
          write(*,*) 'io_init_yorp: Error: ID for a test particle ',
     :      'is larger than NTP = ', ntp
          call util_exit(1)
        endif
        fg_id(id) = tmp
      enddo

      close(7)

      write(*,*) ' '
      if (nGAUSS.gt.GAUSSMAX) then
        write(*,*) 'io_init_yorp: Error: number of Gaussian spheres ',
     :    'nGAUSS.gt.GAUSSMAX = ', GAUSSMAX
        stop
      endif

      if (nGDATA.gt.GDATAMAX) then
        write(*,*) 'io_init_yorp: Error: number of f- and g-functions ',
     :    'data nGDATA.gt.GDATAMAX = ', GDATAMAX
        stop
      endif

      dt = dt*365.25d0
      dtout = dtout*365.25d0
      dobliq = dobliq/DEGRAD

      do i = 0, nGDATA-1
        obliq_func(i) = 0.d0 + dobliq*i
      enddo
c
c  read f- and g-functions from directory fg_functions/
c
      do j = 1,nGAUSS

        write(fgfile,*) j
        fgfile = trim(fgdir) // trim(adjustl(fgfile)) // '.y'
        open(unit=7, file=fgfile, status='old', form='formatted',
     :    iostat=ierr)
        if (ierr.ne.0) then
          write(*,*) 'io_init_yorp: Error opening file ', fgfile
          stop
        endif

        i = 0
5       continue
          read(7, 10, err=30,end=30) str
          if (str(2:2).ne.'#') then     !  drop comments
            read(str,*,err=30,end=30) dummy, f_func(i,j), g_func(i,j),
     :        dummy
            if (abs(g_func(i,j)).lt.1.d-21) then       ! set exact zeros for too small values!
              g_func(i,j) = 0.d0
            endif
            i=i+1
          endif
        goto 5

30      continue
        close(7)

        if (i.ne.nGDATA) then
          write(*,*) 'io_init_yorp: Error: number of data in file ',
     :      fgfile(1:length(fgfile)), ' differs from nGDATA = ', nGDATA
          stop
        endif
 
      enddo   ! nGAUSS

      return
      end     ! io_init_yorp.f

c-----------------------------------------------------------------


