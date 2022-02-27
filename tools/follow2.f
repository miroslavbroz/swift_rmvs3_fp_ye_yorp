c converts swift binary (integer,real,real*8) file to an ascii file
c write all PLs and TPs: id, time, a, e, inc, peri, node, m
c Miroslav Broz, miroslav.broz@email.cz, Apr 8th 2003

      program follow2

c-----------------------------------------------------------------------
      implicit none

c input params
      character*80 outfile,infile
      integer ftype

c temp vars
      integer i,id,istep,ierr,nbod,nleft
      real*8 t,tmax,dr,a,e,inc,capom,omega,capm,apo,peri
      character*80 str,fmt

c functions
      integer io_read_hdr,io_read_line
      integer io_read_hdr_r,io_read_line_r
      integer io_read_hdr_r8,io_read_line_r8

c constants
      real*8 pi
      parameter(pi=3.1415926535d0)

c-----------------------------------------------------------------------

c read input params from stdin
      if (iargc().eq.0) then
        read(*,5) infile
5       format(a)
        read(*,5) outfile
        read(*,*) ftype

c check command-line arguments
      else if (iargc().eq.3) then
        call getarg(1, infile)
        call getarg(2, outfile)
        call getarg(3, str)
        read(str,*,err=13,end=13) ftype
        goto 14
13      continue
        write(*,*) 'Usage: follow2 infile outfile <type [0-2]>'
        stop
14      continue
      else
15      continue
        write(*,*) 'Usage: follow2 infile outfile <type [0-2]>'
        stop
      endif 
      
      if (ftype.eq.2) then
c        write(*,*) ' Reading an real*8 binary file '
      else if (ftype.eq.1) then
c        write(*,*) ' Reading an real*4 binary file '
      else if (ftype.eq.0) then
c        write(*,*) ' Reading an integer*2 binary file '
      else
        write(*,*) 'follow2: Unknown binary file format specified.'
        stop
      endif

c  output format depend on input precision
      if (ftype.eq.2) then
        fmt = '(1x,i8,7(1x,f22.16))'
      else
        fmt = '(1x,i8,1x,f16.10,1x,f14.8,1x,f10.8,4(1x,f10.5))'
      endif

      tmax = 0
      dr = 180.0/pi

      open(unit=1,file=infile,status='old',form='unformatted',
     :  iostat=ierr)
      if (ierr.ne.0) then
        write(*,*) 'follow2: Error opening file ',infile
        stop
      endif
c      open(unit=2,file=outfile,status='new',iostat=ierr)
      open(unit=2,file=outfile,status='unknown',iostat=ierr)
      if (ierr.ne.0) then
        write(*,*) 'follow2: Error opening file ',outfile
        stop
      endif

10    continue
        istep = 0
        if (ftype.eq.2) then
          ierr = io_read_hdr_r8(1,t,nbod,nleft)
        else if (ftype.eq.1) then
          ierr = io_read_hdr_r(1,t,nbod,nleft)
        else if (ftype.eq.0) then
          ierr = io_read_hdr(1,t,nbod,nleft)
        endif
        if (ierr.ne.0) goto 20
        t = t/365.25e6	! Myr
        do i=2,nbod
          if (ftype.eq.2) then
            ierr = io_read_line_r8(1,id,a,e,inc,capom,omega,capm)
          else if (ftype.eq.1) then
            ierr = io_read_line_r(1,id,a,e,inc,capom,omega,capm)
          else if (ftype.eq.0) then
            ierr = io_read_line(1,id,a,e,inc,capom,omega,capm)
          endif
          if (ierr.ne.0) goto 20
          istep = 1
          inc = inc*dr
          capom = capom*dr
          omega = omega*dr
          capm = capm*dr
c          peri = a*(1.0d0-e)
c          apo = a*(1.0d0+e)
c          write(2,1000) id,t,a,e,inc
c1000      format(1x,i4,1x,f16.10,1x,f14.8,1x,f10.8,1x,f10.5)
          write(2,fmt) id,t,a,e,inc,capom,omega,capm
          tmax = t
        enddo

        do i=1,nleft
          if (ftype.eq.2) then
            ierr = io_read_line_r8(1,id,a,e,inc,capom,omega,capm)
          else if (ftype.eq.1) then
            ierr = io_read_line_r(1,id,a,e,inc,capom,omega,capm)
          else if (ftype.eq.0) then
            ierr = io_read_line(1,id,a,e,inc,capom,omega,capm)
          endif
c          write(*,*) 'ierr = ', ierr
c          write(*,*) id,t,a,e,inc,capom,omega,capm
          if (ierr.ne.0) goto 20
          istep = 1
          inc = inc*dr
          capom = capom*dr
          omega = omega*dr
          capm = capm*dr
c          peri = a*(1.0d0-e)
c          apo = a*(1.0d0+e)
c          write(2,1000) id,t,a,e,inc
          write(2,fmt) id,t,a,e,inc,capom,omega,capm
          tmax = t
        enddo
        if (istep.eq.0) goto 20
      goto 10
20    continue

      write(*,*) 'follow2: t_max = ',tmax,' Myr'

      stop
      end

