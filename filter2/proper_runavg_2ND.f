
      subroutine proper_runavg_2ND(t,tstart,tstop,nbod,ntp,istat,
     :  oname,iu)
c
c  running average filter
c  (it produces exactly the same results as off-line version mean2prop)
c
      include '../swift.inc'
      include 'filter.inc'
      include 'proper_2ND.inc'
      include 'cb_flt_2ND.inc'
      include 'cb_meanel_2ND.inc'
      include 'cb_propel_2ND.inc'

      real*8 t,tstart,tstop
      integer nbod,ntp,iu,istat(NTPMAX,NSTAT)
      character*(*) oname

      integer PN
      save PN

c  temporal variables
      integer i,id,ierr
      real*8 tout,capa,e,inc

c  functions
      real*8 arr_avg

      data PN /1/

c=======================================================================
c  ... executable code

      prop_time(PN) = t

      do i = 2, nbod
        id = -i
        prop_elmts(PN,1,id) = mean_elmts(1,id)
        prop_elmts(PN,2,id) = mean_elmts(6,id)
        prop_elmts(PN,3,id) = mean_elmts(7,id)
      enddo

      do i = 1, ntp
        if (istat(i,1).eq.0) then
          prop_elmts(PN,1,i) = mean_elmts(1,i)
          prop_elmts(PN,2,i) = mean_elmts(6,i)
          prop_elmts(PN,3,i) = mean_elmts(7,i)
        endif
      enddo

c  the buffer is filled-up => write the output
      if (t.gt.tstop) then

        tout = t - 0.5d0*prop_win - 0.5d0*FNinibuf*dtfilter

        call io_open(iu,oname,'append','UNFORMATTED',ierr)
        if (prop_writer8) then
          call io_write_hdr_r8(iu,tout,nbod,ntp,istat)
        else
          call io_write_hdr_r(iu,tout,nbod,ntp,istat)
        endif

        do i = 2, nbod
          id = -i
          capa = arr_avg(prop_elmts(1,1,id), PN)
          e    = arr_avg(prop_elmts(1,2,id), PN)
          inc  = arr_avg(prop_elmts(1,3,id), PN)
          if (prop_writer8) then
            call io_write_line_r8(iu,id,capa,e,inc,0.0d0,0.0d0,0.0d0)
          else
            call io_write_line_r(iu,id,capa,e,inc,0.0d0,0.0d0,0.0d0)
          endif
        enddo

        do i = 1, ntp
          if (istat(i,1).eq.0) then
            capa = arr_avg(prop_elmts(1,1,i), PN)
            e    = arr_avg(prop_elmts(1,2,i), PN)
            inc  = arr_avg(prop_elmts(1,3,i), PN)
            if (prop_writer8) then
              call io_write_line_r8(iu,i,capa,e,inc,0.0d0,0.0d0,0.0d0)
            else
              call io_write_line_r(iu,i,capa,e,inc,0.0d0,0.0d0,0.0d0)
            endif
          endif
        enddo

        close(iu)

c  shift the prop_elmts buffer back in time
        call proper_shift_2ND(tstart,tstop,PN,nbod,ntp,
     :    prop_time,prop_elmts,3)

      endif

      PN = PN + 1

      return
      end


