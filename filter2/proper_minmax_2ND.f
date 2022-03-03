
      subroutine proper_minmax_2ND(t,tstart,tstop,nbod,ntp,istat,
     :  oname,iu)
c
c  running minimum (or maximum) filter
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
      integer i,j,id,ierr
      real*8 tout,capa,e,inc,node,peri,varpi

c  functions
      integer arr_minind,arr_maxind
      real*8 zero2pi

      data PN /1/

c=======================================================================
c  ... executable code

      prop_time(PN) = t

      do i = 2, nbod
        id = -i
        prop_elmts(PN,1,id) = mean_elmts(1,id)
        prop_elmts(PN,2,id) = mean_elmts(6,id)
        prop_elmts(PN,3,id) = mean_elmts(7,id)
        varpi = atan2(mean_elmts(3,id), mean_elmts(2,id))
        node = atan2(mean_elmts(5,id), mean_elmts(4,id))
        peri = zero2pi(varpi-node)
        prop_elmts(PN,4,id) = node
        prop_elmts(PN,5,id) = peri
      enddo

      do i = 1, ntp
        if (istat(i,1).eq.0) then
          prop_elmts(PN,1,i) = mean_elmts(1,i)
          prop_elmts(PN,2,i) = mean_elmts(6,i)
          prop_elmts(PN,3,i) = mean_elmts(7,i)
          varpi = atan2(mean_elmts(3,i), mean_elmts(2,i))
          node = atan2(mean_elmts(5,i), mean_elmts(4,i))
          peri = zero2pi(varpi-node)
          prop_elmts(PN,4,i) = node
          prop_elmts(PN,5,i) = peri
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

c  decide whether use minimum or maxium of the orbital element no. prop_ielmt
        if (prop_minmax) then
          do i = 2, nbod
            id = -i
            j = arr_minind(prop_elmts(1,prop_ielmt,id), PN)
            capa = prop_elmts(j, 1, id)
            e    = prop_elmts(j, 2, id)
            inc  = prop_elmts(j, 3, id)
            node = prop_elmts(j, 4, id)
            peri = prop_elmts(j, 5, id)
            if (prop_writer8) then
              call io_write_line_r8(iu,id,capa,e,inc,peri,node,0.0d0)
            else
              call io_write_line_r(iu,id,capa,e,inc,peri,node,0.0d0)
            endif
          enddo

          do i = 1, ntp
            if (istat(i,1).eq.0) then
              j = arr_minind(prop_elmts(1,prop_ielmt,i), PN)
              capa = prop_elmts(j, 1, i)
              e    = prop_elmts(j, 2, i)
              inc  = prop_elmts(j, 3, i)
              node = prop_elmts(j, 4, i)
              peri = prop_elmts(j, 5, i)
              if (prop_writer8) then
                call io_write_line_r8(iu,i,capa,e,inc,peri,node,0.0d0)
              else
                call io_write_line_r(iu,i,capa,e,inc,peri,node,0.0d0)
              endif
            endif
          enddo

c  maximum
        else
          do i = 2, nbod
            id = -i
            j = arr_maxind(prop_elmts(1,prop_ielmt,id), PN)
            capa = prop_elmts(j, 1, id)
            e    = prop_elmts(j, 2, id)
            inc  = prop_elmts(j, 3, id)
            node = prop_elmts(j, 4, id)
            peri = prop_elmts(j, 5, id)
            if (prop_writer8) then
              call io_write_line_r8(iu,id,capa,e,inc,peri,node,0.0d0)
            else
              call io_write_line_r(iu,id,capa,e,inc,peri,node,0.0d0)
            endif
          enddo

          do i = 1, ntp
            if (istat(i,1).eq.0) then
              j = arr_maxind(prop_elmts(1,prop_ielmt,i), PN)
              capa = prop_elmts(j, 1, i)
              e    = prop_elmts(j, 2, i)
              inc  = prop_elmts(j, 3, i)
              node = prop_elmts(j, 4, i)
              peri = prop_elmts(j, 5, i)
              if (prop_writer8) then
                call io_write_line_r8(iu,i,capa,e,inc,peri,node,0.0d0)
              else
                call io_write_line_r(iu,i,capa,e,inc,peri,node,0.0d0)
              endif
            endif
          enddo
        endif
        close(iu)

c  shift the prop_elmts buffer back in time
        call proper_shift_2ND(tstart,tstop,PN,nbod,ntp,
     :    prop_time,prop_elmts,5)

      endif

      PN = PN + 1

      return
      end
