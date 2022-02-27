
c****************************************************************
      subroutine fmft_gsout(output, id)
c****************************************************************
c
c  move full output of the fmft_call subroutine to a 3-D array, wrt. TP or PL id
c
      include '../swift.inc'
      include 'filter.inc'
      include 'proper.inc'

c  g and s frequencies, amplitudes and phases from fmft_call subroutine
      real*8 gsoutput(MAX_FMFT_NFREQ,3)
      common/gs_output/gsoutput

c  this definition should be the same as in proper_fmft subroutine
      real*8 output(-PMAXNPL:PMAXNTP,MAX_FMFT_NFREQ,3)
      integer id

      integer j,k

      if (fmft_propgs) then     ! well, there are many repetitious checks of this :-(
        do j = 1, fmft_nfreq
          do k = 1, 3
            output(id,j,k) = gsoutput(j,k)
          enddo
        enddo
      endif

      return
      end


