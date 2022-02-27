c io_close.f
c Close the files, which remain open all the time.
c Miroslav Broz (miroslav.broz@email.cz), Feb 12th 2008

      subroutine io_close(iub,iuf,iup)

      integer iub,iuf,iup

      close(iub)
      close(iuf)
      if (iup.gt.0) then
        close(iup)

c see proper_fmft.f
        close(iup+1)
        close(iup+2)
      endif

      return
      end

