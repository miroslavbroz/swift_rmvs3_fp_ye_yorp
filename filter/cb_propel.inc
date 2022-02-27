
c  large buffer for averaging (dimensions are: width of the
c  window, number of elements and PLs+TPs);
c  the same definition should be also in proper_shift subroutine

      real*8 prop_elmts(PMAXN,PMAXE,-PMAXNPL:PMAXNTP)
      real*8 prop_time(PMAXN)

      common/proper_elements/prop_elmts,prop_time

