c  common block for osculating elements (passed from io_write_filter_2ND.f
c  /filter_elmts_write_2ND/ subroutine to getacc_yarko.f, yarko_seasonal.f
c  and io_write_proper_2ND.f subroutines)

      real*8 osc_elmts(6,-FMAXNPL:FMAXNTP)
      real*8 t0
      logical upflg_osc
      common/elements_2ND/osc_elmts,t0,upflg_osc

