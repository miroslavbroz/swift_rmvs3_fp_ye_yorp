c*************************************************************************
c CB_FLT.INC
c*************************************************************************
c Include file with /filter/ common block: various parameters
c of filtration process:
c
c  filter_A(),filter_B()	filters A,B
c  dtfilter			the same as dtfilter_param
c  filter_MA,filter_MB		dimension of above filters
c  filter_N			number of filter/decimation cycles
c  filter_D(filter_N)		decimation parameter for each cycle
c  FBUFN			number of samples to store in memory
c  FNinibuf			# of data before the 1st output (initialization)
c  FNdtout			output time-step
c  filter_direct		direct e/inc filtering?
c  filter_write			write binary output file?
c  filter_write_r8		write real*8 binary output file instead of real*4?
c  filter_stat(filter_N)	char array, which filter to use in the cycle?
c  outfilterfile		output filename
c
c Author: Miroslav Broz
c Date: Apr 8th 2003

      real*8 filter_A(-FMAXM:FMAXM),filter_B(-FMAXM:FMAXM),dtfilter
      integer filter_MA,filter_MB,filter_N,filter_D(FMAXC),FBUFN,
     &  FNinibuf,FNdtout
      logical filter_direct,filter_write,filter_write_r8
      integer filter_plid,filter_p,filter_q
      character*10 filter_stat  ! have to be the last in /filter/
      character*80 outfilterfile

      common/filter/filter_A,filter_B,dtfilter,filter_MA,filter_MB,
     &  filter_N,filter_D,FBUFN,FNinibuf,FNdtout,filter_direct,
     &  filter_write,filter_write_r8,filter_plid,filter_p,filter_q,
     &  filter_stat,outfilterfile


