c**********************************************************************
c PROPER.INC
c**********************************************************************
c Include file for proper-elements subroutines.
c
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Apr 8th 2003
c
c Revisions:
c Nov 18th 2003: additonal parameters in proper-filter input file,
c   namely required range of g and s frequencies and possible
c   (vectorial) addition of output frequencies

c  this requires `filter.inc' to be included also

      integer PMAXN, PMAXE, PMAXNPL, PMAXNTP
      parameter (PMAXE = 6)		! max # of elements to store
      parameter (PMAXNPL=FMAXNPL)	! max # of PLs and TPs is the same as for mean-elmts filter
      parameter (PMAXNTP=FMAXNTP)

c  this have to be EXACTLY the same as in fmft.c!!!
      parameter (PMAXN=2048)		! max # of samples in buffer, have to be >=2^N due to FMFT!

c  other parameteres in input file (and common block /proper/)
      integer prop_ftype, prop_ielmt, prop_plid, prop_p, prop_q,
     :  prop_sigmadot
      real*8 prop_win, prop_dt, prop_sigma, prop_dsigma, prop_varpi,
     :  prop_dvarpi, prop_Omega, prop_dOmega
      logical prop_writer8, prop_minmax, prop_mmel(3)

c  FMFT filter parameters
      integer MAX_FMFT_NFREQ, MAX_FMFT_NDATA, MAX_FMFT_FLAG,
     :  MAX_FMFT_NFREQKNOWN

c  these 2 params have to be the same as in fmft.c!!!
      parameter (MAX_FMFT_NFREQ = 100)		! max # of calculated frequencies
      parameter (MAX_FMFT_NDATA = PMAXN)	! max # of data points in buffer

      parameter (MAX_FMFT_NFREQKNOWN = 20)	! max # of (known) plan. frequencies (2 be discarded)
      parameter (MAX_FMFT_FLAG = 3)		! size of FMFT output buffer

      integer fmft_flag, fmft_nfreq, fmft_ndata, fmft_nfreqknown
      real*8 fmft_minf, fmft_maxf
      real*8 fmft_freqknown(MAX_FMFT_NFREQKNOWN), fmft_freqtoler,
     &  fmft_freqsmall
      real*8 fmft_ming, fmft_maxg, fmft_mins, fmft_maxs
      logical fmft_writegs, fmft_vectadd, fmft_chkrst, fmft_alsopl,
     &  fmft_propgs
      character*80 outproperfile

c   large-type first due to a 4-byte padding:
      common /proper/
     &  prop_win, prop_dt,
     &  prop_sigma, prop_dsigma, prop_varpi,
     &  prop_dvarpi, prop_Omega, prop_dOmega,
     &  fmft_minf, fmft_maxf, fmft_freqtoler, fmft_freqsmall,
     &  fmft_ming, fmft_maxg, fmft_mins, fmft_maxs, fmft_freqknown,
     &  prop_ftype,
     &  prop_writer8,
     &  prop_ielmt, prop_minmax, 
     &  prop_mmel,
     &  prop_plid, prop_p, prop_q, prop_sigmadot,
     &  fmft_flag, fmft_nfreq, fmft_ndata, fmft_nfreqknown,
     &  fmft_writegs, fmft_vectadd, fmft_chkrst, fmft_alsopl,
     &  fmft_propgs,
     &  outproperfile


