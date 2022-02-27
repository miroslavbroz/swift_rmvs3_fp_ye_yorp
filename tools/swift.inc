c*************************************************************************
c                        SWIFT.INC
c*************************************************************************
C Include file for SWIFT
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 3/7/93
c
c Modified: Miroslav Broz, miroslav.broz@email.cz
c Date: Oct 10th 2001


        implicit NONE    ! you got it baby

c...   Version of Swift
       real*8 VER_NUM
       parameter(VER_NUM=2.0d0)

c...   Maximum array size
       integer  NPLMAX, NTPMAX
c       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun 
       parameter  (NPLMAX = 51)   ! max number of planets, including the Sun 
c       parameter  (NTPMAX = 1001) ! max number of test particles
       parameter  (NTPMAX = 30001) ! max number of test particles

c...   Size of the test particle integer status flag
	integer NSTATP            ! Number of status parameters
        parameter  (NSTATP = 3)  
	integer NSTAT            ! Number of status parameters
        parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ planet  

c...   Size of the test particle integer status flag
	integer NSTATR    
        parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR

c...   convergence criteria for danby
	real*8 DANBYAC , DANBYB
	parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

c...    loop limits in the Laguerre attempts
	integer NLAG1, NLAG2
	parameter(NLAG1 = 50, NLAG2 = 400)

c...    A small number
	real*8 TINY
        PARAMETER(TINY=4.D-15)

c...    trig stuff
	real*8 PI,TWOPI,PIBY2,DEGRAD
	parameter (PI = 3.14159265358979D0)
	parameter (TWOPI = 2.0D0 * PI)
	parameter (PIBY2 = PI/2.0D0)
	parameter (DEGRAD = 180.0D0 / PI)

c  number of steps in step_kdk_pl2, size of PL coordinates buffer
      integer INTMAX
      parameter (INTMAX = 3)

c...  Coefficients of the 2nd order integrator of type SBAB
      real*8 c_koef(INTMAX-1),d_koef(INTMAX)
      data c_koef /0.5d0,
     :             0.5d0/
      data d_koef /0.16666666666666666667d0,
     :             0.66666666666666666667d0,
     :             0.16666666666666666667d0/

c  functions
      integer ibset
      logical btest

c-------------------------------------------------------------------------

