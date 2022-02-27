c yorp.inc
c Include file for YORP routines.
c Miroslav Broz (miroslav.broz@email.cz), Jan 19th 2010

      integer GAUSSMAX, GDATAMAX
      parameter(GAUSSMAX = 200)		! number of Gaussian spheres
      parameter(GDATAMAX = 180)		! number of data points in f_ and g_functions/ files

c  YORP parameters
      logical gauss_rnd
      integer nGAUSS, nGDATA, fg_id(NTPMAX)
      real*8 c_YORP, a_0, R_0, rho_0, omega_0, obliq_func(0:GDATAMAX),
     :  f_func(0:GDATAMAX, GAUSSMAX), g_func(0:GDATAMAX, GAUSSMAX)

c  input/output YORP parameters in common block /yorp/
      common /yorp/
     :  c_YORP, a_0, R_0, rho_0, omega_0, obliq_func, f_func, g_func,
     :  nGAUSS, nGDATA, fg_id, gauss_rnd


