c yarko.inc
c Include file for Yarkovsky routines.
c Miroslav Broz (miroslav.broz@email.cz), Jan 19th 2010

      integer MAXk_seasonal, k_seasonal            
      parameter(MAXk_seasonal = 7)      ! max number of iterations in calculation of the seasonal Yarkovsky effect
      parameter(k_seasonal = 7)		! actual number of interations to use

c  thermal parameters
      real*8 R(NTPMAX), rho_bulk(NTPMAX), rho_surf(NTPMAX),
     :  K_th(NTPMAX), C_th(NTPMAX), A_Bond(NTPMAX), eps_IR(NTPMAX)

c  pre-calculated Yarkovsky parameters
      real*8 lambda0(NTPMAX), Tstar40(NTPMAX), rho0(NTPMAX),
     :  lambdas0(NTPMAX), Rdash0(NTPMAX), thetas(NTPMAX), rhos0(NTPMAX),
     :  koefsin(MAXk_seasonal,NTPMAX), koefcos(MAXk_seasonal,NTPMAX)
      double complex psi(NTPMAX)

c  input/output thermal parameters in common block /yarko/
      common /yarko/
     :  R, rho_bulk, rho_surf, K_th, C_th, A_Bond, eps_IR,
     :  psi, lambda0, Tstar40, rho0, lambdas0, Rdash0,
     :  thetas, rhos0, koefsin, koefcos


