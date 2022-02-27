c yarko_omega.f
c Precalculate Yarkovsky parameters, which depend on omega.
c Miroslav Broz (miroslav.broz@email.cz), Jan 19th 2010

      subroutine yarko_omega(ntp,istat)

      include '../swift.inc'
      include 'const.inc'
      include 'yarko.inc'
      include 'spin.inc'

c ... Inputs:
      integer ntp
      integer istat(NTPMAX,NSTAT)

      integer i
      real*8 mass,Gamma,ls,Rdash,theta
      complex*8 sinz,cosz,z

      do i = 1,ntp
        if (istat(i,1).eq.0) then

          mass = 4./3.*pi*R(i)**3*rho_bulk(i)
          Gamma = sqrt(rho_surf(i)*C_th(i)*K_th(i))
          ls = sqrt(K_th(i)/(rho_surf(i)*C_th(i)*abs(omega(i))))
          Rdash = R(i)/ls
          Rdash0(i) = R(i)/sqrt(K_th(i)/(rho_surf(i)*C_th(i)))

c  cancelled factor exp(-imag(z)) in fraction psi(i) to avoid numerical problems
          z = 1/sqrt(2.)*(1.,-1.)*Rdash
          sinz = (cos(dble(z))+(0.,1.)*sin(dble(z))
     &      -(cos(-dble(z))+(0.,1.)*sin(-dble(z)))*exp(2.*imag(z)))
     &      /(0.,2.)
          cosz = (cos(dble(z))+(0.,1.)*sin(dble(z))
     &      +(cos(-dble(z))+(0.,1.)*sin(-dble(z)))*exp(2.*imag(z)))
     &      /2.
          psi(i) = ((z**2-3.)*sinz+3.*z*cosz)/(sinz-z*cosz)

          lambda0(i) = Gamma*sqrt(abs(omega(i)))/(eps_IR(i)*sigma
     :      *sqrt(2.)*Rdash)
          theta = sqrt(2.)*Rdash*lambda0(i)
          lambdas0(i) = Gamma/(eps_IR(i)*sigma*sqrt(2.))
          Tstar40(i) = (1.-A_Bond(i))*S0/(eps_IR(i)*sigma)
          rho0(i) = -4./9.*eps_IR(i)*sigma*pi*R(i)**2/(mass*cv)
     &      *day**2/AU

          thetas(i) = Gamma/(eps_IR(i)*sigma*Tstar40(i)**(3./4.))/
     &      sqrt(day/k_gauss)
          rhos0(i) = pi*(R(i)**2)*S0/(mass*cv)*(1.-A_Bond(i))

        endif  ! istat
      enddo    ! ntp

      return
      end      ! yarko_omega.f

c---------------------------------------------------------------------


