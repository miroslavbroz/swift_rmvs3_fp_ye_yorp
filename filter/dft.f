c dft.f
c Discrete Fourier Transform, amplitude of a given frequency.
c Miroslav Broz (miroslav.broz@email.cz), Aug 18th 2008

      subroutine dft(t, x, n, f, ampl, phase)

      implicit none

      integer n
      real*8 t(n), x(n), f, ampl, phase

      integer i
      real*8 A_f, B_f

      A_f = 0.d0
      B_f = 0.d0

      do i = 1, n
        A_f = A_f + 2.d0*x(i) * cos(f*t(i))
        B_f = B_f + 2.d0*x(i) * sin(f*t(i))
      enddo

      A_f = A_f/n
      B_f = B_f/n

      ampl = sqrt(A_f**2 + B_f**2)
      phase = -atan2(B_f, A_f)

      return
      end


