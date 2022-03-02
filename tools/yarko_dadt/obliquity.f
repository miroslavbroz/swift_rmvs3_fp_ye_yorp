
      real*8 function obliquity(a1,a2,a3,b1,b2,b3,s1,s2,s3)

      real*8 a1,a2,a3,b1,b2,b3,s1,s2,s3,obliq
      real*8 axb1,axb2,axb3,norm

      axb1 = (a2*b3 - a3*b2)
      axb2 = (a1*b3 - a3*b1)
      axb3 = (a1*b2 - a2*b1)

      norm = sqrt(axb1**2 + axb2**2 + axb3**2)
      obliquity = acos((s1*axb1 + s2*axb2 + s3*axb3) / norm)

      return
      end

