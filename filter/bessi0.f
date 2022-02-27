c****************************************************************
      double precision function bessi0(x)
c****************************************************************
c
c  Modified Bessel function
c  from Numerical Recipes, p. 273, function bessi0
c  http://cfatab.harvard.edu/nr/bookf.html
c
      implicit none
      double precision x
      double precision ax,p1,p2,p3,p4,p5,p6,p7,
     &  q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      data p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     &  1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     &  0.225319d-1,-0.157565d-1,0.916281d-2,-0.2057706,
     &  0.2635537d-1,-0.1647633d-1,0.392377d-1/

      if (abs(x).lt.3.75) then
        y=(x/3.75)**2
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=abs(x)
        y=3.75/ax
        bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4
     &     +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      endif
      return
      end

