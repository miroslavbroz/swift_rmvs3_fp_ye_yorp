C
      SUBROUTINE INTEP(XP,P,X,F,N,IER)
c      SUBROUTINE INTEP(N,IER)
C      Purpose:
C     To interpolate a function value P for a given argument value XP
C        using a table of N values (X,F).
C     This is a Spline Interpolation scheme based on Hermite polynomials
C     The source is U.S. Airforce Surveys in Geophysics No 272.
C
C     Usage:
C     For random values of XP
C      CALL INTEP (XP,P,X,F,N,IER)
C     or after the first call to INTEP with monotonically increasing
C     or decreasing values of XP consistent with the X vector
C      CALL EINTEP (XP,P,X,F,N,IER)
C
C     Description of parameters:
C     XP   The chosen argument value.
C     P    The resultant interpolated value.
C     X    The vector of independent values.
C     F    The vector of function or dependent values.
C     N    The number of points in the (X,P) vectors.
C     IER  The resultant error parameter.
C     Remarks:
C     If XP is beyond either extreme in the vector X the value of F
C           at that extreme is adopted and IER set ot 2.
C
c      IMPLICIT REAL*8 (A-H,O-Z)
c      REAL*8 LP1,LP2,L1,L2
c      DIMENSION F(1),X(1)
      implicit none
      integer N,IER
      real*8 XP,P,X(N),F(N)
      integer I,IO,IUP,N1
      real*8 LP1,LP2,L1,L2,FP1,FP2,XPI,XPI1

      IER=1
      IO=1
      IUP=0
      IF (X(2).LT.X(1)) IUP=1
      N1=N-1
      IF ((XP.GE.X(N).AND.IUP.EQ.0).OR.(XP.LE.X(N).AND.IUP.EQ.1)) THEN
    5   P=F(N)
        IER=2
        GO TO 6
      ELSE IF ((XP.LE.X(1).AND.IUP.EQ.0).OR.
     :  (XP.GE.X(1).AND.IUP.EQ.1)) THEN
        P=F(1)
    6   IER=2
      END IF
      if (IER.eq.2) then
c        write(*,*) 'intep: Warning: IER = ', IER
        return
      endif
      ENTRY EINTEP (XP,P,X,F,N,IER)
      DO 1 I=IO,N
        IF (XP.LT.X(I).AND.IUP.EQ.0) GO TO 2
        IF (XP.GT.X(I).AND.IUP.EQ.1) GO TO 2
    1 CONTINUE
      GO TO 5
    2 I=I-1
      IF (I.EQ.IO-1) GO TO 4
      IO=I+1
      LP1=1.D0/(X(I)-X(I+1))
      LP2=1.D0/(X(I+1)-X(I))
      IF (I.EQ.1) FP1=(F(2)-F(1))/(X(2)-X(1))
      IF (I.EQ.1) GO TO 3
      FP1=(F(I+1)-F(I-1))/(X(I+1)-X(I-1))
    3 IF (I.GE.N1) FP2=(F(N)-F(N-1))/(X(N)-X(N-1))
      IF (I.GE.N1) GO TO 4
      FP2=(F(I+2)-F(I))/(X(I+2)-X(I))
    4 XPI1=XP-X(I+1)
      XPI=XP-X(I)
      L1=XPI1*LP1
      L2=XPI*LP2
      P=F(I)*(1D0-2D0*LP1*XPI)*L1*L1+F(I+1)*(1D0-2D0*LP2*XPI1)
     :  *L2*L2+FP2*XPI1*L2*L2+FP1*XPI*L1*L1
      RETURN
      END

