c**********************************************************************
      integer function length(s)
c**********************************************************************
c
c  length of string s without spaces
c
      implicit none

      character s*(*)

      integer i

      do 10 i=len(s),1,-1
        if (s(i:i).ne.' ') goto 20
10    continue
20    length=i
      return
      end

