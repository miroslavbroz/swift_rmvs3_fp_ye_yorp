c*************************************************************************
c                            UTIL_VERSION.F
c*************************************************************************
c Prints version of Swift
c
c             NO IO!!!!!!
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/21/94
c Last revision: 
c
c Modified: Miroslav Broz, miroslav.broz@email.cz
c Remarks: additional info about the version of swift_rmvsy code
c Date: Feb 26th 2002

      subroutine util_version

      include '../swift.inc'
      include '../version.inc'

c  temporal variables
      integer ierr

c  functions
      integer length

      character*32 DRIVER
      common /drivername/ DRIVER
      data DRIVER /"swift_rmvsy"/       ! default driver name

c-----
c...  Executable code 

      write(*,990) DRIVER(1:length(DRIVER)),
     &  VERSION(1:length(VERSION))
990   format(
     &  '---------------------------------------------------',/,
     &  a,
     &  ' version ',a,/,
     &  'Author: Miroslav Broz, Charles University',/,
     &  '        e-mail: miroslav.broz@email.cz',/,
     &  '        http://sirrah.troja.mff.cuni.cz/yarko-site/',/,
     &  '---------------------------------------------------',/,
     &  ' ')

      write(*,1000) VER_NUM
 1000 format('************* SWIFT: Version ',f3.1,' *************')
      write(*,*) ' '
      write(*,*) 'Authors:'
      write(*,*) '   Martin Duncan: Queen''s University '
      write(*,*) '   Hal Levison: Southwest Research Institute '
      write(*,*) ' '
      write(*,*) ' Please address any comments or questions to:'
      write(*,*) '   Hal Levison '
      write(*,*) '   Geophysical, Astrophysical, & Planetary Sciences'
      write(*,*) '   Southwest Research Institute'
      write(*,*) '   1050 Walnut St.'
      write(*,*) '   Suite 429 '
      write(*,*) '   Boulder, Co  80302 '
      write(*,*) '   (303) 546-0290 '
      write(*,*) '   Fax: (303) 546-9687 '
      write(*,*) '   (D)  swri::levison '
      write(*,*) '   (I)  hal@gort.space.swri.edu '
      write(*,*) ' '

      return

      end  ! util_exit

c---------------------------------------------------


