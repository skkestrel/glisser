c*************************************************************************
c                            UTIL_HILLS.F
c*************************************************************************
c This subroutine calculates the hill's sphere for the planets using heliocentric distance
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c             Output:
c                  r2hill       ==>  the SQUARE of the planet's hill's sphere 
c                                    (real array)
c
c
c Remarks: 
c Authors:  Hal Levison , KZ
c Date:    2/19/93
c Last revision: 1/6/97

      subroutine util_hills2(nbod,mass,xh,yh,zh,vxh,vyh,vzh,r2hill) 

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Outputs
      real*8 r2hill(nbod)

c...  Internals
      integer i
      real*8 mu,energy,ap,rhil,r,v2

c-----
c...  Executable code 

      do i=2,nbod
         if(mass(i).ne.0.0d0) then
            mu = mass(1)*mass(i)/(mass(1)+mass(i))
            r2hill(i) = xh(i)*xh(i) + yh(i)*yh(i) + zh(i)*zh(i)
         else
            r2hill(i) = 0.0d0
         endif
      enddo
      
      r2hill(1) = 0.0
      
      return
      end                       ! util_hills2

c---------------------------------------------------
