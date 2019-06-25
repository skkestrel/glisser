c   Gererate initial position and velocity for test particles
c   All of them have the same e and i.  The user supplies a range in a. 
c   The tp are equally spaced in a.  The rest of the angles are chosen 
c   at random.

      include 'swift.inc'

      real*8 SMASSYR
      parameter(SMASSYR=TWOPI*TWOPI)

      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

      real*8 gm,a,e,inc,capom,omega,capm,dr,amin,amax,da, ad, delta
      real*8 jmpran,rstat(NTPMAX,NSTATR),q, se, si, rayl

      integer istat(NTPMAX,NSTAT),itp
      integer ntp,ialpha,i,j,iseed,iuflg

      character*80 intpfile,datafile

      external jmpran, rayl

      data dr/1.7453292519943294d-2/                              

 1    write(*,'(a)') '# Units Menu: '
      write(*,'(a)') '#       0 ==> Solar masses, and AU '
      write(*,'(a)') '#       1 ==> AU, and Years '
      read(*,*) iuflg
      if( (iuflg.ne.0) .and. (iuflg.ne.1) ) goto 1

      write(*,'(a)') '# Input number particles '
      read(*,*) ntp
      if (ntp .gt. NTPMAX) then
         ntp = NTPMAX
         print *, '# restricting to', ntp, 'particles.'
      end if
c      ntp = 50
      amin = 5.0d0
      amax = 50.0d0
      write (6, '(a)') '# Input a_min and a_max: '
      read (5, *) amin, amax

c Estimate Delta
      ad = (amax**0.5d0 + amin**0.5d0)/2.d0
 100  continue
      delta = (amax - amin)/ad/dble(ntp)
      a = amin
      do i = 1, ntp
         a = a + delta*a**0.5d0
      end do
      if (dabs(a-amax) .gt. 0.01d0*(amax - amin)) then
         ad = (a - amin)/delta/dble(ntp)
         goto 100
      end if

      if(iuflg.eq.0) then
         gm = 1.0d0
      else
         gm = SMASSYR
      endif

      se = 0.002d0
      si = 0.001d0

      iseed=135671
      write (6, '(a)') '# Input random number generator seed: '
      read (5, *) iseed
      a = amin
      do i=1,ntp
         ialpha = -1
         a = a + delta*a**0.5d0
         e = rayl(iseed)*se
         inc = rayl(iseed)*si
         capm=jmpran(iseed)*2.*pi
         omega=jmpran(iseed)*2.*pi
         capom=jmpran(iseed)*2.*pi
          
         write(*,*)i,a,e,inc,capom,omega,capm

         call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     &        xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i))
         do j=1,NSTAT
            istat(i,j) = 0
         enddo
         do j=1,NSTATR
            rstat(i,j) = 0.0d0
         enddo
      enddo

      write(*,*) 'Enter name of test particle data file : '
      read(*,1001) intpfile
 1001 format(a)
      call io_dump_tp(intpfile,ntp,xht,yht,zht,
     &     vxht,vyht,vzht,istat,rstat)
      
      stop
      end    !  gen_a

c-------------------------------------------------------------------------

	real*8 function jmpran (i)
c********************************************************************
c
c cette subroutine calcule un nombre aleatoire par la formule de
c recurrence :
c i(k+1) = i(k) * 367379597 + 1 mod(2**31).
c x est un nombre aleatoire reel double precision compris entre 0.
c et 1. x = i/2.**31
c
c lors du premier appel de psalun, il faut donner une valeur
c d'initialisation a i puis ne plus changer i.
c
c********************************************************************
      implicit none
	integer*4
     1	  i, k, mask

	parameter
     1	  (k = 367379597)

	data
     1	  mask	/x'7fffffff'/

	i = i*k + 1
	i = iand(i,mask)
	jmpran = float(i)/2147483648.

	return
	end

      real*8 function rayl(i)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns a random number following a Rayleigh distribution
c P(x) dx \propto x exp(-x^2/2)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J.-M. Petit  Observatoire de Besancon
c Version 1 : March 2018
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     i     : random number generator seed (I4)
c
c OUTPUT
c     rayl  : random number following Rayleigh distribution (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      implicit none

      integer*4 i
      real*8 jmpran

      external jmpran

      rayl = dsqrt(-2.d0*dlog(jmpran(i)))
      return
      end
