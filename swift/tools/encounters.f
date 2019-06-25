c**********************************************************************
c		      ENCOUNTERS.F
c**********************************************************************
c

     
	include '../11/swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX)
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT)
	integer nbod,ntp,i,j
	integer iflgchk,iact

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 tfrac,rstat(NTPMAX,NSTATR)

	real*8 rmin,rmax,rmaxu,qmin,gm,rplsq(NPLMAX),j2rp2,j4rp4
        real*8 a,e,inc,capom,omega,capm,peri,apo
        integer ialpha,nwhy(-4:NPLMAX)
        logical*2 lclose

	character*80 outfile,fopenstat

c Get parameters 
	call io_init_param('dump_param.dat',t0,tstop,dt,dtout,
     &       dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,
     &       outfile,fopenstat)

c Get data for sun and planetary orbits.
	call io_init_pl('dump_pl.dat',lclose,iflgchk,nbod,mass,xh,yh,
     &                  zh,vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
	call io_init_tp('dump_tp.dat',ntp,xht,yht,zht,
     &            vxht,vyht,vzht,istat,rstat)


c 
       open(unit=7,file='encounters.out')

       do i=-4,nbod
          nwhy(i) = 0
       enddo

       tfrac = t0/tstop
       if(t0.lt.10000.0) then
          write(*,*) ' Time = ',t0
          write(7,*) ' Time = ',t0
       else
          write(*,10) t0
          write(7,10) t0
 10       format(' Time = ',1p1e13.5)
       endif
       write(*,*) ' Fraction done = ',tfrac
       write(7,*) ' Fraction done = ',tfrac

       iact = 0
       do i=1,ntp
          if(istat(i,1).eq.0) then
             iact = iact + 1
          endif
       enddo
       write(*,*) iact,' out of ',ntp,' still active '
       write(7,*) iact,' out of ',ntp,' still active '

       write(7,*) ' '
       write(7,*) ' Planet : '
       write(7,*) '       #          a                e             i '
       do i=2,nbod
          gm = mass(1) + mass(i)
          call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),vyh(i),
     &          vzh(i),gm,ialpha,a,e,inc,capom,omega,capm)
          inc = inc*180.0d0/PI
          write(7,1000) i,a,e,inc
 1000     format(5x,i4,3(5x,f10.4))
       enddo

       write(7,*) ' '
       write(7,*) ' Particles with close encounters '
       write(7,1011) 
 1011  format('#Part stat fate last pl01 pl02 pl03 pl04 ...')

       iact = 0
       do i=1,ntp
	  lclose = .false.
	  do j = 2, nbod
	     if (istat(i,j+2) .gt. 0) lclose = .true.
	  end do
	  if (lclose) then
	     iact = iact + 1
	     write (7, 1021) i, (istat(i,j), j=1,nbod+2)
	  end if
       enddo
       if (iact .eq. 0) then
	  write (7, '(''NO CLOSE ENCOUNTER FOR ANY PARTICLE'')')
       end if
       write (7, *)
 1021  format (i5,13(1x,i4))

       stop
       end                      ! encounters.f
c---------------------------------------------------------------------






