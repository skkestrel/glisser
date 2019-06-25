C This progam converts from heliocentric osculating elements to
C heliocentric cartisian coordinates to barycentric cartesian
C coordinates to barycentric osculating elements and writes back these
C elements.
C
C This is done only for the test particles. We keep the heliocentric
C coordinates for the pla,ets.

	include 'swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),dr,peri
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT)
        real*8 rstat(NTPMAX,NSTATR)
	integer nbod,ntp,ierr,ifol
	integer iflgchk,iu,nleft,i,id
        integer io_read_hdr,io_read_line
        integer io_read_hdr_r,io_read_line_r

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tmax

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose
        real*8 a,e,inc,capom,omega,capm,j2rp2,j4rp4
        real*8 elh,elk,elp,elq,apo
	real*8 xcm, ycm, zcm, vxcm, vycm, vzcm, mtot, gm
	integer ialpha, iuo
	real*4 ttmp
	integer*2 nleft2,nbod2

	character*80 outfile,inparfile,inplfile,intpfile,fopenstat

c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	read(*,999) inparfile
	call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
	write(*,*) ' '
	write(*,*) 'Enter name of planet data file : '
	read(*,999) inplfile
999 	format(a)
	call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &       vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
	call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)

        iu = 20
	iuo = 21

        dr = 180.0/PI

        if(btest(iflgchk,0)) then
           write(*,*) ' Reading an integer*2 binary file '
        else if(btest(iflgchk,1)) then
           write(*,*) ' Reading an real*4 binary file '
        else
           write(*,*) ' ERROR: no binary file format specified '
           write(*,*) '        in param file '
           stop
        endif

        open(unit=iu, file=outfile, status='old',form='unformatted')
	call io_open(iuo,'tp_bary.dat','unknown','UNFORMATTED',ierr)
	if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in io_write_frame: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
	endif

        tmax = t0
 1      continue
             if(btest(iflgchk,0))  then ! bit 0 is set
                ierr = io_read_hdr(iu,t,nbod,nleft) 
             else
                ierr = io_read_hdr_r(iu,t,nbod,nleft) 
             endif

             if(ierr.ne.0) goto 2
	     ttmp = t
	     nbod2 = nbod
	     nleft2 = nleft
	     write(iuo) ttmp,nbod2,nleft2

	     mtot = mass(1)
	     xcm = 0.d0
	     ycm = 0.d0
	     zcm = 0.d0
	     vxcm = 0.d0
	     vycm = 0.d0
	     vzcm = 0.d0
             do i=2,nbod
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
                if(ierr.ne.0) goto 2
		call io_write_line_r(iuo,id,a,e,inc,capom,omega,capm)
		mtot = mtot + mass(-id)
		gm = mass(1) + mass(-id)
		ialpha = -1
		call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     &		  xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i))
		xcm = xcm + mass(-id)*xh(i)
		ycm = ycm + mass(-id)*yh(i)
		zcm = zcm + mass(-id)*zh(i)
		vxcm = vxcm + mass(-id)*vxh(i)
		vycm = vycm + mass(-id)*vyh(i)
		vzcm = vzcm + mass(-id)*vzh(i)
             enddo
	     xcm = xcm/mtot
	     ycm = ycm/mtot
	     zcm = zcm/mtot
	     vxcm = vxcm/mtot
	     vycm = vycm/mtot
	     vzcm = vzcm/mtot
c	     write (6, *) xcm, ycm, zcm, vxcm, vycm, vzcm

	     gm = mass(1)
             do i=1,nleft
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
                if(ierr.ne.0) goto 2
		ialpha = -1
		call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     &		  xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i))
		xht(i) = xht(i) - xcm
		yht(i) = yht(i) - ycm
		zht(i) = zht(i) - zcm
		vxht(i) = vxht(i) - vxcm
		vyht(i) = vyht(i) - vycm
		vzht(i) = vzht(i) - vzcm
		ialpha = -1
		call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &		  vzht(i),mtot,ialpha,a,e,inc,capom,omega,capm)
c		write (6, *) id
c		write (6, *) xht(i), yht(i), zht(i)
c		write (6, *) a, e, inc
		call io_write_line_r(iuo,id,a,e,inc,capom,omega,capm)
             enddo

        goto 1

 2      continue

        write(*,*) ' Tmax = ',tmax

        stop
        end
