c converts swift output binary file to ascii file.
c Version suited for readpl_bin, i.e. get info from readpl binary format
c planet file.

	include 'swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),dr,peri
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT), idpl(NPLMAX)
        real*8 rstat(NTPMAX,NSTATR)
	integer nbod,ntp,ierr,ifol,istep
	integer iflgchk,iu,nleft,i,id
        integer io_read_hdr,io_read_line
        integer io_read_hdr_r,io_read_line_r

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tmax, tinit, tend, t1

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose
        real*8 a,e,inc,capom,omega,capm,j2rp2,j4rp4
        real*8 elh,elk,elp,elq,apo
	real*8 xcm, ycm, zcm, vxcm, vycm, vzcm, mtot, gm
	integer ialpha, iuo, nprint, niter

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
	call io_init_readpl(inplfile, nbod, mass, rplsq)

c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
	call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)

	write (*,*) ' '
	write (6, '(a,$)') 'Start time: '
	read (5, *) tinit
	write (6, '(a,$)') 'End time:   '
	read (5, *) tend
	write (6, '(a,$)') 'Time step:  '
	read (5, *) dt

        iu = 20

        dr = 180.0/PI

        if(btest(iflgchk,0)) then
           write(*,*) ' Reading an integer*2 binary file '
        else if(btest(iflgchk,1)) then
           write(*,*) ' Reading a real*4 binary file '
        else
           write(*,*) ' ERROR: no binary file format specified '
           write(*,*) '        in param file '
           stop
        endif

        open(unit=iu, file=outfile, status='old',form='unformatted')
        open(unit=7,file='slice_b_all.out')

        write(7,'(a)') '# id         t            a         e   '
	1 //'    inc      capom     omega     capm   '
	1 //'    peri       apo'

        tmax = t0
	t1 = tinit
 1      continue
             if(btest(iflgchk,0))  then ! bit 0 is set
                ierr = io_read_hdr(iu,t,nbod,nleft) 
             else
                ierr = io_read_hdr_r(iu,t,nbod,nleft) 
             endif

             if(ierr.ne.0) goto 2

	     if (t .lt. t1-0.000001d0*dt) then
		do i=2,nbod+nleft
		   if(btest(iflgchk,0))  then ! bit 0 is set
		      ierr = io_read_line(iu,id,a,e,inc,capom,omega
     &			,capm) 
		   else
		      ierr = io_read_line_r(iu,id,a,e,inc,capom,omega
     &			,capm) 
		   endif
		   if(ierr.ne.0) goto 2
		end do
		goto 1
	     end if

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
		idpl(i) = id
		mtot = mtot + mass(-id)
		gm = mass(1) + mass(-id)
c		print *, mtot, mass(-id), gm
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

             do i=2,nbod
		xh(i) = xh(i) - xcm
		yh(i) = yh(i) - ycm
		zh(i) = zh(i) - zcm
		vxh(i) = vxh(i) - vxcm
		vyh(i) = vyh(i) - vycm
		vzh(i) = vzh(i) - vzcm
		ialpha = -1
		call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),
     &		  vyh(i),vzh(i),mtot,ialpha,a,e,inc,capom,omega,
     &		  capm)
		inc = inc*dr
		capom = capom*dr
		omega = omega*dr
		capm = capm*dr
		peri = a*(1.0d0-e)
		apo = a*(1.0d0+e)
		write(7,1000) idpl(i), t,a,e,inc,capom,omega,capm,peri,
	1	  apo
             enddo

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
		call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),
     &		  vyht(i),vzht(i),mtot,ialpha,a,e,inc,capom,omega,
     &		  capm)
		elh = e*cos(omega+capom)
		elk = e*sin(omega+capom)
		elp = sin(inc/2.0)*cos(capom)
		elq = sin(inc/2.0)*sin(capom)
		inc = inc*dr
		capom = capom*dr
		omega = omega*dr
		capm = capm*dr
		peri = a*(1.0d0-e)
		apo = a*(1.0d0+e)
		write(7,1000) id, t,a,e,inc,capom,omega,capm,peri,apo
		tmax = t
             enddo

	     t1 = t1 + dt
	     if (t1 .gt. tend+0.1*dt) goto 2
        goto 1

 2      continue

        write(*,*) ' Tmax = ',tmax

 1000	format(1x,i4,1x,e15.7,1x,f10.4,1x,f7.5,4(1x,f9.4),2(1x,f10.4))

        stop
        end

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c                         IO_INIT_READPL
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c IO_INIT_READPL reads in the data for the Sun and planets.
c Simplified version that reads only the masses of Sun and planets.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c INPUT
c     infile: File name to read from (character*80)
c
c OUTPUT
c     nbod  : number of massive bodies (I4)
c     mass  : mass of bodies (NPLMAX*R8)
c     rplsq : min distance^2 that a tp can get from pl (NPLMAX*R8)
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      subroutine io_init_readpl(infile, nbod, mass, rplsq)

      include 'swift.inc'

      character*(*) infile
      real*8 mass(NPLMAX), rplsq(NPLMAX), t
      integer*4 nbod, rn, ip, id

      open(unit=77,file=infile,status='old',err=100,
     $  form='unformatted',access='direct',recl=60)
      rn = 1
      read (77, rec=rn) mass(1)
      rplsq(1) = 0.0d0
      rn = rn + 1
      read (77, rec=rn, err=111) t, nbod
      nbod = nbod + 1
      do ip = 2, nbod
         rn = rn + 1
         read (77, rec=rn) id, mass(ip)
         rplsq(ip)=(2.2e-7)**2
      end do
      close(unit=77)

      return

 100  continue
      print *, 'ERROR while opening', infile, '. Cowardly aborting.'
      stop

 111  continue
      print *, 'ERROR while reading', infile, '. Cowardly aborting.'
      stop

      end
