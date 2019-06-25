c This program reads in the binary output of a SWIFT integration and
c convert it into ASCII in a way that is useable by swift_rmvs3_readpl.
c Output is as follow:
c
c t, nob
c i, mass, a, e, inc, omega, capom, capM
c
c The last line is repeated nbod times.
c
c i is actually the absolute vale of the identifier of the planet in
c the binary file.

	include 'swift.inc'

	real*8 mass(NPLMAX),dr,peri
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT)
        real*8 rstat(NTPMAX,NSTATR)
	integer nbod,ntp,ierr,ifol,istep, nbod_init
	integer iflgchk,iu,nleft,i,id, j
        integer io_read_hdr,io_read_line
        integer io_read_hdr_r,io_read_line_r

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tmax(NPLMAX), tinit, tend, t1

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX), a_min(NPLMAX)
        logical*2 lclose
        real*8 a,e,inc,capom,omega,capm,j2rp2,j4rp4, tabr(7, NPLMAX)
        real*8 elh,elk,elp,elq,apo, aend(NPLMAX), a_max(NPLMAX)
	real*8 xcm, ycm, zcm, vxcm, vycm, vzcm, mtot, gm
	integer ialpha, iuo, nprint, niter, tabi(NPLMAX), npl

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
	nbod_init = nbod

	write (*,*) ' '
	write (6, '(a,$)') 'Start time: '
	read (5, *) tinit
	write (6, '(a,$)') 'End time:   '
	read (5, *) tend
	write (6, '(a,$)') 'Time step:  '
	read (5, *) dt

        iu = 20

        dr = 180.0/PI

	do i = 2, nbod_init
	   a_max(i) = 0.
	   a_min(i) = 100.
	end do

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
        open(unit=7,file='readpl.out')
	write (7, *) mass(1)

        write(*,*) '1 2 3  4    5     6    7    8    9 '
        write(*,*) 't,a,e,inc,capom,omega,capm,peri,apo'

	t1 = tinit
 1      continue
             if(btest(iflgchk,0))  then ! bit 0 is set
                ierr = io_read_hdr(iu,t,nbod,nleft) 
             else
                ierr = io_read_hdr_r(iu,t,nbod,nleft) 
             endif

             if(ierr.ne.0) goto 2

	     if (t .lt. t1-0.001*dt) then
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

	     npl = 0

             do i=2,nbod
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
                if(ierr.ne.0) goto 2
		if (e .lt. 1.d0) then
		   npl = npl + 1
		   tabi(npl) = -id
		   tabr(1, npl) = mass(-id)
		   tabr(2, npl) = a
		   tabr(3, npl) = e
		   tabr(4, npl) = inc
		   tabr(5, npl) = omega
		   tabr(6, npl) = capom
		   tabr(7, npl) = capm
		   tmax(-id) = t
		   aend(-id) = a
		   if (a > a_max(-id)) a_max(-id) = a
		   if (a < a_min(-id)) a_min(-id) = a
		end if
             enddo

	     write (7, 1001) t, npl
	     do i = 1, npl
		write (7, 1000) tabi(i), (tabr(j,i), j=1,7)
	     end do

             do i=1,nleft
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
             enddo

	     t1 = t + dt
	     if (t1 .gt. tend+0.1*dt) goto 2
        goto 1

 2      continue
	close (iu)
	close (7)

	open (unit=7, file='lastpl.out', status='unknown')
	write (7, '(1x, i4)') nbod_init
	do i = 2, nbod_init
	   write (7, 1002) i, tmax(i), aend(i), a_max(i), a_min(i)
	end do
	close (7)

 1000	format (1x,i4,1xe15.7,1x,1pg13.7,1x,0pf12.10,4(1x,f9.6))
 1001	format (1xe15.7,1x,i4)
 1002	format (1x, i4, 1x, f12.1, 3(1x, g12.4))

        stop
        end
