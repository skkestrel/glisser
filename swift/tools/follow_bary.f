c converts binary file to ascii file

	include 'swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),dr,peri
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT)
        real*8 rstat(NTPMAX,NSTATR)
	integer nbod,ntp,ierr,ifol,istep
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
	call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &       vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
	call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)

        iu = 20

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

        write(*,*) ' Input the particle number to follow (no planet) '
        read(*,*) ifol
        write(*,*) ' Following particle ',ifol

	write (*, *) ' Input the print frequency '
	read (*, *) nprint

        open(unit=iu, file=outfile, status='old',form='unformatted')
        open(unit=7,file='follow_b.out')

        write(*,*) '1 2 3  4    5     6    7    8    9 '
        write(*,*) 't,a,e,inc,capom,omega,capm,peri,apo'

	niter = 0
        tmax = t0
 1      continue
             if(btest(iflgchk,0))  then ! bit 0 is set
                ierr = io_read_hdr(iu,t,nbod,nleft) 
             else
                ierr = io_read_hdr_r(iu,t,nbod,nleft) 
             endif

             if(ierr.ne.0) goto 2

	     niter = niter + 1

	     mtot = mass(1)
	     xcm = 0.d0
	     ycm = 0.d0
	     zcm = 0.d0
	     vxcm = 0.d0
	     vycm = 0.d0
	     vzcm = 0.d0
             istep = 0
             do i=2,nbod
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
                if(ierr.ne.0) goto 2
		mtot = mtot + mass(-id)
		gm = mass(1) + mass(-id)
		if (dabs(e-1.d0) .lt. TINY) then
		   ialpha = 0
		else if (e .lt. 1.d0) then
		   ialpha = -1
		else
		   ialpha = 1
		end if
		call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     &		  xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i))
		xcm = xcm + mass(-id)*xh(i)
		ycm = ycm + mass(-id)*yh(i)
		zcm = zcm + mass(-id)*zh(i)
		vxcm = vxcm + mass(-id)*vxh(i)
		vycm = vycm + mass(-id)*vyh(i)
		vzcm = vzcm + mass(-id)*vzh(i)
                if (id.eq.ifol) then
                   istep = 1
		   if  (mod(niter,nprint) .eq. 0) then
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
		      if (e .lt. 1.d0) then
			 write(7,1000) t,a,e,inc,capom,omega,capm,peri
     &			   ,apo
		      else
			 peri = dsqrt(xh(i)**2 + yh(i)**2 + zh(i)**2)
			 write(7,1002) t,a,e,inc,capom,omega,capm,peri
		      end if
		      tmax = t
		   end if
                endif
             enddo
	     xcm = xcm/mtot
	     ycm = ycm/mtot
	     zcm = zcm/mtot
	     vxcm = vxcm/mtot
	     vycm = vycm/mtot
	     vzcm = vzcm/mtot

             do i=1,nleft
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
                if(ierr.ne.0) goto 2
                if (id.eq.ifol) then
                   istep = 1
		   if  (mod(niter,nprint) .eq. 0) then
		      if (dabs(e-1.d0) .lt. TINY) then
			 ialpha = 0
		      else if (e .lt. 1.d0) then
			 ialpha = -1
		      else
			 ialpha = 1
		      end if
		      call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     &			   xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i))
		      xht(i) = xht(i) - xcm
		      yht(i) = yht(i) - ycm
		      zht(i) = zht(i) - zcm
		      vxht(i) = vxht(i) - vxcm
		      vyht(i) = vyht(i) - vycm
		      vzht(i) = vzht(i) - vzcm
		      call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),
     &			   vyht(i),vzht(i),mtot,ialpha,a,e,inc,capom,omega,
     &			   capm)
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
		      if (e .lt. 1.d0) then
			 write(7,1000) t,a,e,inc,capom,omega,capm,peri
     &			   ,apo
		      else
			 peri = dsqrt(xht(i)**2 + yht(i)**2 + zht(i)**2)
			 write(7,1002) t,a,e,inc,capom,omega,capm,peri
		      end if
		      tmax = t
		   end if
                endif
             enddo
             if(istep.eq.0) goto 2     ! did not find particle this times step

        goto 1

 2      continue

        write(*,*) ' Tmax = ',tmax

 1000	format(1x,e15.7,1x,f10.4,1x,f7.5,4(1x,f9.4),2(1x,f10.4))
 1002	format(1x,e15.7,1x,f10.4,1x,g9.3,3(1x,f9.4),2(1x,g11.4))

        stop
        end
