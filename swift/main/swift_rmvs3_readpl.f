c**********************************************************************
c		      SWIFT_RMVS3_READPL.F
c**********************************************************************
c
c                 INCLUDES CLOSE ENCOUNTERS
c                 To run, need 3 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c                   test particle file like   tp.in
c
c Authors:  Hal Levison \& Martin Duncan
c Date:    8/25/94
c Last revision: 20/10/1998. Modified by Morby in order
c to read the planetesimals from an orbital file.
c look at read(*,77) for file format. Angles in radians.
c semimajor axis in AU and timescale in years!
c dtout and dtdump are automatically set equal to the
c timestep of the planetesimals datafile (which can be variable).
c code suitable for the case where the number of the planetesimals
c and their masses change with time. Mass of the sun is fixed forever.
c Needs to be linked with rmvs3_step_readpl and with the ususal swift library.
c Be sure that NPLMAX in swift.inc is large enough to allow the
c maximum number of planetesimals.  
     
	include '../swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
	
	real*8 xhend(NPLMAX),yhend(NPLMAX),zhend(NPLMAX)
	real*8 vxhend(NPLMAX),vyhend(NPLMAX),vzhend(NPLMAX)
	real*8 ta,tb

	integer istat(NTPMAX,NSTAT),i1st,i1sthill,n,nstep
	integer nbod,ntp,nleft,ip,ipt,ipl,ialpha,iplprec,ii
	integer iflgchk,iub,iuj,iud,iue,ipass
        real*8 rstat(NTPMAX,NSTATR)

	real*8 a0(NPLMAX),a1(NPLMAX),e0(NPLMAX),e1(NPLMAX)
	real*8 ri0(NPLMAX),ri1(NPLMAX),omega0(NPLMAX)
	real*8 omega1(NPLMAX),capom0(NPLMAX),capom1(NPLMAX)
	real*8 capm0(NPLMAX),capm1(NPLMAX)
	real*8 rms(NPLMAX)
	real*8 aa(NPLMAX),da(NPLMAX),ea(NPLMAX),de(NPLMAX)
	real*8 ria(NPLMAX),di(NPLMAX),oma(NPLMAX)
	real*8 dom(NPLMAX),coma(NPLMAX),dcom(NPLMAX)
	real*8 cma(NPLMAX),freq(NPLMAX)
	real*8 Cmfin,corr, eoff
	real*8 gm,a,e,inc,capom,omega,capm

	real*8 t0,tstop,dt,dtout,dtdump,dtorig
	real*8 t,tout,tdump,tfrac

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose 

	character*80 outfile,inparfile,inplfile,intpfile,fopenstat
        common /i1sthill_c/i1sthill


c-----
c...    Executable code 

c       initialise a0 to 0 because it will be the test to check whether
c       a planetesimal already existed at the beginning of the timestep
	do ip=1,NPLMAX
	   a0(ip)=0.d0
	enddo

c...    print version number
        call util_version

c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	call get_command_argument(1, inparfile)
	call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     $     iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)
	dtorig=dt

C Here, we assume that the J2 and J4 moments of the Sun are null
	j2rp2 = 0.d0
	j4rp4 = 0.d0

c Prompt and read name of planet data file
	write(*,*) ' '
	write(*,*) 'Enter name of planet data file : '
	call get_command_argument(2, inplfile)
999 	format(a)
	open(unit=77,file=inplfile,status='old')
	read(77,*)rms(1)
	mass(1)=rms(1)
 99	read(77,*)t,nbod
c	write(*,*)t,nbod
c       nbod is the number of planetesimals. Swift wants to add the sun
	nbod=nbod+1
	if (t.lt.t0) then
	   do ip=2,nbod
	      read(77,*)ipl,rms(ipl),a0(ipl),e0(ipl),ri0(ipl),
     $           omega0(ipl),capom0(ipl),capm0(ipl)
	      rplsq(ipl)=(0.5)**2
	   enddo
	   goto 99
	else
	   if(t.gt.t0)then
	      write(*,*)'time ',t0,' not found in planetary datafile'
	      stop
	   endif
	   if(nbod.gt.NPLMAX)then
	      write(*,*)'NPLMAX=',NPLMAX,'is not large enough'
	      stop
	   endif
c	write(*,*)'!!!! this code multiplies the masses by 10!!!!!'
c Initializes the a0 array
	   do ip = 1, NPLMAX
	      a0(ip) = 0.
	   end do
	   do ip=2,nbod
	      read(77,*)ipl,rms(ipl),a0(ipl),e0(ipl),ri0(ipl),
     $   	   omega0(ipl),capom0(ipl),capm0(ipl)
c             if(ipl.ne.91.and.ipl.ne.92)then
c             mass(ip)=rms(ipl)*100.
c	      else
	      mass(ip)=rms(ipl)
c	      endif
	      a=a0(ipl)
	      e=e0(ipl)
	      inc=ri0(ipl)
	      capom=capom0(ipl)
	      omega=omega0(ipl)
	      capm=capm0(ipl)
	      gm = mass(1)+mass(ip)
	      ialpha = -1
	      call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     $   	   xh(ip),yh(ip),zh(ip),
     $   	   vxh(ip),vyh(ip),vzh(ip))
	   enddo
	end if
	write(*,*)t,nbod
	ta=t0


c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	call get_command_argument(3, intpfile)
	call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)

c Initialize initial time and times for first output and first dump
	t = t0

        iub = 20
        iuj = 30
        iud = 40
        iue = 60

c...    Do the initial io write
        if(btest(iflgchk,0))  then ! bit 0 is set
           call io_write_frame(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &         xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
        endif
        if(btest(iflgchk,1))  then ! bit 1 is set
           call io_write_frame_r(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &         xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
        endif
        if(btest(iflgchk,2))  then    ! bit 2 is set
	   eoff = 0.d0
           call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &          vyh,vzh,iue,fopenstat, eoff)
        endif
        if(btest(iflgchk,3))  then    ! bit 3 is set
           call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
        endif

c...    must initize discard io routine
        if(btest(iflgchk,4))  then ! bit 4 is set
           call io_discard_write(0,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &          vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     &          'discard.out',fopenstat,nleft)
        endif

        nleft = ntp
c Sun's coordinates (once forever)
	   xh(1)=0.d0
	   yh(1)=0.d0
	   zh(1)=0.d0
	   vxh(1)=0.d0
	   vyh(1)=0.d0
	   vzh(1)=0.d0
	   xhend(1)=0.d0
	   yhend(1)=0.d0
	   zhend(1)=0.d0
	   vxhend(1)=0.d0
	   vyhend(1)=0.d0
	   vzhend(1)=0.d0

c***************here's the big loop *************************************
        write(*,*) ' ************** MAIN LOOP ****************** '

	ipass=0
	do while (t .le. tstop.and. (nleft.gt.0) )
	   ipass=ipass+1
	   read(77,*,end=111)tb,nbod
	   write(*,*)tb,nbod
c	   if(tb.gt.ta)then
c       add the Sun to the number of massive bodies
	   nbod=nbod+1
c       forces the output to follow the timestep in the planetesimals
c       datafile; forces dt such that there is an integer number of 
c       steps in dtout
	   dtout=tb-ta
	   tout = tb
	   nstep=dtout/dtorig
c	   if(dtout.eq.0..or.nstep.eq.0)write(*,*)'error:',dtout,nstep
	   dt=dtout/nstep
c	   write(*,*)'dt=',dt
	   n=1
	   ip=1
	   iplprec=1
	   i1st = 0
	   i1sthill = 0
	   do ipt = 1, NPLMAX
	      a1(ipt) = 0.
	   end do
	   do ipt=2,nbod
	      read(77,*)ipl,rms(ipl),a1(ipl),e1(ipl),ri1(ipl),
     $   	   omega1(ipl),capom1(ipl),capm1(ipl)
c       computes initial value (at t=ta) and time derivatives of the
c       planetesimals' osculating elements. It rejects planetesimals that 
c       don't survive at t=tb. Sets the masses equal to those at t=tb.
c       A planetesimal that is created at t=tb is not considered in this
c       timestep. It will be considered in the following one
	      if(a0(ipl).ne.0)then
		 ip=ip+1
c		 if(ipl.ne.91.and.ipl.ne.92)then
c		    mass(ip)=rms(ipl)*100.
c		 else
		 mass(ip)=rms(ipl)
c		 endif
		 aa(ip)=a0(ipl)
		 da(ip)=(a1(ipl)-a0(ipl))/dtout
		 ea(ip)=e0(ipl)
		 de(ip)=(e1(ipl)-e0(ipl))/dtout
		 ria(ip)=ri0(ipl)
		 di(ip)=(ri1(ipl)-ri0(ipl))/dtout
		 oma(ip)=omega0(ipl)
		 dom(ip)=omega1(ipl)-omega0(ipl)
		 if(dom(ip).gt.pi)dom(ip)=dom(ip)-2.d0*pi
		 if(dom(ip).lt.-pi)dom(ip)=dom(ip)+2.d0*pi
		 dom(ip)=dom(ip)/dtout
		 coma(ip)=capom0(ipl)
		 dcom(ip)=capom1(ipl)-capom0(ipl)
		 if(dcom(ip).gt.pi)dcom(ip)=dcom(ip)-2.d0*pi
		 if(dcom(ip).lt.-pi)dcom(ip)=dcom(ip)+2.d0*pi
		 dcom(ip)=dcom(ip)/dtout
		 cma(ip)=capm0(ipl)
c       guess the right mean motion frequency. Here a must be in AU and t in
c       years.
		 freq(ip)=2.d0*pi*sqrt(1.+mass(ip)/mass(1))
     $   	      /aa(ip)**1.5
		 Cmfin=cma(ip)+freq(ip)*dtout
		 Cmfin=mod(Cmfin,2.*pi)
		 corr=capm1(ipl)-Cmfin
		 if(corr.gt.pi)corr=corr-2.d0*pi
		 if(corr.lt.-pi)corr=corr+2.d0*pi
c       corrects the mean motion frequency in order to have the good final 
c       position of the planetesimal
		 freq(ip)=freq(ip)+corr/dtout
	      endif
	   enddo
	   do ipl = 1, NPLMAX
	      a0(ipl)=a1(ipl)
	      e0(ipl)=e1(ipl)
	      ri0(ipl)=ri1(ipl)
	      omega0(ipl)=omega1(ipl)
	      capom0(ipl)=capom1(ipl)
	      capm0(ipl)=capm1(ipl)
	   end do
	   nbod=ip
	   write (*,*) 'Actual number of planets: ', nbod-1
	   
c this loop looks useless during the integration but it is not
c since the number and order of planetesimals might have changed
c so that the planetesimal pointed by the index ip might have changed
	   do ip=2,nbod
	      a=aa(ip)
	      e=ea(ip)
	      inc=ria(ip)
	      capom=coma(ip)
	      omega=oma(ip)
	      capm=cma(ip)
	      gm = mass(1)+mass(ip)
	      ialpha = -1
	      call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     $   	   xh(ip),yh(ip),zh(ip),
     $   	   vxh(ip),vyh(ip),vzh(ip))
	   enddo
	   
	   do while ( (n .le. nstep))
	      
	      do ip=2,nbod
		 a=aa(ip)+da(ip)*(t+dt-ta)
		 e=ea(ip)+de(ip)*(t+dt-ta)
		 inc=ria(ip)+di(ip)*(t+dt-ta)
		 capom=coma(ip)+dcom(ip)*(t+dt-ta)
		 omega=oma(ip)+dom(ip)*(t+dt-ta)
		 capm=cma(ip)+freq(ip)*(t+dt-ta)
		 gm = mass(1)+mass(ip)
		 ialpha = -1
		 call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     $   	      xhend(ip),yhend(ip),zhend(ip),
     $   	      vxhend(ip),vyhend(ip),vzhend(ip))
	      enddo
	      
	      
	      call rmvs3_step_readpl(i1st,t,nbod,ntp,mass,j2rp2
     $   	   ,j4rp4,xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht
     $   	   ,vyht,vzht,istat,rstat,dt,xhend,
     $   	   yhend,zhend,vxhend,vyhend,vzhend)
	      
	      t = t + dt
	      n=n+1
c the new step subroutine does not update the positions of the planets!
	      do ip=2,nbod
		 xh(ip)=xhend(ip)
		 yh(ip)=yhend(ip)
		 zh(ip)=zhend(ip)
		 vxh(ip)=vxhend(ip)
		 vyh(ip)=vyhend(ip)
		 vzh(ip)=vzhend(ip)
	      enddo

	      if(btest(iflgchk,4))  then ! bit 4 is set
		 call discard(t,dt,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     $   	      xht,yht,zht,vxht,vyht,vzht,rmin,rmax,rmaxu,
     $   	      qmin,lclose,rplsq,istat,rstat)
		 call io_discard_write(1,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     $   	      vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     $   	      'discard.out',fopenstat,nleft)
	      else
		 nleft = ntp
	      endif
	   enddo
	   t=tb
	   
c       if it is time, output orb. elements, 
	   
	   if(ipass.eq.1)then
	      ipass=0
	      if(btest(iflgchk,0))  then ! bit 0 is set
		 call  io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,
     $   	      vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     $   	      iub,fopenstat)
	      endif
	      if(btest(iflgchk,1))  then ! bit 1 is set
		 call  io_write_frame_r(t,nbod,ntp,mass,xh,yh,zh,vxh,
     $   	      vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     $   	      iub,fopenstat)
	      endif
	      
	      
c       If it is time, do a dump
	      
	      tfrac = (t-t0)/(tstop-t0)
	      write(*,998) t,tfrac,nleft
 998	      format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     $   	   ': Number of active tp =',i4)
	      call io_dump_pl('dump_pl.dat',nbod,mass,xh,yh,zh,
     $   	   vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	      call io_dump_tp('dump_tp.dat',ntp,xht,yht,zht,
     $   	   vxht,vyht,vzht,istat,rstat)
	      call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     $   	   dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
	      
	      if(btest(iflgchk,2))  then ! bit 2 is set
		 call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,
     $   	      xh,yh,zh,vxh,vyh,vzh,iue,fopenstat, eoff)
	      endif
	      if(btest(iflgchk,3))  then ! bit 3 is set
		 call anal_jacobi_write(t,nbod,ntp,mass,xh,yh,zh,vxh,
     $   	      vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,
     $   	      iuj,fopenstat)
	      endif
	   endif
	   
	   ta=tb
c	   endif
	enddo   
c********** end of the big loop from time 't0' to time 'tstop'
	
c       Do a final dump for possible resumption later 
	
 111	call io_dump_pl('dump_pl.dat',nbod,mass,xh,yh,zh,
     $       vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	call io_dump_tp('dump_tp.dat',ntp,xht,yht,zht,
     $       vxht,vyht,vzht,istat,rstat)
	call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     $       dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
	
        call util_exit(0)
        end			! swift_rmvs3.f
c---------------------------------------------------------------------


