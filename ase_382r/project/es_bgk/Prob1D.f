	program Prob1D
	data pi/3.1415926536/
	real phi(500,-10:10,-10:10,-10:10), phif(-10:10,-10:10,-10:10)
	real phi1(500,-10:10,-10:10,-10:10)	
	real nd(500), ux(500), T(500), vy(500), wz(500)
	real ndm(500,11), uxm(500,11), vym(500,11), Tm(500,11)
	real pm(500,11), qm(500,11), tauxxm(500,11)
	real p(500), q(500), tauxx(500), norm, ndf, uf, Tf, timestamp(11)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	common /coll/ phif
c
c	 phi - scaled density normalized distribution function
c	 nd - scaled number density
c	 ux - scaled x gas velocity, similarly vy, wz
c	 T  - scaled temperature
c 	 p  - scaled pressure
c	 q  - scaled heat flux - not yet implemented
c	 tauxx - scaled normal shear stress - not yet implemented
c	 ndm, uxm, Tm are matrices to get time-dependent spatial profiles of variables
c        Need to initialize arrays tauxxm and qxm for time-dependent profiles of tauxx and q
c
	open(unit=11, file='1Dinp.txt',status='old')
	open(unit=12, file='phistart.txt',status='unknown')
	open(unit=13, file='phimid.txt',status='unknown')
	open(unit=14, file='phiend.txt',status='unknown')
	open(unit=15, file='Dprof.txt',status='unknown')
	open(unit=16, file='Uprof.txt',status='unknown')
	open(unit=17, file='Tprof.txt',status='unknown')
	open(unit=18, file='Pprof.txt',status='unknown')
	open(unit=19, file='Qprof.txt',status='unknown')
	open(unit=20, file='TAUprof.txt',status='unknown')
c       
c	Read space, velocity, and time scaled discretization steps
c
c 	The variables are as follows:
c
c 	nspace 		- number of space steps
c 	iv*max/min 	- max/min normalized velocity space index
c 	alphax 		- scaled deltaX
c 	betav 		- scaled deltaV
c 	deltat 		- scaled deltaT
c
c 	ndf			- scaled freestream density
c 	uf 			- scaled freestream velocity
c 	Tf 			- scaled freestream Temperature
c
c 	ntstep 		- number of time steps total
c 	nsplot 		- spatial index location to use for dist. output
c 	npr 		- print interval (in units of time steps)
c
	read(11,*) nspace, ivxmin, ivxmax, ivymin, ivymax, ivzmin, ivzmax
	if (nspace.gt.500) then
	  write(*,9003) nspace
 9003	  format('# space points ',i3,'> max value 500. Reset!')
	  stop
	endif
	read(11,*) alphax, betav, deltat
c
c	 Read scaled freestreeam velocity
c
	read(11,*) ndf, uf, Tf
	write(*,3900) alphax, betav, deltat, ndf, uf, Tf
	read(11,*) ntstep, nsplot, npr
	write(*,9000) npr
 9000	format('Print interval = ',i3)
!	if (ntstep.gt.10) then
!	npr=int(ntstep/10)	! Compute profile printing interval
!	else
!	  npr=1
!	endif
	ipcount=int(ntstep/npr) ! Check if print interval will lead to print array bounds being exceeded
	if(ipcount.gt.11) then
	   npr=int(ntstep/10)
	   write(*,9001)npr
 9001	   format('Print interval recomputed. New value = ', i3)
	endif
	ipcount=npr ! Initialize counter for print interval to print initial values
	ipindx=1 ! Initialize print matrix column index
c 
c       Initialize distribution function. Freestream density and temperature are typically
c	reference values (unity), but to preserve generality they are free inputs.
c
! Initialize free-stream distribution - shifted Maxwellian 
	call Maxwell_1(betav,ndf,uf,Tf,phif) ! phif is free-stream phi
! Set up flow field
	do ns=1,nspace
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax
		phi(ns,i,j,k)=phif(i,j,k)
	      enddo
	    enddo
	  enddo
	enddo
	
! Write out initial 2-D slice of distribution function at selected physical space point
! (1-D problem in space so phi(i,j,0)=phi(i,0,k))
	do i=ivxmin,ivxmax
	  do j=ivymin,ivymax
	      write(12,2000)i,j,phi(nsplot,i,j,0)
 2000	      format(2i4,2x,1pe15.6)
	  enddo
	enddo
! Initialization complete. Begin solution of 1D BGK equation.
! Solution via time splitting;  collide then convect.
	do ntime=1,ntstep
! Compute properties for Krook collision operator	  
	  call density(betav,nd)
	  call xvelocity(betav,nd,ux)
	  call yvelocity(betav,nd,vy)
	  call zvelocity(betav,nd,wz)	
	  call Temp(betav,nd,ux,T)
	  call pressure(betav,nd,T,p)
	  call heatflux(betav,nd,ux,q)
	  call tau11(betav,nd,ux,p,tauxx)
! write properties into property matrices for time dependent profiles	  
	  if (ipcount.eq.npr) then
	    do ns=1,nspace
	      ndm(ns,ipindx)=nd(ns)
	      uxm(ns,ipindx)=ux(ns)
c	      vym(ns,ipindx)=vy(ns)
	      Tm(ns,ipindx)=T(ns)
	      pm(ns,ipindx)=p(ns)
	      qm(ns,ipindx)=q(ns)
	      tauxxm(ns,ipindx)=tauxx(ns)
	    enddo
	      write(*,9002)ipindx,ntime
 9002	      format('ipindx =',i3,5x,'ntime =',i3)
	    ipcount=1 ! Reset print interval counter after printing
	    ipindx=ipindx+1 ! Increment print column index
	  else
	    ipcount=ipcount+1 ! Increment print interval counter
	  endif
	  ! Ellipsoidal Statistical collision term (ES-BGK)
	  call ESBGKcoll_1(nd, ux, vy, wz, T, betav, deltat)
	  ! Krook collision term (BGK)
	  ! call Krookcoll_1(nd,ux,vy,wz,T,betav,deltat)
	  if(ntime.eq.(ntstep/2)) then
	    ! Write out intermediate distribution function
	    do i=ivxmin,ivxmax
	      do j=ivymin,ivymax
		write(13,2000)i,j,phi(nsplot,i,j,0)
	      enddo
	    enddo
	  endif
	  call convect(alphax,betav,deltat)
	enddo
! Write out final distribution function at point chosen for plot
	do i=ivxmin,ivxmax
	  do j=ivymin,ivymax
	    write(14,2000)i,j,phi(nsplot,i,j,0)
	  enddo
	enddo
! Compute and property profiles at end of computation
	  call density(betav,nd)
	  call xvelocity(betav,nd,ux)
	  call yvelocity(betav,nd,vy)
	  call zvelocity(betav,nd,wz)	
	  call Temp(betav,nd,ux,T)
	  call pressure(betav,nd,T,p)
	  call heatflux(betav,nd,ux,q)
	  call tau11(betav,nd,ux,p,tauxx)
! write properties into property matrices for time dependent profiles	  
	  do ns=1,nspace
	    ndm(ns,ipindx)=nd(ns)
	    uxm(ns,ipindx)=ux(ns)
c	    vym(ns,ipindx)=vy(ns)
	    Tm(ns,ipindx)=T(ns)
	    pm(ns,ipindx)=p(ns)
	    qm(ns,ipindx)=q(ns)
	    tauxxm(ns,ipindx)=tauxx(ns)
	  enddo	
! Write out property profiles
	do ipr=1,ipindx-1
	  timestamp(ipr)=float(ipr-1)*npr*deltat
	enddo
	timestamp(ipindx)=float(ntstep)*deltat
	write(15,3900) alphax, betav, deltat, ndf, uf, Tf
3900	format('alpha(x) = ',f6.3,2x,'beta(v) = ',f6.3,2x,
     c  'delta (t) =', f6.3/'Free stream: density = ',f6.3,2x,
     c  'velocity = ',f6.3,2x,'temperature = ',f6.3)
	write(16,3900) alphax, betav, deltat, ndf, uf, Tf
	write(17,3900) alphax, betav, deltat, ndf, uf, Tf
	write(18,3900) alphax, betav, deltat, ndf, uf, Tf
	write(19,3900) alphax, betav, deltat, ndf, uf, Tf
	write(20,3900) alphax, betav, deltat, ndf, uf, Tf
	write(15,4001)(timestamp(ipr),ipr=1,ipindx)
	write(16,4001)(timestamp(ipr),ipr=1,ipindx)
	write(17,4001)(timestamp(ipr),ipr=1,ipindx)
	write(18,4001)(timestamp(ipr),ipr=1,ipindx)
	write(19,4001)(timestamp(ipr),ipr=1,ipindx)
	write(20,4001)(timestamp(ipr),ipr=1,ipindx)
 4001	format('x',',',11('t= ',f8.3,",",2x))
	do ns=1,nspace
	   xx=float(ns-1)*alphax
	  write(15,4002)xx,(ndm(ns,ipr),ipr=1,ipindx)
	  write(16,4002)xx,(uxm(ns,ipr),ipr=1,ipindx)
	  write(17,4002)xx,(Tm(ns,ipr),ipr=1,ipindx)
	  write(18,4002)xx,(pm(ns,ipr),ipr=1,ipindx)
	  write(19,4002)xx,(qm(ns,ipr),ipr=1,ipindx)
	  write(20,4002)xx,(tauxxm(ns,ipr),ipr=1,ipindx)
 4002	  format(12(f12.6,',',3x))
	enddo
	stop
	end program

ccccccccccccccccccccc density cccccccccccccccccccccccccccccc
	subroutine density(betav,dens)
	real phi(500,-10:10,-10:10,-10:10),dens(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	betacub=betav**3
	do ns=1,nspace
	  dens(ns)=0.
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	     do k=ivzmin,ivzmax
			dens(ns)=dens(ns)+phi(ns,i,j,k)
	     enddo
	    enddo
	  enddo
	  dens(ns)=dens(ns)*betacub
	enddo
	return
	end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccc xvelocity ccccccccccccccccccccccccccc
	subroutine xvelocity(betav,dens,xvel)
	real phi(500,-10:10,-10:10,-10:10),dens(500),xvel(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	betacub=betav**3	
	do ns=1,nspace
	  xvel(ns)=0.
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax
			xvel(ns)=xvel(ns)+phi(ns,i,j,k)*betav*float(i)
	      enddo
	    enddo
	  enddo
	  xvel(ns)=xvel(ns)*betacub/dens(ns)
	enddo
	return
	end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccc Temperature ccccccccccccccccccccccccc
	subroutine Temp(betav,dens,xvel,Tcal)
! This routine has been modified to account for the fact that eta_ref=sqrt(k*T_ref/m) NOT  sqrt(2*k*T_ref/m)
	real phi(500,-10:10,-10:10,-10:10),dens(500),xvel(500),Tcal(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	betasq=betav*betav
	betacub=betasq*betav
	do ns=1,nspace
	  Tcal(ns)=0.
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax
	        Cxsq=(betav*float(i)-xvel(ns))**2
			Cysq=betasq*float(j*j)
			Czsq=betasq*float(k*k)
			Csq=Cxsq+Cysq+Czsq
			Tcal(ns)=Tcal(ns)+phi(ns,i,j,k)*Csq
	      enddo
	    enddo
	  enddo
	  Tcal(ns)=Tcal(ns)*betacub/(3.*dens(ns))
	enddo
	return
	end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccc yvelocity cccccccccccccccccccccccc
	subroutine yvelocity(betav,dens,yvel)
	real phi(500,-10:10,-10:10,-10:10),dens(500),yvel(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	betacub=betav**3	
	do ns=1,nspace
	  yvel(ns)=0.
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax
			yvel(ns)=yvel(ns)+phi(ns,i,j,k)*betav*float(j)
	      enddo
	    enddo
	  enddo
	  yvel(ns)=yvel(ns)*betacub/dens(ns)
	enddo
	return
	end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccc zvelocity ccccccccccccccccccccccccc
	subroutine zvelocity(betav,dens,zvel)
	real phi(500,-10:10,-10:10,-10:10),dens(500),zvel(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	betacub=betav**3	
	do ns=1,nspace
	  zvel(ns)=0.
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax
			zvel(ns)=zvel(ns)+phi(ns,i,j,k)*betav*float(k)
	      enddo
	    enddo
	  enddo
	  zvel(ns)=zvel(ns)*betacub/dens(ns)
	enddo
	return
	end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccc Pressure ccccccccccccccccccccccccc
	subroutine pressure(betav, dens, temp, press)
! This routine has been modified to account for the fact that eta_ref=sqrt(k*T_ref/m) NOT  sqrt(2*k*T_ref/m)
	real phi(500,-10:10,-10:10,-10:10),dens(500),temp(500),press(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	do ns=1,nspace
	  press(ns)=dens(ns)*temp(ns)
	enddo
	return
	end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccc heat flux ccccccccccccccccccccccccc
	subroutine heatflux(betav,dens,xvel,qn)
	! this only calculates heat flux in the x direction (b/c 1D flow)
	real phi(500,-10:10,-10:10,-10:10),dens(500),xvel(500),qn(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	betasq=betav*betav
	betacub=betasq*betav
	do ns=1,nspace
	  qn(ns)=0.
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax
	      	Cx=(betav*float(i)-xvel(ns))
	        Cxsq=Cx**2
			Cysq=betasq*float(j*j)
			Czsq=betasq*float(k*k)
			Csq=Cxsq+Cysq+Czsq
			qn(ns)=qn(ns)+phi(ns,i,j,k)*Csq*Cx
	      enddo
	    enddo
	  enddo
	  qn(ns)=qn(ns)*betacub*dens(ns)/2.0
	enddo
	return
	end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccc shear stress ccccccccccccccccccccccccc
	subroutine tau11(betav,dens,xvel,press,tau)
	real phi(500,-10:10,-10:10,-10:10),dens(500),xvel(500)
	real press(500),tau(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	betasq=betav*betav
	betacub=betasq*betav
	do ns=1,nspace
	  tau(ns)=0.
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax
	        Cxsq=(betav*float(i)-xvel(ns))**2
			Cysq=betasq*float(j*j)
			Czsq=betasq*float(k*k)
			Csq=Cxsq+Cysq+Czsq
			tau(ns)=tau(ns)+phi(ns,i,j,k)*Cxsq
	      enddo
	    enddo
	  enddo
	  tau(ns)= -1.0*tau(ns)*betacub*dens(ns)
	  tau(ns)=tau(ns) + press(ns)
	enddo
	return
	end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccc krook collision operator ccccccccccccc
	subroutine Krookcoll_1(dens,xvel,yvel,zvel,temp,betav,deltat)
! Calculations assuming Maxwell molecules - scaled frequency=scaled density
	real phi(500,-10:10,-10:10,-10:10)
	real normfac
	real dens(500),xvel(500), yvel(500),zvel(500),temp(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	data pi/3.1415926536/
! Compute change in phi from a local Maxwellian at each point in space
	do ns=1,nspace
	  normfac=dens(ns)/sqrt((2.*pi*temp(ns))**3)
	  T_denom=0.5/temp(ns)
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax
		Cxsq=(betav*float(i)-xvel(ns))**2
		Cysq=(betav*float(j)-yvel(ns))**2		
		Czsq=(betav*float(k)-zvel(ns))**2
		Csq=Cxsq+Cysq+Czsq
		phieq=normfac*exp(-Csq*T_denom)
		! This line implements the actual krook operator
		delphi=deltat*dens(ns)*(phieq-phi(ns,i,j,k))
		phi(ns,i,j,k)=phi(ns,i,j,k)+delphi
	      enddo
	    enddo
	  enddo
	enddo
	return
	end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccc ESBGK collision operator ccccccccccccc
	subroutine ESBGKcoll_1(dens,xvel,yvel,zvel,temp,betav,deltat)
! Calculations assuming K_hat = 1
	real phi(500,-10:10,-10:10,-10:10)
	real eps_ij(3,3), M_ij(3,3), Ms_ij(3,3), CiCj(3,3)
	real normfac, Kf, lambda, Pr, trace_M, contract
	real dens(500),xvel(500), yvel(500),zvel(500),temp(500)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	data pi/3.1415926536/
	data kboltz/1.38064852e-23/ ! is this even necessary with normalization?

	betacub=betav**3
	Pr = 2.0/3.0
	lambda = -0.5
! Compute change in phi from a local Maxwellian at each point in space
	do ns=1,nspace
		Kf = 1.0 !temp(ns)*Pr
		normfac=dens(ns)/sqrt((2.*pi*temp(ns))**3)

		! zero out M_ij
		do l=1,3
			do m=1,3
				M_ij(l,m) = 0
			enddo
		enddo
	  
	  ! calculate M_ij and Ms_ij
	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax

	      	! there are 9 C_i*C_j terms
	      	CiCj(1,1) = (betav*float(i)-xvel(ns))**2
	      	CiCj(1,2) = (betav*float(i)-xvel(ns))
     &					*(betav*float(j)-yvel(ns))
	      	CiCj(1,3) = (betav*float(i)-xvel(ns))
     &     				*(betav*float(k)-zvel(ns))
	      	CiCj(2,1) = (betav*float(j)-yvel(ns))
     &     				*(betav*float(i)-xvel(ns))
	      	CiCj(2,2) = (betav*float(j)-yvel(ns))**2
	      	CiCj(2,3) = (betav*float(j)-yvel(ns))
     &     				*(betav*float(k)-zvel(ns))
	      	CiCj(3,1) = (betav*float(k)-zvel(ns))
     &     				*(betav*float(i)-xvel(ns))
	      	CiCj(3,2) = (betav*float(k)-zvel(ns))
     &     				*(betav*float(j)-yvel(ns))
	      	CiCj(3,3) = (betav*float(k)-zvel(ns))**2

	      	! these terms are summed up in M_ij
	      	do l=1,3
	      		do m=1,3
	      			M_ij(l,m) = M_ij(l,m) + phi(ns,i,j,k)*CiCj(l,m)
	      		enddo
	      	enddo

	      enddo
	    enddo
	  enddo
	  	! multiply the betav^3 factor that comes from the dV term
		do l=1,3
			do m=1,3
				M_ij(l,m) = M_ij(l,m)*betacub
			enddo
		enddo

		! Calculate Ms_ij
		trace_M = M_ij(1,1) + M_ij(2,2) + M_ij(3,3)
		do l=1,3
			do m=1,3
				Ms_ij(l,m) = M_ij(l,m)
				if (l .eq. m) then
					Ms_ij(l,m) = Ms_ij(l,m) - trace_M/3.0
				endif
			enddo
		enddo

		! calculate eps_ij
		do l=1,3
			do m=1,3
				eps_ij(l,m) = -lambda*Ms_ij(l,m)/(temp(ns)**2)
				if (l .eq. m) then
					eps_ij(l,m) = eps_ij(l,m) + 1.0/temp(ns)
				endif
			enddo
		enddo

	  do i=ivxmin,ivxmax
	    do j=ivymin,ivymax
	      do k=ivzmin,ivzmax

	      	! there are 9 C_i*C_j terms
	      	CiCj(1,1) = (betav*float(i)-xvel(ns))**2
	      	CiCj(1,2) = (betav*float(i)-xvel(ns))
     &					*(betav*float(j)-yvel(ns))
	      	CiCj(1,3) = (betav*float(i)-xvel(ns))
     &     				*(betav*float(k)-zvel(ns))
	      	CiCj(2,1) = (betav*float(j)-yvel(ns))
     &     				*(betav*float(i)-xvel(ns))
	      	CiCj(2,2) = (betav*float(j)-yvel(ns))**2
	      	CiCj(2,3) = (betav*float(j)-yvel(ns))
     &     				*(betav*float(k)-zvel(ns))
	      	CiCj(3,1) = (betav*float(k)-zvel(ns))
     &     				*(betav*float(i)-xvel(ns))
	      	CiCj(3,2) = (betav*float(k)-zvel(ns))
     &     				*(betav*float(j)-yvel(ns))
	      	CiCj(3,3) = (betav*float(k)-zvel(ns))**2

			! form the eps_ij*Ci*Cj contraction
			contract = 0.0
			do l=1,3
				do m=1,3
					contract = contract + eps_ij(l,m)*CiCj(l,m)
				enddo
			enddo

			! This line implements the actual ESBGK operator
	      	psi=normfac*exp(-0.5*contract)
			delphi=deltat*dens(ns)*Kf*(psi-phi(ns,i,j,k))
			phi(ns,i,j,k)=phi(ns,i,j,k)+delphi
		  enddo
		enddo
	  enddo


	enddo
	return
	end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccc convection step ccccccccccccccccccccccccc
	subroutine convect(alphax,betav,deltat)
! First order upwind scheme	

	real phi(500,-10:10,-10:10,-10:10), phif(-10:10,-10:10,-10:10)
	real phi1(500,-10:10,-10:10,-10:10)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	common /props/ phi
	common /coll/ phif
	cfl=betav*deltat/alphax
! First and last space points in domain are set by boundary conditions
	do ns=2,nspace-1
	  do j=ivymin,ivymax
	    do k=ivzmin,ivzmax
	      do i=ivxmin,-1
		! For i<0 upwind is on the right (ns+1)
		delphi=cfl*float(i)*(phi(ns+1,i,j,k)-phi(ns,i,j,k))
		phi1(ns,i,j,k)=phi(ns,i,j,k)-delphi
	      enddo
	      ! No convection for i=0
	      phi1(ns,0,j,k)=phi(ns,0,j,k)
	      do i=1,ivxmax
		! For i>0 upwind is on the left (ns-1)
		delphi=cfl*float(i)*(phi(ns,i,j,k)-phi(ns-1,i,j,k))
		phi1(ns,i,j,k)=phi(ns,i,j,k)-delphi
	      enddo
	    enddo
	  enddo
	enddo
! Left side boundary condition - for i>0 specular reflection of incoming (i<0) phi
! Normal convective update for i<0	
! We overwrite phi values at first space point for i>0 using -i values
	do j=ivymin,ivymax
	  do k=ivzmin,ivzmax
	    do i=ivxmin,0
		! For i<0 upwind is on the right (ns+1)
	      delphi=cfl*float(i)*(phi(2,i,j,k)-phi(1,i,j,k))
	      phi1(1,i,j,k)=phi(1,i,j,k)-delphi
	    enddo
	    do i=1,ivxmax
	      phi1(1,i,j,k)=phi(1,-i,j,k)
	    enddo
	  enddo
	enddo
! Right side boundary condition - freestream inflow for i<0
! Normal convective update for i>0	
	do j=ivymin,ivymax
	  do k=ivzmin,ivzmax
	    do i=ivxmin,0
	      phi1(nspace,i,j,k)=phif(i,j,k)
	    enddo
	    do i=1,ivxmax
		! For i>0 upwind is on the left (ns-1)
	      delphi=cfl*float(i)*(phi(nspace,i,j,k)-phi(nspace-1,i,j,k))
	      phi1(nspace,i,j,k)=phi(nspace,i,j,k)-delphi
	    enddo

	  enddo
	enddo
	
! Update the distribution function at all points in space after convection step
	do ns=1,nspace
	  do j=ivymin,ivymax
	    do k=ivzmin,ivzmax
	      do i=ivxmin,-1
		phi(ns,i,j,k)=phi1(ns,i,j,k)
	      enddo
	      ! No change computed for i=0, because there is no convection
	      do i=1,ivxmax
		phi(ns,i,j,k)=phi1(ns,i,j,k)		
	      enddo
	    enddo
	  enddo
	enddo
	
	return
	end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	

cccccccccccccccccccccccccccc Maxwellian distribution ccccccccccc 
c
c calculates a maxwellian distribution for initial state and returns it
c
	subroutine Maxwell_1(betav,ndf,uf,Tf,phif)
! This Maxwellian is defined using eta_ref=sqrt(k*T_ref/m) instead of sqrt(2*k*T_ref/m)
	common nspace,ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax
	data pi/3.1415926536/
	real norm, ndf, uf, Tf, phif(-10:10,-10:10,-10:10)
	T_denom=0.5/Tf
	norm=ndf/sqrt((2.*pi*Tf)**3)	
	do i=ivxmin,ivxmax
	  do j=ivymin,ivymax
	    do k=ivzmin,ivzmax
		Cxsq=(betav*float(i)-uf)**2
		Cysq=(betav*float(j))**2		
		Czsq=(betav*float(k))**2
		Csq=Cxsq+Cysq+Czsq
		phif(i,j,k)= norm*exp(-Csq*T_denom)
	    enddo
	  enddo
	enddo
	return
	end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
