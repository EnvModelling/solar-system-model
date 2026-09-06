	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>program to solve for the orbits of the planets
	!>Solar System Model (SSM): 
	!>Solves equations of motion in 3-D:
	!> <br><br>
	!> This program uses Newton's law of gravitation 
	!> which states that the force between two masses, 
	!> \f$m_i\f$ and \f$m_j\f$, separated by vector \f$\vec{r_{i,j}}\f$ is <br>
	!> \f{equation*}{ \vec{F}_{i,j}=-G\frac{m_i m_j}{r_{i,j}^3}\vec{r_{i,j}} \f}
	!><br><br>
	!> The coordinate system is cartesian based (i.e. x,y,z) 
	!> so Newton's Law of gravitation is written using this system. 
	!> The vector \f$\vec{r_{i,j}}\f$ can be written down in column form as: <br><br>
	!> \f{equation*}{ \vec{r_{i,j}}=\left(\begin{array}{c} x_j-x_i 
	!> \\y_j-y_i \\z_j-z_i \end{array} \right) \f}
	!> <br><br>
	!> Hence, the components of the force in the x,y and z direction are: <br>
	!>\f{equation*}{ \left(\begin{array}{c}
    !>  Fx_{i,j}\\
    !>  Fy_{i,j}\\
    !>  Fz_{i,j}
    !>  \end{array}\right)
    !> =-G\frac{m_i m_j}{\left({\left(x_j-x_i\right)^2 +
    !> \left(y_j-y_i\right)^2 +\left(z_j-z_i\right)^2}\right)^{3/2}}
    !>  \left(\begin{array}{c} x_j-x_i \\y_j-y_i \\z_j-z_i \end{array} \right) \f} <br>
    !> Newton's 2nd Law states for net force, \f$ \vec{F}\f$ , 
    !> applied to object mass, m causing the 
    !> object to accelerate with acceleration, a: <br>
    !> \f{equation*}{ \vec{F}=m\vec{a}\f} <br>
    !> Hence, to solve for the motion of each object (i.e. find the \f$x_i,y_i,z_i\f$ 
    !> positions 
    !> with time) we must solve three differential equations for each dimension: <br>
    !> \f{eqnarray*}{
    !> Fx_{i,j} &=& m_i\frac{d^2 x_i}{dt^2}\\
    !> Fy_{i,j} &=& m_i\frac{d^2 y_i}{dt^2}\\
    !> Fz_{i,j} &=& m_i\frac{d^2 z_i}{dt^2}
    !> \f} <br> <br>
    !> The code solves 1st order differential equations 
    !> (those having a first order derivative). Clearly the equations to solve are 
    !> 2nd order differential equations, but we turn them into 1st order differential 
    !> equations by specifying another 3 equations for the 3 components of velocity, 
    !> \f$vx_i,vy_i,vz_i\f$, which give the following 6 equations to solve for each 
    !> body in the solar system: <br><br>
    !> \f{eqnarray*}{
    !> Fx_{i,j} &=& m_i\frac{d vx_i}{dt}\\
    !> Fy_{i,j} &=& m_i\frac{d vy_i}{dt}\\
    !> Fz_{i,j} &=& m_i\frac{d vz_i}{dt}\\
    !> vx_i &=& \frac{d x_i}{dt}\\
    !> vy_i &=& \frac{d y_i}{dt}\\
    !> vz_i &=& \frac{d z_i}{dt}
    !> \f}
	!> compile using the Makefile (note requires netcdf) and then run using: <br>
	!> ./main.exe namelist.in
	!> <br><br>
	!> (namelist used for initialisation).
	!> <br><br>




	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>variables for the solar system model
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Module so variables can be easily passed to derivative function        !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	module consts_and_vars
		use, intrinsic :: iso_fortran_env
		use numerics_type
		implicit none

		integer(i4b), parameter :: n_bodies=10
		integer(i4b) :: neq=n_bodies*6,nn_interact=n_bodies
		real(wp), parameter :: c=2.99792458e8_wp ! speed of light

		real(wp), dimension(n_bodies) :: gm, &
		                 x,y,z, ux,uy,uz, meandist, lx,ly,lz,lt
		real(wp), dimension(n_bodies*6) :: yinit1
		integer(i4b), dimension(n_bodies) :: inds1, inds, &
		                                     interactions,interactions2

		logical, dimension(n_bodies) :: interact=.true.
		logical :: general_relativity=.true.
		real(wp), parameter :: G=6.6719e-11_wp, m=1.98855e30_wp
		real(wp) :: tt=0.e0_wp
	end module consts_and_vars
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main program of Solar System Model (SSM):
	program solar_system

	use netcdf
	use numerics_type
	use numerics, only : dvode
	use consts_and_vars, only : n_bodies, neq,gm,x,y,z,ux,uy,uz, yinit1, &
	                        meandist, inds1,inds,G,m,interactions,tt, &
	                        interact,nn_interact,interactions,interactions2, &
	                        general_relativity, lx,ly,lz,lt

	implicit none

	integer(i4b) :: i, j, allocatestatus
	real(wp), dimension(6*n_bodies) :: yinit, ydot, ysol, atol,rtol1
	real(wp), dimension(6*n_bodies) :: rpar
	integer(i4b), dimension(6*n_bodies) :: ipar

	! Defaults are in the same units as the namelist:
	! dt [days], tfinal [years]. They are converted to SI once after reading.
	real(wp) :: tout, dt=5e0_wp, rtol, tfinal, interval_io, tt_last
	integer(i4b) :: itol, itask, istate, iopt, ng, lrw, liw, mf, mflag
	real(wp), allocatable, dimension(:) :: rwork
	integer(i4b), allocatable, dimension(:) :: iwork

	! netcdf stuff
	character(len=100) :: outputfile01='output02.nc'
	integer(i4b), parameter :: nsav=10000
	real(wp), dimension(3,2,n_bodies) :: pos
	real(wp), dimension(3,2,n_bodies,nsav) :: pos_save
	real(wp), dimension(nsav) :: tt_save
	integer(i4b) :: ncid, varid, x_dimid, i_dimid, j_dimid, &
	                a_dimid, icur,isav, gm_varid, interact_varid
	integer(i4b), dimension(n_bodies) :: interact_int
	logical :: run_forward_in_time=.true.

	! namelist
	character(len=200) :: nmlfile=' '
	namelist /initial_state/ x,y,z,ux,uy,uz, gm
	namelist /run_vars/ outputfile01,run_forward_in_time, general_relativity, &
	                   dt,tfinal, interval_io, interact

	external solar01, jsolar01


	interactions=(/(i,i=1,n_bodies)/)

	! Initial data for the solar system (taken from JPL ephemeris)
	! the product of G and m for the bodies in the solar system
	!      gm=(/G*m/1.d9, &
	! Note: defaults are kept in the same km / km s-1 / km3 s-2 units as the
	! namelist, and are converted to SI once after the namelist is read.
	gm=(/1.327124400e11_wp, &
	  22032.09e0_wp, 324858.63e0_wp, 398600.440e0_wp, &
	  42828.3e0_wp, 126686511e0_wp, 37931207.8e0_wp, &
	  5793966e0_wp, 6835107e0_wp, 872.4e0_wp/)

	! The positions (x,y,z) and the velocities (vx,vy,va) of all the planets
	x=(/0.e0_wp, 1.563021412664830e+07_wp, -9.030189258080004e+07_wp, &
	 -1.018974476358996e+08_wp, -2.443763125844157e+08_wp, &
	 -2.35165468275322006e+08_wp, -1.011712827283427e+09_wp, &
	  2.934840841770302e+09_wp, 4.055112581124043e+09_wp, &
	  9.514009594170194e+08_wp/)

	y=(/0.e0_wp, 4.327888220902108e+07_wp, 5.802615456116644e+07_wp, &
	 1.065689158175689e+08_wp, 4.473211564076996e+07_wp, &
	 7.421837640432589e+08_wp, -1.077496255617324e+09_wp, &
	 6.048399137411513e+08_wp, -1.914578873112663e+09_wp, &
	 -4.776029500570151e+09_wp/)

	z=(/0.e0_wp, 2.102123103174893e+06_wp, 6.006513603716755e+06_wp, &
	 -3.381951053601424e+03_wp, 6.935657388967808e+06_wp, &
	 2.179850895804323e+06_wp, 5.901251900068215e+07_wp, &
	 -3.576451387567792e+07_wp, -5.400973716179796e+07_wp, &
	 2.358627841705075e+08_wp/)

	ux=(/0.e0_wp, -5.557001175482630e+01_wp, -1.907374632532257e+01_wp, &
	 -2.201749257051057e+01_wp, -3.456935754608896e+00_wp, &
	 -1.262559929908801e+01_wp, 6.507898648442419e+00_wp, &
	 -1.433852081777671e+00_wp, 2.275119229131818e+00_wp, &
	 5.431808363374300e+00_wp/)

	uy=(/0.e0_wp, 1.840863017229157e+01_wp, -2.963461693326599e+01_wp, &
	 -2.071074857788741e+01_wp, -2.176307370133160e+01_wp, &
	 -3.332552395475581e+00_wp, -6.640809674126991e+00_wp, &
	 6.347897341634990e+00_wp, 4.942356914027413e+00_wp, &
	 -2.387056445508962e-02_wp/)

	uz=(/0.e0_wp, 6.602621285552567e+00_wp, 6.946391255404438e-01_wp, &
	 1.575245213712245e-03_wp, -3.711433859326417e-01_wp, &
	 2.962741332356101e-01_wp, -1.434198106014633e-01_wp, &
	 4.228261484335974e-02_wp, -1.548950389954096e-01_wp, &
	 -1.551877289694926e+00_wp/)

	tfinal=5.e2_wp


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! namelist                                                               !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call getarg(1,nmlfile)
	open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
	read(8,nml=initial_state)

	! Convert initial-state namelist values from km / km s-1 / km3 s-2 to SI.
	x=x*1e3_wp
	y=y*1e3_wp
	z=z*1e3_wp
	ux=ux*1e3_wp
	uy=uy*1e3_wp
	uz=uz*1e3_wp
	gm=gm*1e9_wp

	read(8,nml=run_vars)
	dt=dt*86400.e0_wp
	tfinal=tfinal*365.25e0_wp*86400.e0_wp
	close(8)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	nn_interact=count(interact) ! number of true elements
	j=1
	do i=1,n_bodies
	  if(interact(i)) then
	     interactions2(j)=interactions(i)
	     j=j+1
	  end if
	end do

	if(.not.run_forward_in_time) then
	  ux=-ux
	  uy=-uy
	  uz=-uz
	endif

	! mean distance from the sun
	meandist=(/7e8_wp, 5.79e10_wp, 1.082e11_wp, &
	      1.496e11_wp, 2.279e11_wp, 7.783e11_wp, &
	      1.426e12_wp, 2.871e12_wp, 4.497e12_wp, 5.914e12_wp/)

	inds1=1

	! angular momenta per unit mass
	lx=uy*z-uz*y
	ly=uz*x-ux*z
	lz=ux*y-uy*x
	lt=sqrt(lx*lx+ly*ly+lz*lz)


	! open netcdf file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check( nf90_create(outputfile01, nf90_clobber, ncid) )
	call check( nf90_def_dim(ncid, "times", nf90_unlimited, x_dimid) )
	call check( nf90_def_dim(ncid, "n_bodies", n_bodies, i_dimid) )
	call check( nf90_def_dim(ncid, "n_dims", 3 , j_dimid) )
	! close and free up any buffers
	! The file is now deliberately kept open, so the old close/open pair is
	! no longer required here.

	! open file for writing:
	! nf90_create has already opened the file for writing.
	! define mode
	! nf90_create also leaves the file in define mode.



	! define a variable TIME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check( nf90_def_var(ncid, "time", nf90_float, & ! was nf90_double
	                  (/x_dimid/), varid) )
	! get id to a_dimid
	call check( nf90_inq_varid(ncid, "time", a_dimid) )
	! units
	call check( nf90_put_att(ncid, a_dimid, "units", "s") )


	! define a variable POS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check( nf90_def_var(ncid, "pos", nf90_float, &
	                  (/j_dimid, i_dimid,x_dimid/), varid) )
	! get id to a_dimid
	call check( nf90_inq_varid(ncid, "pos", a_dimid) )
	! units
	call check( nf90_put_att(ncid, a_dimid, "units", "m") )


	! define a variable VEL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check( nf90_def_var(ncid, "vel", nf90_float, &
	                  (/j_dimid, i_dimid,x_dimid/), varid) )
	! get id to a_dimid
	call check( nf90_inq_varid(ncid, "vel", a_dimid) )
	! units
	call check( nf90_put_att(ncid, a_dimid, "units", "m s-1") )


	! define static gravitational parameter GM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check(nf90_def_var(ncid, "gm", nf90_float, (/i_dimid/), gm_varid))
	call check(nf90_put_att(ncid, gm_varid, "units", "m3 s-2"))
	call check(nf90_put_att(ncid, gm_varid, "long_name", &
	                        "gravitational parameter"))

	! define the gravity-source mask!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check(nf90_def_var(ncid, "gravity_source", nf90_int, &
	                        (/i_dimid/), interact_varid))
	call check(nf90_put_att(ncid, interact_varid, "long_name", &
	                        "1 if body is enabled in interact"))

	! define run metadata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check(nf90_put_att(ncid, nf90_global, "title", &
	                        "Solar System Model output"))
	call check(nf90_put_att(ncid, nf90_global, "initial_epoch", &
	                        "2014-02-02 20:18:00"))
	call check(nf90_put_att(ncid, nf90_global, "body_order", &
	                        "Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,"// &
	                        "Uranus,Neptune,Pluto"))
	call check(nf90_put_att(ncid, nf90_global, "namelist", trim(nmlfile)))
	call check(nf90_put_att(ncid, nf90_global, "run_forward_in_time", &
	                        merge(1_i4b,0_i4b,run_forward_in_time)))
	call check(nf90_put_att(ncid, nf90_global, "general_relativity", &
	                        merge(1_i4b,0_i4b,general_relativity)))
	call check(nf90_put_att(ncid, nf90_global, "dt_days", dt/86400._wp))
	call check(nf90_put_att(ncid, nf90_global, "tfinal_years", &
	                        tfinal/(365.25_wp*86400._wp)))
	call check(nf90_put_att(ncid, nf90_global, "interval_io_years", interval_io))

#ifdef heliocentric
	call check(nf90_put_att(ncid, nf90_global, "heliocentric", 1_i4b))
#else
	call check(nf90_put_att(ncid, nf90_global, "heliocentric", 0_i4b))
#endif

	! END DEFINE MODE. 
	! TELL netCDF WE ARE DONE DEFINING METADATA.
	call check( nf90_enddef(ncid) )

	! Write variables that do not change with time.
	interact_int=merge(1_i4b,0_i4b,interact)
	call check( nf90_put_var(ncid, gm_varid, gm) )
	call check( nf90_put_var(ncid, interact_varid, interact_int) )
	call check( nf90_sync(ncid) )

	! The original version closed the file here.  It is now kept open so that
	! buffered trajectory writes do not repeatedly open and close the file.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set-up variables for the ode solver                                    !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,n_bodies
	  yinit(1+(i-1)*6:6*i)=(/x(i), y(i), z(i), ux(i), uy(i), uz(i)/)
	end do
	yinit1=yinit


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! set up variables to pass to dvode                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	tt=0.e0_wp
	ysol=yinit
	tout=min(tt+dt,tfinal)
	!      tfinal=365.25d0*86400.d0*1.d6
	!itol=2 ! both relative and absolute error convergence
	itol=4 ! both relative and absolute error convergence
	!      rtol=1.d-8
	!rtol1=1.e-6_wp
	!rtol1(7:12)=1e-10_wp
	rtol1=1e-12_wp
	atol(:)=1.e-6_wp ! nearest micron and micron / second (i.e. negligible)
	itask=1
	istate=1
	iopt=1 ! optional input
	mf=10 !22 !
	lrw=22+9*neq+2*neq**2
	liw=30*neq

	allocate(rwork(lrw), stat=allocatestatus)
	if (allocatestatus /= 0) stop "*** not enough memory ***"

	allocate(iwork(liw), stat=allocatestatus)
	if (allocatestatus /= 0) stop "*** not enough memory ***"

	iwork(6)=10000 ! max steps
	iwork(7)=10  ! max messages printed per problem
	iwork(5)=12 ! 5   ! order
	rwork(5)=0.e0_wp ! initial time-step
	rwork(6)=dt   ! max time-step
	rwork(7)=0.e0_wp ! min time-step allowed
	rwork(14)=2.e0_wp ! tolerance scale factor
	mflag=0
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!CALL XSETF(MFLAG)


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! solve system of equations calling dvode                                !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	icur=1
	isav=1
	istate=1
	interval_io=interval_io*365.25_wp*86400._wp ! put in seconds
	tt_last=-interval_io

	do while (tt.lt.tfinal)
	  do while(tt.lt.tout)
	    call dvode(solar01,neq,ysol,tt,tout,itol,rtol1,atol,itask,istate, &
	               iopt,rwork,lrw,iwork,liw,jsolar01,mf,rpar,ipar)
	    istate=2
	  enddo

	  ! Do not request an output time beyond the requested end time.
	  tout=min(tt+dt,tfinal)

	  tt_save(isav)=tt
	  pos_save(1:3,1:2,1:n_bodies,isav:isav)= &
	      reshape(ysol(1:neq),(/3,2,n_bodies,1/))

	  if(tt-tt_last >= interval_io) then
	     ! position in run:
	     print *,'Model is up to here: ',tt/(365.25*86400), &
	             ' done of ',tfinal/(365.25*86400)
	     tt_last=tt
	  endif

	  if (isav.eq.nsav.or.tt.ge.tfinal) then

	    ! open the netcdf file (do not overwrite if it exists)
	    ! The file is already open; keeping it open avoids repeated open/close calls.

	    ! write netcdf variables
	    call check(nf90_inq_varid(ncid, "time", varid))
	    call check(nf90_put_var(ncid, varid, tt_save(1:isav), &
	                            start=(/ICUR/)))

	    call check(nf90_inq_varid(ncid, "pos", varid))
	    call check(nf90_put_var(ncid, varid, &
	                            pos_save(1:3,1,1:n_bodies,1:isav), &
	                            start=(/1,1,ICUR/)))

	    call check(nf90_inq_varid(ncid, "vel", varid))
	    call check(nf90_put_var(ncid, varid, &
	                            pos_save(1:3,2,1:n_bodies,1:isav), &
	                            start=(/1,1,ICUR/)))

	    ! close the netCDF file:
	    ! Keep the file open between buffered writes; just flush it to disk here.
	    call check( nf90_sync(ncid) )

	    icur=icur+isav
	    isav=0
	  end if

	  isav=isav+1
	end do

	! close the netCDF file:
	call check( nf90_close(ncid) )
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end program solar_system
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>provides derivatives for the RHS of the equations
	!>@param[inout] neq: number of odes to solve
	!>@param[inout] tt: time variable
	!>@param[inout] y, ydot: solution and derivative
	!>@param[inout] rpar, ipar: real and integer work arrays
	subroutine solar01(neq, tt, y, ydot, rpar, ipar)

	use consts_and_vars, only : gm, n_bodies,inds,interactions, yinit1, &
	                           interact,interactions2,nn_interact,c, &
	                           general_relativity,lx,ly,lz,lt
	use numerics_type
	implicit none

	real(wp), intent(inout) :: tt
	real(wp), intent(inout), dimension(neq) :: y, ydot
	integer(i4b), intent(inout) :: neq
	real(wp), intent(inout),dimension(6*n_bodies) :: rpar
	integer(i4b), intent(inout),dimension(6*n_bodies) :: ipar

	! locals
	integer(i4b) :: i,n,k,j
	real(wp), dimension(n_bodies-1) :: r2
	logical, dimension(n_bodies) :: interact2
	real(wp) :: lx1,ly1,lz1, lt1, r2_sun

	ydot=0._wp


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Loop over all bodies                                                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,n_bodies
	  interact2=interact

	  ! set the rate of change of position to the y-values in the 2nd order diff eq
	  ! This is done before checking the interaction list so source-free bodies
	  ! continue to move with their existing velocity.
	  ydot((i-1)*6+1)=y((i-1)*6+4)
	  ydot((i-1)*6+2)=y((i-1)*6+5)
	  ydot((i-1)*6+3)=y((i-1)*6+6)

	  n=count(interactions2(1:nn_interact).ne.i)

	  if (n.gt.0) then
	    ! Do not interact with self:
	    if (i.eq.interactions(i)) interact2(i)=.false.

	    ! inds is the body that i interacts with:
	    inds(1:n)=pack(interactions,interact2)

	    ! inds indexes the body that body i interacts with.
	    ! however, the x,y,z's are separated by 6


	    ! square of distance between this planet and the other objects:
	    r2(1:n)=(y((i-1)*6+1)-y((inds(1:n)-1)*6+1))**2 + &
	            (y((i-1)*6+2)-y((inds(1:n)-1)*6+2))**2 + &
	            (y((i-1)*6+3)-y((inds(1:n)-1)*6+3))**2

	    ! inverse square law between body i and the rest of them
	    ! dux/dt:
	    ydot((i-1)*6+4)=ydot((i-1)*6+4)- &
	      sum(gm(inds(1:n))/r2(1:n) * &
	      (y((i-1)*6+1)-y((inds(1:n)-1)*6+1))/sqrt(r2(1:n)))

	    ! duy/dt:
	    ydot((i-1)*6+5)=ydot((i-1)*6+5)- &
	      sum(gm(inds(1:n))/r2(1:n) * &
	      (y((i-1)*6+2)-y((inds(1:n)-1)*6+2))/sqrt(r2(1:n)))

	    ! duz/dt:
	    ydot((i-1)*6+6)=ydot((i-1)*6+6)- &
	      sum(gm(inds(1:n))/r2(1:n) * &
	      (y((i-1)*6+3)-y((inds(1:n)-1)*6+3))/sqrt(r2(1:n)))
	  end if


	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! extra term for general relativity - do not apply to sun's motion                  !
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  if (general_relativity .and. (i.gt.1)) then
	    ! terms for General Relativity - 
	    ! http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node116.html
	    n=1
	    ! note that l^2=G*M0*r
	    ! The correction is now explicitly calculated relative to body 1 (the Sun),
	    ! rather than whichever Newtonian interaction happened to be first in inds.
	    r2_sun=(y((i-1)*6+1)-y(1))**2 + &
	           (y((i-1)*6+2)-y(2))**2 + &
	           (y((i-1)*6+3)-y(3))**2

	    ! dux/dt:
	    ydot((i-1)*6+4)=ydot((i-1)*6+4)- &
	      3.e0_wp*gm(1)**2/c**2/r2_sun**2 * &
	      (y((i-1)*6+1)-y(1))

	    ! duy/dt:
	    ydot((i-1)*6+5)=ydot((i-1)*6+5)- &
	      3.e0_wp*gm(1)**2/c**2/r2_sun**2 * &
	      (y((i-1)*6+2)-y(2))

	    ! duz/dt:
	    ydot((i-1)*6+6)=ydot((i-1)*6+6)- &
	      3.e0_wp*gm(1)**2/c**2/r2_sun**2 * &
	      (y((i-1)*6+3)-y(3))
	  end if
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





	  ! set the rate of change of position to the y-values in the 2nd order diff eq
	  ! These three assignments are performed near the top of the loop now so that
	  ! a body with no active gravity sources still moves ballistically.
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef heliocentric
	! adjust positions and accelerations, so the sun remains at (0,0,0)
	ydot(7:n_bodies*6)=ydot(7:n_bodies*6)- &
	      reshape(spread(ydot(1:6),2,n_bodies-1),(/6*(n_bodies-1)/))
	ydot(1:6)=0.d0
#endif

	end subroutine solar01
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>jacobian needs to be defined - does not need to do anything
	!>@param[inout] neq: number of odes to solve
	!>@param[inout] tt: time variable
	!>@param[inout] y, ydot: solution and derivative
	!>@param[inout] rpar, ipar: real and integer work arrays
	subroutine jsolar01(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

		use consts_and_vars, only : n_bodies
		use numerics_type

		implicit none
		real(wp), intent(in) :: t
		real(wp), dimension(neq), intent(inout) :: y
		real(wp), dimension(nrpd, neq), intent(inout) :: pd
		integer(i4b), intent(inout) :: neq, ml, mu, nrpd
		real(wp), intent(inout),dimension(6*n_bodies) :: rpar
		integer(i4b), intent(inout),dimension(6*n_bodies) :: ipar

	end subroutine jsolar01
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>helper routine for netcdf   
	!>@param[in] status of netcdf file pointer
	subroutine check(status)
		use netcdf
		use numerics_type
		integer(i4b), intent(in) :: status

		if(status /= nf90_noerr) then
			print *, trim(nf90_strerror(status))
			stop "stopped"
		end if
	end subroutine check
