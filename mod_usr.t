!===============================================================================
! 1D CAK wind module for launching a spherically symmetric line-driven wind
! from the stellar surface
!
! Coded up by Flo for his KU Leuven PhD thesis 2018/2022 using a version of
! Nico's routine of his KU Leuven master thesis 2017/2018
!-------------------------------------------------------------------------------
! Options for wind to be specified in usr.par file:
!   ifrc = 0  : radially streaming CAK wind
!   ifrc = 1  : finite disk corrected CAK wind
!   ifrc = 2  : finite disk + opacity cut-off
!
! (April 2019) -- Flo
!   > major cleaning and deleting unnecessary defined variables
!   > included opacity cut-off option (for LDI sims)
!
! (August 2019) -- Flo
!   > fixed problems with unit conversion in original code of Nico
!   > extensive testing + further cleaning + commenting
!
! (September 2019) -- Flo
!   > last fixes in line-force computation
!   > removed effective gravity from output (not interesting, only for testing)
!
! (July 2020) -- Flo
!   > renaming of my_unit_mass/lum/ggrav in order to compile without conflict as
!     AMRVAC v.3 has a new particle module (with protected var my_unit_mass)
!
! (August 2020) -- Flo
!   > implemented exact dv/dr algorithm for non-uniform grids + some minor
!     performance improvements
!   > changed CAK source computation from CGS to unitless + saved unitless to
!     be in line with other unitless output
!
! (February 2021) -- Flo
!   > implemented special outer boundary conditions as well, some kinks occured
!     in gcak + fdfac due to abrupt change in dv/dr in outer boundary
!     essentially the AMRVAC 'cont' takes zero gradient but constant slope
!     extrapolation alleviates the problem
!
!===============================================================================

module mod_usr

  ! Include a physics module
  use mod_HD

  implicit none

  ! The usual suspects
  real(8), parameter :: msun=1.989d33, lsun=3.827d33, rsun=6.96d10, &
                        Ggrav=6.67d-8, kappae=0.34d0

  ! Stellar parameters: luminosity, mass, radius, surface density, eff. temp.
  real(8) :: lstar, mstar, rstar, rhobound, twind

  ! Unit quantities that are handy: gravitational constant, luminosity, mass
  real(8) :: my_unit_ggrav, my_unit_lum, my_unit_mass

  !
  ! Wind parameters: CAK alpha, Gayley Qbar + Qmax, beta power velocity law,
  !                  Eddington gamma, escape speed, CAK + fd mass-loss rate,
  !                  terminal wind speed, sound speed
  real(8) :: alpha, Qbar, Qmax, beta, gammae, vesc, mdot, mdotfd, vinf, asound

  ! Dimensionless variables of relevant variables
  real(8) :: dlstar, dmstar, drstar, drhobound, dtwind, dkappae, dvesc, &
             dvinf, dmdot, dasound, dclight, dGgrav, dgammae

  ! log(g), eff. log(g), scale height, year (s), mean molecular weight
  real(8) :: logg, logge, heff, year=3.15567d7, mumol

  ! Time-step after first iteration taking into account CAK radiation force
  real(8) :: new_timestep

  ! New variables to store in conservative variables array 'w'
  integer :: my_gcak, my_fdf

  ! Wind options and avoid magic numbers for wind option
  integer :: ifrc
  integer, parameter :: radstream = 0, finite_disk = 1, finite_disk_cutoff = 2

  character(len=8) :: todayis

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system: 1D spherical
    call set_coordinate_system("spherical")

    ! Read in the required stellar + wind params from the usr.par file
    call usr_params_read(par_files)

    !
    ! Choose independent normalization units, only 3 have to be specified:
    !     (length,temp,ndens) or (length,vel,ndens)
    ! Numberdensity chosen such that unit density becomes boundary density
    ! IMPORTANT: simulations cannot be run in CGS units <=> unit quantities = 1
    !
    unit_length        = rstar                                        ! cm
    unit_temperature   = twind                                        ! K
    unit_numberdensity = rhobound/((1.0d0+4.0d0*He_abundance)*mp_cgs) ! g cm^-3

    ! Activate the physics module
    call HD_activate()

    ! Initialize user (global) quantities after initializing AMRVAC
    usr_set_parameters => initglobaldata_usr

    ! Initial conditions
    usr_init_one_grid => initial_conditions

    ! Special Boundary conditions
    usr_special_bc => special_bound

    ! CAK force computation
    usr_source => CAK_source

    ! Adjusted timestep to account for CAK force
    usr_get_dt => special_dt

    ! Define user custom variables to store extra computations in output
    my_gcak = var_set_extravar("gcak", "gcak")
    my_fdf  = var_set_extravar("fdfac", "fdfac")

  end subroutine usr_init

!===============================================================================

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ mstar, lstar, rstar, twind, rhobound, alpha, &
                          Qbar, Qmax, beta, ifrc

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    end do

    ! Scale to cgs units
    lstar = lstar * lsun
    mstar = mstar * msun
    rstar = rstar * rsun

  end subroutine usr_params_read

!===============================================================================

  subroutine initglobaldata_usr
    !
    ! Compute some quantities of interest (in CGS) before making unitless
    !

    ! Stellar structure
    gammae = kappae * lstar/(4.0d0*dpi * Ggrav * mstar * const_c)
    logg   = log10(Ggrav * mstar/rstar**2.0d0)
    logge  = logg * (1.0d0 - gammae)**0.5d0
    mumol  = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    asound = (twind * kB_cgs/(mumol * mp_cgs))**0.5d0
    heff   = asound**2.0d0 / 10.d0**logge

    ! Wind quantities in CAK theory
    Qmax   = Qmax * Qbar
    vesc   = (2.0d0 * Ggrav * mstar * (1.0d0 - gammae)/rstar)**0.5d0
    vinf   = vesc * (alpha/(1.0d0 - alpha))**0.5d0
    mdot   = lstar/const_c**2.0d0 * alpha/(1.0d0 - alpha) &
              * (Qbar * gammae/(1.0d0 - gammae))**((1.0d0 - alpha)/alpha)
    mdotfd = mdot/(1.0d0 + alpha)**(1.0d0/alpha)

    ! Make all relevant variables dimensionless
    call make_dimless_vars()

    if (mype == 0) then
      call date_and_time(todayis)
      print*, 'MPI-AMRVAC simulation ran on ', &
                    todayis(7:8), '/', todayis(5:6), '/', todayis(1:4)
      print*
      print*, '======================'
      print*, '   Unity quantities   '
      print*, '======================'
      print*, 'unit length        = ', unit_length
      print*, 'unit density       = ', unit_density
      print*, 'unit velocity      = ', unit_velocity
      print*, 'unit numberdensity = ', unit_numberdensity
      print*, 'unit pressure      = ', unit_pressure
      print*, 'unit temperature   = ', unit_temperature
      print*, 'unit time          = ', unit_time
      print*
      print*, '==============================================='
      print*, '   Stellar and wind parameters in CGS units    '
      print*, '==============================================='
      print*, 'L/Lsun                 = ', lstar/lsun
      print*, 'M/Msun                 = ', mstar/msun
      print*, 'R/Rsun                 = ', rstar/rsun
      print*, 'Twind                  = ', twind
      print*, 'Mean molecular weight  = ', mumol
      print*, 'log(g)                 = ', logg
      print*, 'eff. log(g)            = ', logge
      print*, 'eff. scale height heff = ', heff
      print*, 'heff/Rstar             = ', heff/rstar
      print*, 'Eddington gamma        = ', gammae
      print*
      print*, 'adiabatic gamma = ', hd_gamma
      print*, 'alpha           = ', alpha
      print*, 'Qbar            = ', Qbar
      print*, 'Qmax/Qbar       = ', Qmax/Qbar
      print*, 'beta            = ', beta
      print*, 'asound          = ', asound
      print*, 'eff. vesc       = ', vesc
      print*, 'vinf            = ', vinf
      print*
      print*, 'wind option            = ', ifrc
      print*, '   0 : radial stream CAK '
      print*, '   1 : CAK + fd          '
      print*, '   2 : CAK + cut-off     '
      print*, 'surface density        = ', rhobound
      print*, 'analytic Mdot CAK      = ', mdot * (year/msun)
      print*, '... with FD correction = ', mdotfd * (year/msun)
      print*
      print*, '========================================'
      print*, '    Dimensionless AMRVAC quantities     '
      print*, '========================================'
      print*, 'Extra computed unit quantities:'
      print*, '   unit Lum  = ', my_unit_lum
      print*, '   unit Mass = ', my_unit_mass
      print*, '   unit Grav = ', my_unit_ggrav
      print*, 'Lstar        = ', dlstar
      print*, 'Mstar        = ', dmstar
      print*, 'Rstar        = ', drstar
      print*, 'Twind        = ', dtwind
      print*, 'Edd. gamma   = ', dgammae
      print*, 'rhobound     = ', drhobound
      print*, 'Mdot         = ', dmdot
      print*, 'alpha        = ', alpha
      print*, 'Qbar         = ', Qbar
      print*, 'Qmax         = ', Qmax/Qbar
      print*, 'kappae       = ', dkappae
      print*, 'beta         = ', beta
      print*, 'asound       = ', dasound
      print*, 'eff. vesc    = ', dvesc
      print*, 'vinf         = ', dvinf
      print*, 'clight       = ', dclight
      print*, 'Ggrav        = ', dGgrav
      print*
    endif

  end subroutine initglobaldata_usr

!===============================================================================

  subroutine initial_conditions(ixI^L,ixO^L,w,x)
    !
    ! Initial conditions start from beta velocity law
    !

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: sfac

    ! Small offset (asound/vinf) to avoid starting at terminal wind speed
    sfac = 1.0d0 - 1.0d-3**(1.0d0/beta)

    where (x(ixI^S,1) >= drstar)
       ! Set initial velocity field to beta law
       w(ixI^S,mom(1)) = dvinf * ( 1.0d0 - sfac * drstar / x(ixI^S,1) )**beta

       ! Set initial density
       w(ixI^S,rho_) = dmdot / (4.0d0*dpi * x(ixI^S,1)**2.0d0 * w(ixI^S,mom(1)))
    endwhere

    ! Convert hydro vars to conserved to let AMRVAC do computations
    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

!===============================================================================

  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)
    !
    ! Modified boundary values at inner (star) and outer radial boundary
    !

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: i

    select case (iB)
    case(1)

      ! Convert hydro vars to primitive
      call hd_to_primitive(ixI^L,ixI^L,w,x)

      w(ixB^S,rho_) = drhobound

      ! Radial velocity field (constant slope extrapolation)
      do i = ixBmax1,ixBmin1,-1
        if (i == ixBmax1) then
          w(i,mom(1)) = w(i+1,mom(1)) - (w(i+2,mom(1)) - w(i+1,mom(1))) &
                            * (x(i+1,1) - x(i,1))/(x(i+2,1) - x(i+1,1))
        else
          w(i,mom(1)) = w(i+1,mom(1))
        endif

      enddo

      ! Prohibit ghosts to be supersonic, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), dasound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -dasound)

      ! Convert hydro vars back to conserved to let AMRVAC do computations
      call hd_to_conserved(ixI^L,ixI^L,w,x)

    case(2)
      ! Constant slope extrapolation of all hydro vars

      ! Convert hydro vars to primitive
      call hd_to_primitive(ixI^L,ixI^L,w,x)

      w(ixB^S,rho_) = w(ixBmin1-1,rho_) * (x(ixBmin1-1,1) / x(ixB^S,1))**2.0d0

      do i = ixBmin1,ixBmax1
        w(i,mom(1)) = w(i-1,mom(1)) + (w(ixBmin1-1,mom(1)) - w(ixBmin1-2,mom(1)))
      enddo

      ! Convert hydro vars back to conserved to let AMRVAC do computations
      call hd_to_conserved(ixI^L,ixI^L,w,x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!===============================================================================

  subroutine CAK_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    !
    ! Compute the analytical CAK line-force
    !

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L, iw^LIM
    real(8), intent(in)    :: qdt, qtC, qt
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: dvdr_up(ixO^S), dvdr_down(ixO^S), dvdr_cent(ixO^S), dvdr(ixO^S)
    real(8) :: gcak(ixO^S), geff(ixO^S)
    real(8) :: vr(ixI^S), rho(ixI^S)
    real(8) :: beta_fd(ixO^S), fdfac(ixO^S), timedum(ixO^S)
    real(8) :: taum(ixO^S), taumfac(ixO^S)
    real(8) :: fac, fac1, fac2
    integer :: jx^L, hx^L

    !========================================================================
    ! Convert to primitives

    ! Define the time-centred, local velocity from the local momentum
    vr(ixI^S) = wCT(ixI^S,mom(1)) / wCT(ixI^S,rho_)

    ! Time-centred density
    rho(ixI^S) = wCT(ixI^S,rho_)

    !========================================================================

    !
    ! Make new indices covering whole grid by increasing +1 (j) and decreasing
    ! by -1 (h). Special, fancy syntax that AMRVAC understands
    !
    jx^L=ixO^L+kr(1,^D);
    hx^L=ixO^L-kr(1,^D);

    ! Get dv/dr on non-uniform grid according to Sundqvist & Veronis (1970)
    ! Forward difference
    dvdr_up(ixO^S)   = (x(ixO^S,1) - x(hx^S,1)) * vr(jx^S) &
                        / ((x(jx^S,1) - x(ixO^S,1)) * (x(jx^S,1) - x(hx^S,1)))

    ! Backward difference
    dvdr_down(ixO^S) = -(x(jx^S,1) - x(ixO^S,1)) * vr(hx^S) &
                        / ((x(ixO^S,1) - x(hx^S,1)) * (x(jx^S,1) - x(hx^S,1)))

    ! Central difference
    dvdr_cent(ixO^S)  = (x(jx^S,1) + x(hx^S,1) - 2.0d0*x(ixO^S,1)) * vr(ixO^S) &
                        / ((x(ixO^S,1) - x(hx^S,1)) * (x(jx^S,1) - x(ixO^S,1)))

    ! Total gradient (in CAK this has to be >0, otherwise stagnant flow)
    dvdr(ixO^S) = abs(dvdr_down(ixO^S) + dvdr_cent(ixO^S) + dvdr_up(ixO^S))

    ! Finite disk factor parameterisation (Owocki & Puls 1996)
    beta_fd(ixO^S) = (1.0d0 - vr(ixO^S) / (x(ixO^S,1) * dvdr(ixO^S))) &
                      * (drstar/x(ixO^S,1))**2.0d0

    ! Check the finite disk array and determine finite disk factor
    if (ifrc == 0) then
      fdfac(ixO^S) = 1.0d0
    else
      where (beta_fd >= 1.0d0)
        fdfac = 1.0d0/(1.0d0 + alpha)
      elsewhere (beta_fd < -1.0d10)
        fdfac = abs(beta_fd)**alpha / (1.0d0 + alpha)
      elsewhere (abs(beta_fd) > 1.0d-3)
        fdfac = (1.0d0 - (1.0d0 - beta_fd)**(1.0d0 + alpha)) &
                    / (beta_fd*(1.0d0 + alpha))
      elsewhere
        fdfac = 1.0d0 - 0.5d0*alpha*beta_fd &
                        * (1.0d0 + 1.0d0/3.0d0 * (1.0d0 - alpha)*beta_fd)
      endwhere
    endif

    ! Fill up the empty finite disk correction array variable
    w(ixO^S,my_fdf) = fdfac(ixO^S)

    ! Calculate CAK line-force
    fac1 = 1.0d0/(1.0d0 - alpha) * dkappae * dlstar*Qbar/(4.0d0*dpi * dclight)
    fac2 = 1.0d0/(dclight * Qbar * dkappae)**alpha
    fac  = fac1 * fac2

    gcak(ixO^S) = fac/x(ixO^S,1)**2.0d0 * (dvdr(ixO^S)/rho(ixO^S))**alpha

    ! Based on wind option proceed now
    if (ifrc == radstream) then
      gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S)
    else if (ifrc == finite_disk) then
      gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S)
    else if (ifrc == finite_disk_cutoff) then
      taum(ixO^S) = dkappae * dclight * Qmax * rho(ixO^S)/dvdr(ixO^S)

      taumfac(ixO^S) = ((1.0d0 + taum(ixO^S))**(1.0d0 - alpha) - 1.0d0) &
                          / taum(ixO^S)**(1.0d0 - alpha)

      gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S) * taumfac(ixO^S)
    else
      stop 'Error in wind option, take a valid ifrc {0,1,2}'
    endif

    ! Fill the empty CAK array variable
    w(ixO^S,my_gcak) = gcak(ixO^S)

    ! Effective gravity computation
    geff(ixO^S) = - dGgrav * dmstar * (1.0d0 - dgammae)/x(ixO^S,1)**2.0d0

    ! Update conservative vars: w = w + qdt*gsource
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) &
                      + qdt * (gcak(ixO^S) + geff(ixO^S))*wCT(ixO^S,rho_)

    ! Define a new time-step corrected for continuum and line-acceleration
    timedum(ixO^S) = (x(jx^S,1) - x(ixO^S,1)) / (gcak(ixO^S) + geff(ixO^S))
    new_timestep   = 0.3d0 * minval( abs(timedum(ixO^S)) )**0.5d0

  end subroutine CAK_source

!===============================================================================

  subroutine special_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    !
    ! After first iteration assign the new time-step of computation CAK force
    !

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: dx^D, x(ixI^S,1:ndim)
    real(8), intent(in)    :: w(ixI^S,1:nw)
    real(8), intent(inout) :: dtnew

    if (it >= 1) then
      dtnew = new_timestep
    endif

  end subroutine special_dt

!===============================================================================

  subroutine make_dimless_vars()
    !
    ! Normalize quantities in use to unit quantities defined and computed
    ! These quantities are actually used by AMRVAC in its computations!
    !

    ! From the AMRVAC unit variables compute some extra relevant for us
    my_unit_ggrav = unit_density * unit_time**2.0d0
    my_unit_lum   = unit_density * unit_length**5.0d0 / unit_time**3.0d0
    my_unit_mass  = unit_density * unit_length**3.0d0

    drhobound = rhobound/unit_density
    dlstar    = lstar/my_unit_lum
    dmstar    = mstar/my_unit_mass
    drstar    = rstar/unit_length
    dtwind    = twind/unit_temperature
    dmdot     = mdot/(my_unit_mass/unit_time)
    dasound   = asound/unit_velocity
    dclight   = const_c/unit_velocity
    dvesc     = vesc/unit_velocity
    dvinf     = dvesc * (alpha/(1.0d0 - alpha))**0.5d0
    dkappae   = kappae * unit_density * unit_length
    dGgrav    = Ggrav * my_unit_ggrav
    dgammae   = dkappae * dlstar/(4.d0*dpi * dGgrav * dmstar * dclight)

  end subroutine make_dimless_vars

end module mod_usr
