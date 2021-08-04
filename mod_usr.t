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
!   ifrc = 2  : finite disk + opacity cut-off (from SO13 dN/dq distribution)
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
!   > implemented exact dv/dr algorithm for non-uniform grids
!   > changed CAK source computation from CGS to unitless + saved unitless to
!     be in line with other unitless output
!
! (February 2021) -- Flo
!   > implemented special outer boundary conditions, some kinks occured in gcak
!     + fdfac due to abrupt change in dv/dr in outer boundary
!     essentially the AMRVAC 'cont' takes zero gradient but constant slope
!     extrapolation alleviates the problem
!   > put effective gravity force in usr_gravity and removed from usr_source
!   > determination radiation timestep (only CAK line force now) in special_dt
!     by using gcak slot of nwextra in w-array
!===============================================================================

module mod_usr

  ! Include a physics module
  use mod_HD

  implicit none

  ! The usual suspects
  real(8), parameter :: msun=1.989d33, lsun=3.827d33, rsun=6.96d10
  real(8), parameter :: Ggrav=6.67d-8, kappae=0.34d0, sigmaSB=5.67e-5

  ! Unit quantities that are handy
  real(8) :: my_unit_ggrav, my_unit_lum, my_unit_mass

  ! Extra input parameters
  real(8) :: mstar, rstar, rhobound, twind, alpha, Qbar, Qmax, beta
  integer :: ifrc

  ! Additionally useful stellar and wind parameters
  real(8) :: lstar, gammae, vesc, mdot, vinf, asound, logg, logge, heff, mumol

  ! Dimensionless variables of relevant variables
  real(8) :: dlstar, dmstar, drstar, drhobound, dtwind, dkappae, dvesc
  real(8) :: dvinf, dmdot, dasound, dclight, dGgrav, dgammae

  ! Extra variables to store in conservative variables array 'w'
  integer :: my_gcak, my_fdf

contains

  !======================================================================
  ! This routine should set user methods, and activate the physics module
  !======================================================================
  subroutine usr_init()

    call set_coordinate_system("spherical")
    call usr_params_read(par_files)

    ! Choose normalisation units: (length,temp,ndens) or (length,vel,ndens)
    ! numberdensity chosen such that unit density becomes boundary density
    unit_length        = rstar                                        ! cm
    unit_temperature   = twind                                        ! K
    unit_numberdensity = rhobound/((1.0d0+4.0d0*He_abundance)*mp_cgs) ! cm^-3

    call HD_activate()

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initial_conditions
    usr_special_bc     => special_bound
    usr_gravity        => effective_gravity
    usr_source         => line_force
    usr_get_dt         => special_dt

    my_gcak = var_set_extravar("gcak", "gcak")
    my_fdf  = var_set_extravar("fdfac", "fdfac")

  end subroutine usr_init

  !========================================================
  ! Read in the usr.par file with the problem specific list
  !========================================================
  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ mstar, rstar, twind, rhobound, alpha, Qbar, Qmax, &
                         beta, ifrc

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    enddo

    ! Scale to cgs units
    mstar = mstar * msun
    rstar = rstar * rsun

  end subroutine usr_params_read

  !====================================================================
  ! Compute some quantities of interest (in CGS) before making unitless
  !====================================================================
  subroutine initglobaldata_usr

    ! Stellar structure
    lstar  = 4.0d0*dpi * rstar**2.0d0 * sigmaSB * twind**4.0d0
    gammae = kappae * lstar/(4.0d0*dpi * Ggrav * mstar * const_c)
    logg   = log10(Ggrav * mstar/rstar**2.0d0)
    logge  = logg + log10(1.0d0 - gammae)
    mumol  = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    asound = sqrt(twind * kB_cgs/(mumol * mp_cgs))
    heff   = asound**2.0d0 / 10.0d0**logge

    ! Wind quantities in CAK theory
    Qmax   = Qmax * Qbar
    vesc   = sqrt(2.0d0 * Ggrav * mstar * (1.0d0 - gammae)/rstar)
    vinf   = vesc * sqrt(alpha/(1.0d0 - alpha))
    mdot   = lstar/const_c**2.0d0 * alpha/(1.0d0 - alpha) &
              * (Qbar * gammae/(1.0d0 - gammae))**((1.0d0 - alpha)/alpha)

    call make_dimless_vars()

  end subroutine initglobaldata_usr

  !==========================================================================
  ! Initial conditions start from spherical 1d beta law and mass conservation
  !==========================================================================
  subroutine initial_conditions(ixI^L,ixO^L,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: sfac

    ! Small offset (asound/vinf) to avoid starting at terminal wind speed
    sfac = 1.0d0 - 1.0d-3**(1.0d0/beta)

    where (x(ixI^S,1) >= drstar)
       w(ixI^S,mom(1)) = dvinf * ( 1.0d0 - sfac * drstar / x(ixI^S,1) )**beta

       w(ixI^S,rho_) = dmdot / (4.0d0*dpi * x(ixI^S,1)**2.0d0 * w(ixI^S,mom(1)))
    endwhere

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

  !==================================================================
  ! Special user boundary conditions at inner + outer radial boundary
  !==================================================================
  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: i

    select case (iB)
    case(1)

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

      call hd_to_conserved(ixI^L,ixI^L,w,x)

    case(2)
      ! Constant extrapolation of all

      call hd_to_primitive(ixI^L,ixI^L,w,x)

      w(ixB^S,rho_) = w(ixBmin1-1,rho_) * (x(ixBmin1-1,1) / x(ixB^S,1))**2.0d0

      do i = ixBmin1,ixBmax1
        w(i,mom(1)) = w(i-1,mom(1)) + (w(ixBmin1-1,mom(1)) - w(ixBmin1-2,mom(1)))
      enddo

      call hd_to_conserved(ixI^L,ixI^L,w,x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  !=======================================================================
  ! Extra source using the analytical CAK line force in Gayley's formalism
  !=======================================================================
  subroutine line_force(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L, iw^LIM
    real(8), intent(in)    :: qdt, qtC, qt
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: vr(ixI^S), rho(ixI^S)
    real(8) :: dvdr_up(ixO^S), dvdr_down(ixO^S), dvdr_cent(ixO^S), dvdr(ixO^S)
    real(8) :: gcak(ixO^S), beta_fd(ixO^S), fdfac(ixO^S)
    real(8) :: taum(ixO^S), taumfac(ixO^S)
    real(8) :: fac, fac1, fac2
    integer :: jx^L, hx^L

    ! Define time-centred, radial velocity from the radial momentum and density
    vr(ixI^S)  = wCT(ixI^S,mom(1)) / wCT(ixI^S,rho_)
    rho(ixI^S) = wCT(ixI^S,rho_)

    ! Index +1 (j) and index -1 (h) in radial direction; kr(dir,dim)=1, dir=dim
    jx^L=ixO^L+kr(1,^D);
    hx^L=ixO^L-kr(1,^D);

    ! Get dv/dr on non-uniform grid according to Sundqvist & Veronis (1970)
    ! Forward difference
    dvdr_up(ixO^S) = (x(ixO^S,1) - x(hx^S,1)) * vr(jx^S) &
                     / ((x(jx^S,1) - x(ixO^S,1)) * (x(jx^S,1) - x(hx^S,1)))

    ! Backward difference
    dvdr_down(ixO^S) = -(x(jx^S,1) - x(ixO^S,1)) * vr(hx^S) &
                        / ((x(ixO^S,1) - x(hx^S,1)) * (x(jx^S,1) - x(hx^S,1)))

    ! Central difference
    dvdr_cent(ixO^S) = (x(jx^S,1) + x(hx^S,1) - 2.0d0*x(ixO^S,1)) * vr(ixO^S) &
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

    ! Calculate CAK line force
    fac1 = 1.0d0/(1.0d0 - alpha) * dkappae * dlstar*Qbar/(4.0d0*dpi * dclight)
    fac2 = 1.0d0/(dclight * Qbar * dkappae)**alpha
    fac  = fac1 * fac2

    gcak(ixO^S) = fac/x(ixO^S,1)**2.0d0 * (dvdr(ixO^S)/rho(ixO^S))**alpha

    ! Based on wind option do corrections
    if (ifrc == 0) then
      gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S)
    elseif (ifrc == 1) then
      gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S)
    elseif (ifrc == 2) then
      taum(ixO^S) = dkappae * dclight * Qmax * rho(ixO^S)/dvdr(ixO^S)

      taumfac(ixO^S) = ((1.0d0 + taum(ixO^S))**(1.0d0 - alpha) - 1.0d0) &
                          / taum(ixO^S)**(1.0d0 - alpha)

      gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S) * taumfac(ixO^S)
    else
      call mpistop('Error in wind option, take a valid ifrc=0,1,2')
    endif

    ! Fill the nwextra slots for output
    w(ixO^S,my_gcak) = gcak(ixO^S)
    w(ixO^S,my_fdf)  = fdfac(ixO^S)

    ! Update conservative vars: w = w + qdt*gsource
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * gcak(ixO^S) * wCT(ixO^S,rho_)

  end subroutine line_force

  !========================================================================
  ! After first iteration the usr_source routine has been called, take now
  ! also timestep from CAK line force into account
  !========================================================================
  subroutine special_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: dx^D, x(ixI^S,1:ndim)
    real(8), intent(in)    :: w(ixI^S,1:nw)
    real(8), intent(inout) :: dtnew

    ! Local variables
    real(8) :: tdum(ixO^S), dt_cak

    ! Get dt from line force that is saved in the w-array in nwextra slot
    tdum(ixO^S) = sqrt( block%dx(ixO^S,1) / abs(w(ixO^S,my_gcak)) )
    dt_cak      = courantpar * minval(tdum(ixO^S))

    if (it >= 1) then
      dtnew = min(dtnew,dt_cak)
    endif

  end subroutine special_dt

  !===================================================================
  ! Combine stellar gravity and continuum electron scattering into an
  ! effective gravity using Eddington's gamma
  !===================================================================
  subroutine effective_gravity(ixI^L,ixO^L,wCT,x,gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim)
    real(8), intent(in)  :: wCT(ixI^S,1:nw)
    real(8), intent(out) :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -dGgrav*dmstar * (1.0d0-dgammae)/x(ixO^S,1)**2.0d0

  end subroutine effective_gravity

  !======================================================================
  ! Normalise relevant quantities to be used in the code + print overview
  !======================================================================
  subroutine make_dimless_vars()

    ! Local variable
    character(len=8) :: todayis

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
    dvinf     = dvesc * sqrt(alpha/(1.0d0 - alpha))
    dkappae   = kappae * unit_density * unit_length
    dGgrav    = Ggrav * my_unit_ggrav
    dgammae   = dkappae * dlstar/(4.d0*dpi * dGgrav * dmstar * dclight)

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
      print*, 'CAK vinf        = ', vinf
      print*, 'FD vinf         = ', 3.0d0 * vesc
      print*
      print*, 'wind option            = ', ifrc
      print*, '   0 : radial stream CAK '
      print*, '   1 : CAK + fd          '
      print*, '   2 : CAK + cut-off     '
      print*, 'surface density        = ', rhobound
      print*, 'analytic Mdot CAK      = ', mdot * (const_years/msun)
      print*, '... with FD correction = ', mdot/(1.0d0 + alpha)**(1.0d0/alpha) * (const_years/msun)
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

  end subroutine make_dimless_vars

end module mod_usr
