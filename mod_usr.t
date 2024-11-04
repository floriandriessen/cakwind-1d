!==============================================================================
! 1D CAK wind module for launching a spherically symmetric line-driven wind
! from the stellar surface
!
! Coded up by Flo for his KU Leuven PhD thesis 2018/2022
!
! Options for wind to be specified in usr.par file:
!   ifrc = 0  : radially streaming CAK wind
!   ifrc = 1  : finite disk corrected CAK wind
!   ifrc = 2  : finite disk + opacity cut-off (Sundqvist & Owocki 2013, dN/dq)
!==============================================================================
module mod_usr

  use mod_hd

  implicit none

  ! Extra input parameters
  integer :: ifrc
  real(8) :: mstar_sol, rstar_sol, twind_cgs, rhosurf_cgs
  real(8) :: alpha, Qbar, Qmax, beta

  ! Dimensionless variables of relevant variables
  real(8) :: lstar, mstar, rstar, rhosurf, twind, kappae, vinf, mdot, gammae
  real(8) :: csound, clight, gmstar

  ! 1-D CAK line force option from ifrc
  integer, parameter :: radstream=0, fdisc=1, fdisc_cutoff=2

  ! Extra variables to store in conservative variables array 'w'
  integer :: igcak_, ifdfac_

contains

  !============================================================================
  ! This routine should set user methods, and activate the physics module.
  !============================================================================
  subroutine usr_init()

    call usr_params_read(par_files)

    ! Choose normalisation units:
    unit_length      = rstar_sol * const_RSun
    unit_temperature = twind_cgs
    unit_density     = rhosurf_cgs

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initial_conditions
    usr_special_bc     => special_bound
    usr_gravity        => effective_gravity
    usr_source         => line_force
    usr_get_dt         => special_dt

    call set_coordinate_system("spherical")
    call hd_activate()

    igcak_  = var_set_extravar("gcak", "gcak")
    ifdfac_ = var_set_extravar("fdfac", "fdfac")

  end subroutine usr_init

  !============================================================================
  ! Read in the usr.par file with a user problem specific list.
  !============================================================================
  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n
    !--------------------------------------------------------------------------

    namelist /star_list/ mstar_sol, rstar_sol, twind_cgs, rhosurf_cgs, &
         alpha, Qbar, Qmax, beta, ifrc

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read

  !============================================================================
  ! Initialise global variables to be used during run. This subroutine is read
  ! at the first iteration and after a restart.
  !============================================================================
  subroutine initglobaldata_usr

    ! Local variables
    character(len=8) :: todayis
    real(8) :: unit_ggrav, unit_lum, unit_mass
    real(8) :: lstar_cgs, mstar_cgs, rstar_cgs, gammae_cgs, vesc_cgs, mdot_cgs
    real(8) :: vinf_cgs, csound_cgs, logg_cgs, logge_cgs, heff_cgs, mumol, vesc
    !--------------------------------------------------------------------------

    mstar_cgs = mstar_sol * const_MSun
    rstar_cgs = rstar_sol * const_RSun

    ! Stellar structure
    lstar_cgs  = 4.0d0*dpi * rstar_cgs**2.0d0 * const_sigma * twind_cgs**4.0d0
    gammae_cgs = const_kappae &
         * lstar_cgs / (4.0d0*dpi * const_G * mstar_cgs * const_c)
    logg_cgs   = log10(const_G * mstar_cgs / rstar_cgs**2.0d0)
    logge_cgs  = logg_cgs + log10(1.0d0 - gammae_cgs)
    mumol      = (1.0d0 + 4.0d0*He_abundance) / (2.0d0 + 3.0d0*He_abundance)
    csound_cgs = sqrt(twind_cgs * kB_cgs / (mumol * mp_cgs))
    heff_cgs   = csound_cgs**2.0d0 / 10.0d0**logge_cgs

    ! Wind quantities in CAK theory
    Qmax     = Qmax * Qbar
    vesc_cgs = sqrt( 2.0d0 * const_G * mstar_cgs * (1.0d0 - gammae_cgs) &
         / rstar_cgs )
    vinf_cgs = vesc_cgs * sqrt(alpha / (1.0d0 - alpha))
    mdot_cgs = lstar_cgs / const_c**2.0d0 * alpha/(1.0d0 - alpha) &
         * (Qbar * gammae_cgs / (1.0d0 - gammae_cgs))**((1.0d0 - alpha)/alpha)

    ! Code units
    unit_ggrav = unit_density * unit_time**2.0d0
    unit_lum   = unit_density * unit_length**5.0d0 / unit_time**3.0d0
    unit_mass  = unit_density * unit_length**3.0d0

    rhosurf = rhosurf_cgs / unit_density
    lstar   = lstar_cgs / unit_lum
    mstar   = mstar_cgs / unit_mass
    rstar   = rstar_cgs / unit_length
    twind   = twind_cgs / unit_temperature
    mdot    = mdot_cgs / (unit_mass / unit_time)
    csound  = csound_cgs / unit_velocity
    clight  = const_c / unit_velocity
    vesc    = vesc_cgs / unit_velocity
    vinf    = vesc * sqrt(alpha / (1.0d0 - alpha))
    kappae  = const_kappae * unit_density * unit_length
    gmstar  = const_G * unit_ggrav * mstar
    gammae  = kappae * lstar / (4.d0*dpi * gmstar * clight)

    if (mype == 0 .and. .not.convert) then
      call date_and_time(todayis)
      print*, 'MPI-AMRVAC simulation ran on ', todayis(7:8), '/', &
           todayis(5:6), '/', todayis(1:4)
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
      print*, 'L/Lsun                 = ', lstar_cgs / const_LSun
      print*, 'M/Msun                 = ', mstar_cgs / const_MSun
      print*, 'R/Rsun                 = ', rstar_cgs / const_RSun
      print*, 'Twind                  = ', twind_cgs
      print*, 'Mean molecular weight  = ', mumol
      print*, 'log(g)                 = ', logg_cgs
      print*, 'eff. log(g)            = ', logge_cgs
      print*, 'eff. scale height heff = ', heff_cgs
      print*, 'heff/Rstar             = ', heff_cgs / rstar_cgs
      print*, 'Eddington gamma        = ', gammae_cgs
      print*
      print*, 'adiabatic gamma = ', hd_gamma
      print*, 'alpha           = ', alpha
      print*, 'Qbar            = ', Qbar
      print*, 'Qmax/Qbar       = ', Qmax/Qbar
      print*, 'beta            = ', beta
      print*, 'asound          = ', csound_cgs
      print*, 'eff. vesc       = ', vesc_cgs
      print*, 'CAK vinf        = ', vinf_cgs
      print*, 'FD vinf         = ', 3.0d0 * vesc_cgs
      print*
      print*, 'wind option            = ', ifrc
      print*, '   0 : radial stream CAK '
      print*, '   1 : CAK + fd          '
      print*, '   2 : CAK + cut-off     '
      print*, 'surface density        = ', rhosurf_cgs
      print*, 'analytic Mdot CAK      = ', mdot_cgs * const_years / const_MSun
      print*, '... with FD correction = ', &
           mdot_cgs / (1.0d0 + alpha)**(1.0d0/alpha) * const_years / const_MSun
      print*
      print*, '========================================'
      print*, '    Dimensionless AMRVAC quantities     '
      print*, '========================================'
      print*, 'Extra computed unit quantities:'
      print*, '   unit Lum  = ', unit_lum
      print*, '   unit Mass = ', unit_mass
      print*, '   unit Grav = ', unit_ggrav
      print*, 'Lstar        = ', lstar
      print*, 'Mstar        = ', mstar
      print*, 'Rstar        = ', rstar
      print*, 'Twind        = ', twind
      print*, 'Edd. gamma   = ', gammae
      print*, 'rhosurface   = ', rhosurf
      print*, 'Mdot         = ', mdot
      print*, 'alpha        = ', alpha
      print*, 'Qbar         = ', Qbar
      print*, 'Qmax         = ', Qmax/Qbar
      print*, 'kappae       = ', kappae
      print*, 'csound       = ', csound
      print*, 'eff. vesc    = ', vesc
      print*, 'vinf         = ', vinf
      print*, 'clight       = ', clight
      print*
    endif

  end subroutine initglobaldata_usr

  !============================================================================
  ! Calculate the initial conditions on the grid for the stellar wind outflow.
  ! After a restart this subroutine will not be accessed.
  !============================================================================
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: sfac
    !--------------------------------------------------------------------------

    ! Small offset (asound/vinf) to avoid starting at terminal wind speed
    sfac = 1.0d0 - 1.0d-3**(1.0d0/beta)

    where (x(ixI^S,1) >= rstar)
       w(ixI^S,mom(1)) = vinf * ( 1.0d0 - sfac * rstar / x(ixI^S,1) )**beta
       w(ixI^S,rho_) = &
            mdot / ( 4.0d0*dpi * x(ixI^S,1)**2.0d0 * w(ixI^S,mom(1)) )
    endwhere

    call hd_to_conserved(ixI^L, ixO^L, w, x)

  end subroutine initial_conditions

  !============================================================================
  ! Non-standard boundary conditions at inner + outer edge boundary.
  !============================================================================
  subroutine special_bound(qt, ixI^L, ixB^L, iB, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: i
    !--------------------------------------------------------------------------

    select case (iB)
    case(1)

      call hd_to_primitive(ixI^L, ixI^L, w, x)

      w(ixB^S,rho_) = rhosurf

      ! Radial velocity field (constant slope extrapolation)
      do i = ixBmax1,ixBmin1,-1
        if (i == ixBmax1) then
          w(i,mom(1)) = w(i+1,mom(1)) - (w(i+2,mom(1)) - w(i+1,mom(1))) &
               * (x(i+1,1) - x(i,1)) / (x(i+2,1) - x(i+1,1))
        else
          w(i,mom(1)) = w(i+1,mom(1))
        endif
      enddo

      ! Prohibit ghosts to be supersonic, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), csound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -csound)

      call hd_to_conserved(ixI^L, ixI^L, w, x)

    case(2)
      ! Constant extrapolation of all

      call hd_to_primitive(ixI^L, ixI^L, w, x)

      w(ixB^S,rho_) = w(ixBmin1-1,rho_) * (x(ixBmin1-1,1) / x(ixB^S,1))**2.0d0

      do i = ixBmin1,ixBmax1
        w(i,mom(1)) = w(i-1,mom(1)) &
             + ( w(ixBmin1-1,mom(1)) - w(ixBmin1-2,mom(1)) )
      enddo

      call hd_to_conserved(ixI^L, ixI^L, w, x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  !============================================================================
  ! Calculate extra source terms for the problem.
  !============================================================================
  subroutine line_force(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L, iw^LIM
    real(8), intent(in)    :: qdt, qtC, qt
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    integer :: jx^L, hx^L
    real(8) :: vr(ixI^S), rho(ixI^S)
    real(8) :: dvdr_up(ixO^S), dvdr_down(ixO^S), dvdr_cent(ixO^S), dvdr(ixO^S)
    real(8) :: ge(ixO^S), gcak(ixO^S), beta_fd(ixO^S), fdfac(ixO^S), &
         tausob(ixO^S)
    !--------------------------------------------------------------------------

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
         * (rstar / x(ixO^S,1))**2.0d0

    ! Check the finite disk array and determine finite disk factor
    select case (ifrc)
    case(radstream)
      fdfac(ixO^S) = 1.0d0

    case(fdisc, fdisc_cutoff)
      where (beta_fd(ixO^S) >= 1.0d0)
        fdfac(ixO^S) = 1.0d0 / (1.0d0 + alpha)
      elsewhere (beta_fd(ixO^S) < -1.0d10)
        fdfac(ixO^S) = abs(beta_fd(ixO^S))**alpha / (1.0d0 + alpha)
      elsewhere (abs(beta_fd(ixO^S)) > 1.0d-3)
        fdfac(ixO^S) = (1.0d0 - (1.0d0 - beta_fd(ixO^S))**(1.0d0 + alpha)) &
             / (beta_fd(ixO^S) * (1.0d0 + alpha))
      elsewhere
        fdfac(ixO^S) = 1.0d0 - 0.5d0*alpha * beta_fd(ixO^S) &
             * (1.0d0 + 1.0d0/3.0d0 * (1.0d0 - alpha)*beta_fd(ixO^S))
      endwhere
    end select

    ! Thomson force
    ge(ixO^S) = kappae * lstar / (4.0d0*dpi * clight * x(ixO^S,1)**2.0d0)

    ! Sobolev optical depth for line ensemble (tau = Qbar * t_r) and CAK force
    select case (ifrc)
    case(radstream, fdisc)
      tausob(ixO^S) = Qbar * kappae * clight * rho(ixO^S) / dvdr(ixO^S)
      gcak(ixO^S) = Qbar/(1.0d0 - alpha) * ge(ixO^S) / tausob(ixO^S)**alpha

    case(fdisc_cutoff)
      tausob(ixO^S) = Qmax * kappae * clight * rho(ixO^S) / dvdr(ixO^S)
      gcak(ixO^S) = Qbar * ge(ixO^S) &
           * ( (1.0d0 + tausob(ixO^S))**(1.0d0 - alpha) - 1.0d0 ) &
           / ( (1.0d0 - alpha) * tausob(ixO^S) )

    case default
      call mpistop("Error in wind option. Take a valid ifrc=0,1,2")
    end select

    gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S)

    ! Fill the nwextra slots for output
    w(ixO^S,igcak_)  = gcak(ixO^S)
    w(ixO^S,ifdfac_) = fdfac(ixO^S)

    ! Update conservative vars: w = w + qdt*gsource
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * gcak(ixO^S) * wCT(ixO^S,rho_)

  end subroutine line_force

  !============================================================================
  ! Adjust the hydrodynamic timestep by taking extra forces into account.
  !============================================================================
  subroutine special_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: dx^D, x(ixI^S,1:ndim)
    real(8), intent(in)    :: w(ixI^S,1:nw)
    real(8), intent(inout) :: dtnew

    ! Local variables
    real(8) :: tdum(ixO^S), dt_cak
    !--------------------------------------------------------------------------

    ! Get dt from line force that is saved in the w-array in nwextra slot
    tdum(ixO^S) = sqrt( block%dx(ixO^S,1) / abs(w(ixO^S,igcak_)) )
    dt_cak      = courantpar * minval(tdum(ixO^S))

    dtnew = min(dtnew, dt_cak)

  end subroutine special_dt

  !============================================================================
  ! Combine stellar gravity and continuum electron scattering into an
  ! effective gravity using Eddington's gamma.
  !============================================================================
  subroutine effective_gravity(ixI^L, ixO^L, wCT, x, gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim)
    real(8), intent(in)  :: wCT(ixI^S,1:nw)
    real(8), intent(out) :: gravity_field(ixI^S,ndim)
    !--------------------------------------------------------------------------

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -gmstar * (1.0d0 - gammae) / x(ixO^S,1)**2.0d0

  end subroutine effective_gravity

end module mod_usr
