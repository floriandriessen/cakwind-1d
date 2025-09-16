!==============================================================================
! 1D CAK wind module for launching a spherically symmetric line-driven wind
! from the stellar surface
!
! Coded up by Flo for his KU Leuven PhD thesis 2018/2022
!
! Options for wind to be specified in usr.par file:
!   ifrc = 0  : radially streaming CAK wind
!   ifrc = 1  : finite disk corrected CAK wind
!   ifrc = 2  : finite disk + opacity cut-off
!==============================================================================
module mod_usr

  use mod_hd

  implicit none

  ! Extra input parameters
  integer :: ifrc
  real(8) :: mstar_sol, rstar_sol, twind_cgs, rhosurf_cgs
  real(8) :: cak_alpha, gayley_qbar, gayley_q0, beta

  ! Dimensionless variables of relevant variables
  real(8) :: lstar, mstar, rstar, rhosurf, twind, kappae, vinf, mdot
  real(8) :: csound, clight, gmstar
  logical :: use_lte_table = .false.

  ! Variables for adaptive surface mass density
  logical :: use_poniatowski_bc = .false.
  real(8) :: timedyn, gammae
  real(8), parameter :: tfloor = 0.8d0, rho_coupling = 2.0d0

  ! 1-D CAK line force option from ifrc
  integer, parameter :: radstream = 0, fdisc = 1, fdisc_cutoff = 2

  ! Extra variables to store in conservative variables array 'w'
  integer :: ige_, igcak_, ifdfac_, ialpha_, iqbar_, iq0_, ike_, ikcak_

contains

  !============================================================================
  ! This routine should set user methods, and activate the physics module.
  !============================================================================
  subroutine usr_init()

    call usr_params_read(par_files)

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initial_conditions
    usr_special_bc     => special_bound
    usr_gravity        => stellar_gravity
    usr_source         => line_force
    usr_get_dt         => special_dt

    call set_coordinate_system("spherical")
    call hd_activate()

    ige_    = var_set_extravar("gelectron", "gelectron")
    igcak_  = var_set_extravar("gcak", "gcak")
    ifdfac_ = var_set_extravar("fdfactor", "fdfactor")
    ialpha_ = var_set_extravar("alpha", "alpha")
    iqbar_  = var_set_extravar("Qbar", "Qbar")
    iq0_    = var_set_extravar("Q0", "Q0")
    ike_    = var_set_extravar("kappae", "kappae")
    ikcak_  = var_set_extravar("kappacak", "kappacak")

  end subroutine usr_init

  !============================================================================
  ! Read in the usr.par file with a user problem specific list.
  !============================================================================
  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n
    !--------------------------------------------------------------------------

    namelist /star_list/ mstar_sol, rstar_sol, twind_cgs, rhosurf_cgs, &
         cak_alpha, gayley_qbar, gayley_q0, beta, ifrc, use_lte_table

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

    use mod_cak_opacity, only: init_cak_table

    ! Local variables
    character(len=8) :: todayis
    real(8) :: unit_ggrav, unit_lum
    real(8) :: lstar_cgs, mstar_cgs, rstar_cgs, vesc_cgs, mdot_cgs
    real(8) :: vinf_cgs, csound_cgs, logg_cgs, logge_cgs, heff_cgs, mumol, vesc
    !--------------------------------------------------------------------------

    mstar_cgs = mstar_sol * const_MSun
    rstar_cgs = rstar_sol * const_RSun

    ! Stellar structure
    lstar_cgs  = 4.0d0*dpi * rstar_cgs**2.0d0 * const_sigma * twind_cgs**4.0d0
    gammae     = const_kappae * lstar_cgs &
         / (4.0d0*dpi * const_G * mstar_cgs * const_c)
    logg_cgs   = log10(const_G * mstar_cgs / rstar_cgs**2.0d0)
    logge_cgs  = logg_cgs + log10(1.0d0 - gammae)
    mumol      = (1.0d0 + 4.0d0*He_abundance) / (2.0d0 + 3.0d0*He_abundance)
    csound_cgs = sqrt(twind_cgs * kB_cgs / (mumol * mp_cgs))
    heff_cgs   = csound_cgs**2.0d0 / 10.0d0**logge_cgs

    ! Wind quantities in CAK theory
    gayley_q0 = gayley_q0 * gayley_qbar
    vesc_cgs  = sqrt( 2.0d0 * const_G * mstar_cgs * (1.0d0 - gammae) &
         / rstar_cgs )
    vinf_cgs  = vesc_cgs * sqrt(cak_alpha / (1.0d0 - cak_alpha))
    mdot_cgs  = lstar_cgs / const_c**2.0d0 * cak_alpha / (1.0d0 - cak_alpha) &
         * (gayley_qbar * gammae / (1.0d0 - gammae))**( (1.0d0 - cak_alpha) &
         / cak_alpha )

    ! Initialise CAK tables
    if (use_lte_table) &
         call init_cak_table("./lte-tables",set_user_tabledir=.true.)

    ! Code units
    unit_length        = rstar_cgs
    unit_temperature   = twind_cgs
    unit_density       = rhosurf_cgs
    unit_numberdensity = unit_density / (mumol * mp_cgs)
    unit_pressure      = unit_numberdensity * kB_cgs * unit_temperature
    unit_velocity      = sqrt(unit_pressure / unit_density)
    unit_time          = unit_length / unit_velocity
    unit_mass          = unit_density * unit_length**3.0d0
    unit_opacity       = unit_length**2.0d0 / unit_mass

    unit_ggrav   = unit_density * unit_time**2.0d0
    unit_lum     = unit_density * unit_length**5.0d0 / unit_time**3.0d0

    rhosurf = rhosurf_cgs / unit_density
    lstar   = lstar_cgs / unit_lum
    mstar   = mstar_cgs / unit_mass
    rstar   = rstar_cgs / unit_length
    twind   = twind_cgs / unit_temperature
    mdot    = mdot_cgs / (unit_mass / unit_time)
    csound  = csound_cgs / unit_velocity
    clight  = const_c / unit_velocity
    vesc    = vesc_cgs / unit_velocity
    vinf    = vesc * sqrt(cak_alpha / (1.0d0 - cak_alpha))
    kappae  = const_kappae / unit_opacity
    gmstar  = const_G * unit_ggrav * mstar
    timedyn = rstar / (3.0d0 * vinf)

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
      print*, 'unit mass          = ', unit_mass
      print*, 'unit opacity       = ', unit_opacity
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
      print*, 'Eddington gamma        = ', gammae
      print*
      print*, 'adiabatic gamma = ', hd_gamma
      print*, 'alpha           = ', cak_alpha
      print*, 'Qbar            = ', gayley_qbar
      print*, 'Q0/Qbar         = ', gayley_q0 / gayley_qbar
      print*, 'beta            = ', beta
      print*, 'asound          = ', csound_cgs
      print*, 'eff. vesc       = ', vesc_cgs
      print*, 'CAK vinf        = ', vinf_cgs
      print*, 'FD vinf         = ', 3.0d0 * vesc_cgs
      print*, 'use_lte_table   = ', use_lte_table
      print*, 'use Luka BC     = ', use_poniatowski_bc
      print*
      print*, 'wind option            = ', ifrc
      print*, '   0 : radial stream CAK '
      print*, '   1 : CAK + fd          '
      print*, '   2 : CAK + cut-off     '
      print*, 'surface density        = ', rhosurf_cgs
      print*, 'analytic Mdot CAK      = ', mdot_cgs * const_years / const_MSun
      print*, '... with FD correction = ', &
           mdot_cgs / (1.0d0 + cak_alpha)**(1.0d0/cak_alpha) &
           * const_years / const_MSun
      print*
      print*, '========================================'
      print*, '    Dimensionless AMRVAC quantities     '
      print*, '========================================'
      print*, 'Extra computed unit quantities:'
      print*, '   unit Lum     = ', unit_lum
      print*, '   unit Mass    = ', unit_mass
      print*, '   unit Grav    = ', unit_ggrav
      print*, '   unit opacity = ', unit_opacity
      print*, 'Lstar        = ', lstar
      print*, 'Mstar        = ', mstar
      print*, 'Rstar        = ', rstar
      print*, 'Twind        = ', twind
      print*, 'rhosurface   = ', rhosurf
      print*, 'Mdot         = ', mdot
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
    real(8) :: sfac, dvdr(ixO^S), ge(ixO^S), gcak(ixO^S), tausob(ixO^S), vsurf
    !--------------------------------------------------------------------------

    ! Fix base speed at 10km/s
    vsurf = 10.0d0 * 1e5 / unit_velocity
    sfac  = 1.0d0 - (vsurf / vinf)**(1.0d0/beta)

    where (x(ixI^S,1) >= rstar)
       w(ixI^S,mom(1)) = vinf * ( 1.0d0 - sfac * rstar / x(ixI^S,1) )**beta
       w(ixI^S,rho_)   = rhosurf * vsurf / w(ixI^S,mom(1)) &
            * (rstar / x(ixI^S,1))**2.0d0
    endwhere

    call hd_to_conserved(ixI^L, ixO^L, w, x)

    ! Initial forces
    ge(ixO^S) = kappae * lstar/(4.0d0*dpi * clight * x(ixO^S,1)**2.0d0)

    dvdr(ixO^S)   = beta * vinf * sfac * rstar / x(ixO^S,1)**2.0d0 &
         * (1.0d0 - sfac * rstar / x(ixO^S,1))**(beta - 1.0d0)
    tausob(ixO^S) = gayley_qbar * kappae * clight * w(ixO^S,rho_) / dvdr(ixO^S)
    gcak(ixO^S)   = gayley_qbar / (1.0d0 - cak_alpha) * ge(ixO^S) &
           / tausob(ixO^S)**cak_alpha

    ! Constant line-statistic parameters at start
    w(ixO^S,ige_)    = ge(ixO^S)
    w(ixO^S,igcak_)  = gcak(ixO^S)
    w(ixO^S,ialpha_) = cak_alpha
    w(ixO^S,iqbar_)  = gayley_qbar
    w(ixO^S,iq0_)    = gayley_q0
    w(ixO^S,ike_)    = kappae
    w(ixO^S,ikcak_)  = gcak(ixO^S) * (4.0d0 * dpi * clight) / lstar

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
    integer :: ir
    real(8) :: wlocal(ixI^S,1:nw), scaleheight, soundspeed(ixI^S)
    !--------------------------------------------------------------------------

    select case (iB)
    case(1)

      call hd_to_primitive(ixI^L, ixI^L, w, x)

      if (use_poniatowski_bc) then
        ! Adaptive lower boundary mass density
        wlocal(ixI^S,1:nw) = w(ixI^S,1:nw)
        call adaptive_surface_density(ixI^L, qt, wlocal, rhosurf, soundspeed)

        ! Compute barometric scale height to get hydrostatic atmosphere
        scaleheight = soundspeed(ixBmax1+1)**2.0d0 &
            * x(ixBmax1+1,1) / (gmstar * (1.0d0 - gammae))

        w(ixB^S,rho_) = rhosurf * exp( -2.0d0*x(ixBmax1+1,1) / scaleheight &
            * (1.0d0 - 1.0d0 / x(ixB^S,1)) )

        ! Radial velocity field (from continuity)
        do ir = ixBmax1,ixBmin1,-1
          w(ir,mom(1)) = w(ir+1,mom(1)) * w(ir+1,rho_) * x(ir+1,1)**2.0d0 &
              / (w(ir,rho_) * x(ir,1)**2.0d0)
        enddo
      else
        ! Standard conditions
        w(ixB^S,rho_) = rhosurf

        ! Radial velocity field (constant slope extrapolation: d^vr/dr^2 = 0)
        do ir = ixBmax1,ixBmin1,-1
          w(ir^%1ixB^S,mom(1)) = 2.0d0*w(ir+1^%1ixB^S,mom(1)) &
              - w(ir+2^%1ixB^S,mom(1))
        enddo
      endif

      ! Prohibit ghosts to be supersonic, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), csound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -csound)

      call hd_to_conserved(ixI^L, ixI^L, w, x)

    case(2)
      ! Constant extrapolation of all

      call hd_to_primitive(ixI^L, ixI^L, w, x)

      w(ixB^S,rho_) = w(ixBmin1-1,rho_) * (x(ixBmin1-1,1) / x(ixB^S,1))**2.0d0

      do ir = ixBmin1,ixBmax1
        w(ir,mom(1)) = w(ir-1,mom(1)) &
             + ( w(ixBmin1-1,mom(1)) - w(ixBmin1-2,mom(1)) )
      enddo

      call hd_to_conserved(ixI^L, ixI^L, w, x)

    case default
      call mpistop("BC not specified")
    end select

  contains

    !==========================================================================
    ! Adaptive boundary condition for mass density at stellar surface.
    !==========================================================================
    subroutine adaptive_surface_density(ixI^L, qt, w, rhobound, soundspeed)

      ! Subroutine arguments
      integer, intent(in)    :: ixI^L
      real(8), intent(in)    :: qt, w(ixI^S,1:nw)
      real(8), intent(inout) :: rhobound
      real(8), intent(out)   :: soundspeed(ixI^S)

      ! Local variables
      integer :: id
      real(8) :: timecurr, timeprev

      integer, save :: nupdate
      real(8), save :: rhosurfav
      !------------------------------------------------------------------------

      ! Dimensionless isothermal sound speed
      soundspeed(ixI^S) = sqrt(tfloor * twind)

      ! Adapt after three dynamical timescales (set by finite-disk CAK model)
      if (qt > 3.0d0 * timedyn) then

        ! Get index of velocity closed to sonic velocity
        id = minloc(abs(w(ixI^S,mom(1)) - soundspeed(ixI^S)), 1)

        ! Compute average density (need to normalise/unnormalise each time)
        timecurr  = qt - timedyn
        timeprev  = qt - dt - timedyn
        rhosurfav = rhosurfav * timeprev
        rhosurfav = rhosurfav + dt * w(id,rho_)
        rhosurfav = rhosurfav / timecurr

        if (qt > nupdate * timedyn) then
          if ( abs(rhobound/(rho_coupling * rhosurfav) - 1) >= 0.2d0) then
            ! print*, "updated boundary density:", &
            !      log10(rhobound * unit_density), &
            !      " to:", log10(rho_coupling * rhosurfav * unit_density)
            rhobound = rho_coupling * rhosurfav
          endif

          nupdate = nupdate + 5
        endif

      else
        rhosurfav = 0.0d0
        nupdate   = 8

        if (it < 10) then
          id = minloc(abs(w(ixI^S,mom(1)) - soundspeed(ixI^S)), 1)
          rhobound  = rho_coupling * w(id,rho_)
        endif
      endif

      ! if (mod(it,10000)==0) print*, 'rhobound', log10(rhobound * unit_density)

    end subroutine adaptive_surface_density

  end subroutine special_bound

  !============================================================================
  ! Calculate extra source terms for the problem.
  !============================================================================
  subroutine line_force(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x)

    use mod_cak_opacity, only: set_cak_opacity

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L, iw^LIM
    real(8), intent(in)    :: qdt, qtC, qt
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    integer :: i, jx^L, hx^L
    real(8) :: vr(ixI^S), rho(ixI^S)
    real(8) :: dvdr_up(ixO^S), dvdr_down(ixO^S), dvdr_cent(ixO^S), dvdr(ixO^S)
    real(8) :: ge(ixO^S), gcak(ixO^S), beta_fd(ixO^S), fdfac(ixO^S), &
         tausob(ixO^S)
    real(8) :: qbar(ixO^S), q0(ixO^S), alpha(ixO^S), kappae(ixO^S)
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

    ! Retrieve line-statistic parameters from local density and temperature
    if (use_lte_table) then
      do i = ixOmin1, ixOmax1
        call set_cak_opacity(rho(i) * unit_density,       &
             tfloor * twind * unit_temperature, alpha(i), &
             qbar(i), q0(i), kappae(i))
      enddo

      ! Make table kappae unitless
      kappae(ixO^S) = kappae(ixO^S) / unit_opacity

    else
      alpha(ixO^S)  = w(ixO^S,ialpha_)
      qbar(ixO^S)   = w(ixO^S,iqbar_)
      q0(ixO^S)     = w(ixO^S,iq0_)
      kappae(ixO^S) = w(ixO^S,ike_)
    endif

    ! Finite disk factor parameterisation (Owocki & Puls 1996)
    beta_fd(ixO^S) = (1.0d0 - vr(ixO^S) / (x(ixO^S,1) * dvdr(ixO^S))) &
         * (rstar / x(ixO^S,1))**2.0d0

    ! Check the finite disk array and determine finite disk factor
    select case (ifrc)
    case(radstream)
      fdfac(ixO^S) = 1.0d0

    case(fdisc, fdisc_cutoff)
      where (beta_fd(ixO^S) >= 1.0d0)
        fdfac(ixO^S) = 1.0d0 / (1.0d0 + alpha(ixO^S))
      elsewhere (beta_fd(ixO^S) < -1.0d10)
        fdfac(ixO^S) = abs(beta_fd(ixO^S))**alpha(ixO^S) &
             / (1.0d0 + alpha(ixO^S))
      elsewhere (abs(beta_fd(ixO^S)) > 1.0d-3)
        fdfac(ixO^S) = &
             (1.0d0 - (1.0d0 - beta_fd(ixO^S))**(1.0d0 + alpha(ixO^S))) &
             / (beta_fd(ixO^S) * (1.0d0 + alpha(ixO^S)))
      elsewhere
        fdfac(ixO^S) = 1.0d0 - 0.5d0*alpha(ixO^S) * beta_fd(ixO^S) &
             * (1.0d0 + (1.0d0 - alpha(ixO^S))/3.0d0 * beta_fd(ixO^S))
      endwhere
    end select

    ! Thomson force
    ge(ixO^S) = kappae(ixO^S) * lstar/(4.0d0*dpi * clight * x(ixO^S,1)**2.0d0)

    ! Sobolev optical depth for line ensemble (tau = Qbar * t_r) and CAK force
    select case (ifrc)
    case(radstream, fdisc)
      tausob(ixO^S) = qbar(ixO^S) * kappae(ixO^S) * clight * rho(ixO^S) &
           / dvdr(ixO^S)
      gcak(ixO^S) = qbar(ixO^S)/(1.0d0 - alpha(ixO^S)) * ge(ixO^S) &
           / tausob(ixO^S)**alpha(ixO^S)

    case(fdisc_cutoff)
      tausob(ixO^S) = q0(ixO^S) * kappae(ixO^S) * clight * rho(ixO^S) &
           / dvdr(ixO^S)
      gcak(ixO^S) = qbar(ixO^S) * ge(ixO^S) &
           * ( (1.0d0 + tausob(ixO^S))**(1.0d0 - alpha(ixO^S)) - 1.0d0 ) &
           / ( (1.0d0 - alpha(ixO^S)) * tausob(ixO^S) )

    case default
      call mpistop("Error in wind option. Take a valid ifrc=0,1,2")
    end select

    gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S)

    ! Fill the nwextra slots for output
    w(ixO^S,ige_)    = ge(ixO^S)
    w(ixO^S,igcak_)  = gcak(ixO^S)
    w(ixO^S,ifdfac_) = fdfac(ixO^S)
    w(ixO^S,ialpha_) = alpha(ixO^S)
    w(ixO^S,iqbar_)  = qbar(ixO^S)
    w(ixO^S,iq0_)    = q0(ixO^S)
    w(ixO^S,ike_)    = kappae(ixO^S)
    w(ixO^S,ikcak_)  = gcak(ixO^S) * (4.0d0*dpi * clight) / lstar

    ! Update conservative vars: w = w + qdt*gsource
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) &
         + qdt * (ge(ixO^S) + gcak(ixO^S)) * wCT(ixO^S,rho_)

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
    real(8) :: tmp(ixO^S), dt_cak
    !--------------------------------------------------------------------------

    ! Get dt from line force that is saved in the w-array in nwextra slot
    tmp(ixO^S) = sqrt(block%dx(ixO^S,1) / abs(w(ixO^S,ige_) + w(ixO^S,igcak_)))
    dt_cak     = courantpar * minval(tmp(ixO^S))

    dtnew = min(dtnew, dt_cak)

  end subroutine special_dt

  !============================================================================
  ! Set the gravitational field of the problem.
  !============================================================================
  subroutine stellar_gravity(ixI^L, ixO^L, wCT, x, gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim)
    real(8), intent(in)  :: wCT(ixI^S,1:nw)
    real(8), intent(out) :: gravity_field(ixI^S,ndim)
    !--------------------------------------------------------------------------

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -gmstar / x(ixO^S,1)**2.0d0

  end subroutine stellar_gravity

end module mod_usr
