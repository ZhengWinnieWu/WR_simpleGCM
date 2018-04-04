
module hs_forcing_mod

!-----------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use     constants_mod, only: KAPPA, CP_AIR, GRAV, PI, SECONDS_PER_DAY, RDGAS

use           fms_mod, only: error_mesg, FATAL, file_exist,       &
                             open_namelist_file, set_domain,      &
			     read_data, check_nml_error,          &
                             mpp_pe, mpp_root_pe, close_file,     &
                             write_version_number, stdlog,        &
                             uppercase,&  !pjk
                             mpp_clock_id,mpp_clock_begin,mpp_clock_end,CLOCK_COMPONENT!,mpp_chksum

use  time_manager_mod, only: time_type, get_time

use  diag_manager_mod, only: register_diag_field, send_data

use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_number_tracers
use   interpolator_mod, only: interpolate_type, interpolator_init, &
                              interpolator, interpolator_end, &
                              CONSTANT, INTERP_WEIGHTED_P, INTERP_LINEAR_P
implicit none
private

!-----------------------------------------------------------------------
!---------- interfaces ------------

   public :: hs_forcing, hs_forcing_init, hs_forcing_end

   type(interpolate_type),save         ::  heating_source_interp
   type(interpolate_type),save         ::  u_interp
   type(interpolate_type),save         ::  v_interp
   type(interpolate_type),save         ::  temp_interp
   type(interpolate_type),save         ::  tau_interp

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

   logical :: no_forcing = .false.
   logical :: surface_forcing_input = .false.

   real :: t_zero=315., t_strat=200., delh=60., delv=10., sigma_b=0.7

   real :: P00 = 1.e5
   real :: epsS=0., epsN=0. !mj eps differ in SH and NH

   real :: ka = -40. !  negative values are damping time in days
   real :: k_strat = -40. ! tjr
   real :: ks =  -4., kf = -1.

   logical :: do_conserve_energy = .true.

   real :: trflux = 1.e-5   !  surface flux for optional tracer
   real :: trsink = -0.     !  damping time for tracer

!------------------- local heating ------------------------------------------------
   character(len=256) :: local_heating_option='' ! Valid options are 'from_file', 'Isidoro', and 'Gaussian'. Local heating not done otherwise.
   character(len=256) :: local_heating_file=''   ! Name of file relative to $work/INPUT  Used only when local_heating_option='from_file'
   real :: local_heating_srfamp=0.0              ! Degrees per day.   Used only when local_heating_option='Isidoro' or 'Gaussian'
   real :: local_heating_constamp=0.0            ! sigma height       Used only when local_heating_option='Gaussian' !mj NOT IMPLEMENTED YET
   real :: local_heating_xwidth=10.              ! degrees longitude  Used only when local_heating_option='Isidoro'
   real :: local_heating_ywidth=10.              ! degrees latitude   Used only when local_heating_option='Isidoro'
   real :: local_heating_xcenter=180.            ! degrees longitude  Used only when local_heating_option='Isidoro'
   real :: local_heating_ycenter=45.             ! degrees latitude   Used only when local_heating_option='Isidoro'
   real :: local_heating_latwidth=0.4            ! radians latitude   Used only when local_heating_option='Gaussian'
   real :: local_heating_sigwidth=0.11           ! sigma height       Used only when local_heating_option='Gaussian'
   real :: local_heating_sigcenter=0.3           ! sigma height       Used only when local_heating_option='Gaussian'
   logical :: polar_heating_option=.false.       ! want to add some heating over the pole? 
   real :: polar_heating_srfamp=0.0              ! Degrees per day    Used only when polar heating_option='true'
   real :: polar_heating_latwidth=0.0            ! radians latitude   Used only when polar_heating_option='true'
   real :: polar_heating_latcenter=0.0           ! radians latitude   Used only when polar_heating_option='true'
   real :: polar_heating_sigwidth=0.0            ! sigma height       Used only when polar_heating_option='true'
   real :: polar_heating_sigcenter=0.0           ! sigma height       Used only when polar_heating_option='true'
   real :: eichelberger_height=12.e3             ! vertical height    Used only when local_heating_option='Eichelberger'
   real :: local_heating_vert_decay=1.e4         ! pascals            Used only when local_heating_option='Isidoro'
!-----------------------------------------------------------------------------------

   logical :: relax_to_specified_wind = .false.
   character(len=256) :: u_wind_file='u', v_wind_file='v' ! Name of files relative to $work/INPUT  Used only when relax_to_specified_wind=.true.

   character(len=256) :: equilibrium_t_option = 'Held_Suarez' ! mj 'Held_Suarez', 'JFV', 'from_file', 'strat_file'
   character(len=256) :: equilibrium_t_file='temp'  ! Name of file relative to $work/INPUT  Used only when equilibrium_t_option='from_file'
   character(len=256) :: equilibrium_tau_option='Held_Suarez' !mj: 'Held-Suarez','JFV','from_file', 'strat_file'
   character(len=256) :: equilibrium_tau_file='tau'  ! Name of file relative to $work/INPUT  Used only when equilibrium_tau_option='from_file' or 'strat_file'

!
!  standard atmosphere (sat), polar vortex (pv) and TOA sponge layer tjr
!
   logical :: pv_sat_flag   = .true.  ! flag for SAT strat + polar vortex
   logical :: sat_only_flag = .false. ! mj flag for SAT strat without  polar vortex
   real :: pv_phi0  = 60.    !  polar vortex edge location (in degrees)
   real :: pv_dphi  = 10.    !  polar vortex edge width (in degrees)
   real :: pv_gamma = -1.e-3 !  polar vortex lapse rate (in degK/m)
!
   logical :: sponge_flag = .true. !flag for sponge at top of model
   real :: sponge_pbottom = 1.e2    !bottom of sponge layer, where damping is zero (Pa)
   real :: sponge_tau_days  = 1.0   !damping time scale for the sponge (days)

   real :: p_tropopause = 0.1       !tropopause pressure divided by reference pressure

   real :: scaife_damp  = 2.5      !damping time scale for the scaife damp (days)
   logical :: scaife_flag = .false. !flag for scaife damping

   real :: sigma_strat1=-1,sigma_strat2=-2 ! levels defining transition from tka to tk_strat tjr

   logical :: sc_flag = .false. ! flag for seasonal cycle in Te through eps and PK02 stratosphere
   real ::  sc_phi0n,sc_phi0s,sc_dphin,sc_dphis !generalizations, for seasonal cycle option, of pv_phi0 and pv_dphi above

!-------------- Jucker et al JGR (2014) background ---------------------
   real :: p_hs=250.e2,p_bd=100.e2   ! boundaries for transition from Held_Suarez to JFV stratosphere
   real :: A_NH_0=15, A_NH_1=45, A_SH_0=25, A_SH_1=60, A_s=15, phi_N=80, phi_S=-80 ! stratospheric polar temperature amplitdues
   real :: T_fact=1.0 ! scale the stratospheric Te profile (0 -> Te=t_strat, 1 -> full Te profile)
   real :: tau_t=40, tau_N_p=20, tau_S_p=20, delta_phi=30, tau_m=5 ! stratospheric relaxation time setup (tropics, northpole, southpole)
   real :: tau_fact=1.0 ! scale the stratospheric tau profile (0 -> tau=tau_t, 1 -> full tau profile)
   logical :: do_seasonal_cycle=.true.  !add seasonal cycle in stratosphere?
   real :: days_per_year=365             !for determining seasonal cycle
!-----------------------------------------------------------------------

   namelist /hs_forcing_nml/  no_forcing, surface_forcing_input,             &
   	    		      t_zero, t_strat, delh, delv,                   &
                              sigma_b, ka, ks, kf, do_conserve_energy,       &
                              trflux, trsink, local_heating_srfamp,          &
                              local_heating_xwidth,  local_heating_ywidth,   &
                              local_heating_xcenter, local_heating_ycenter,  &
                              local_heating_latwidth, local_heating_sigwidth,& !cc
                              polar_heating_option, polar_heating_srfamp,    & !cc
                              polar_heating_sigcenter,polar_heating_latwidth,& !cc
                              polar_heating_latcenter,polar_heating_sigwidth,& !cc
                              local_heating_sigcenter, eichelberger_height,  & !mj !cc
                              local_heating_vert_decay, local_heating_option,&
                              local_heating_file, relax_to_specified_wind,   &
                              u_wind_file, v_wind_file, equilibrium_t_option,&
                              equilibrium_t_file,                            &
                              pv_sat_flag, pv_phi0, pv_dphi, pv_gamma,       &  ! tjr
                              sponge_flag,sponge_pbottom,sponge_tau_days,    &  ! tjr
                              p_tropopause,scaife_damp, scaife_flag,         &  ! hmchen
                              sigma_strat1,sigma_strat2,k_strat,             &  ! tjr
                              sc_flag,sc_phi0n,sc_phi0s,sc_dphin,sc_dphis,   &  !mj
                              epsS,epsN,                                     &  !mj
                              sat_only_flag,                                 &  !mj
                              equilibrium_tau_option,equilibrium_tau_file,   &  !mj
                              p_hs,p_bd,A_NH_0,A_NH_1,A_SH_0,A_SH_1,A_s,     &  !mj
                              phi_N,phi_S,tau_t,tau_N_p,tau_S_p,delta_phi,   &  !mj
                              Te_fact,tau_fact,                              &  !mj
                              tau_m,do_seasonal_cycle,days_per_year             !mj

!-----------------------------------------------------------------------

   character(len=128) :: version='$Id: hs_forcing.f90, 2014/10/24 mj $'
   character(len=128) :: tagname='$Name: riga_201012 $'

   real :: tka, tks, vkf
   real :: scdamp
   real :: trdamp, twopi
   real :: tk_strat

   integer :: id_teq, id_tau, id_tdt, id_udt, id_vdt, id_tdt_diss, id_diss_heat, id_local_heating, id_newtonian_damping
   real    :: missing_value = -1.e10
   real    :: xwidth, ywidth, xcenter, ycenter ! namelist values converted from degrees to radians
   real    :: srfamp, polar_srfamp! local_heating_srfamp converted from deg/day to deg/sec
   character(len=14) :: mod_name = 'hs_forcing'

   logical :: module_is_initialized = .false.
   integer :: id_newt_damp1,id_newt_damp2,id_newt_damp3
!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine hs_forcing ( is, ie, js, je, dt, Time, lon, lat, p_half, p_full, &
                         u, v, t, vor, r, um, vm, tm, vom, rm, udt, vdt, tdt, rdt,&
                         s_geo, & !mj
                         mask, kbot )
!-----------------------------------------------------------------------
use     tracer_manager_mod, only: get_tracer_index, NO_TRACER !mj

!-----------------------------------------------------------------------
   integer, intent(in)                        :: is, ie, js, je
      real, intent(in)                        :: dt
 type(time_type), intent(in)                  :: Time
      real, intent(in),    dimension(:,:)     :: lon, lat
      real, intent(in),    dimension(:,:)     :: s_geo !mj
      real, intent(in),    dimension(:,:,:)   :: p_half, p_full
      real, intent(in),    dimension(:,:,:)   :: u, v, t, vor, um, vm, tm, vom
      real, intent(in),    dimension(:,:,:,:) :: r, rm
      real, intent(inout), dimension(:,:,:)   :: udt, vdt, tdt
      real, intent(inout), dimension(:,:,:,:) :: rdt

      real, intent(in),    dimension(:,:,:), optional :: mask
   integer, intent(in),    dimension(:,:)  , optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(t,1),size(t,2))           :: ps, diss_heat, surface_forcing
   real, dimension(size(t,1),size(t,2),size(t,3)) :: ttnd, utnd, vtnd, teq, tau, pmass
   real, dimension(size(r,1),size(r,2),size(r,3)) :: rst, rtnd, vort
   real, dimension(size(r,1),size(r,2),size(r,3)) :: tst, pvt !mj
   integer :: i, j, k, kb, n, num_tracers
   integer :: n_hum, n_pv !mj
   logical :: used
   real    :: flux, sink, value
   character(len=128) :: scheme, params

!-----------------------------------------------------------------------
     if (no_forcing) return

     if (.not.module_is_initialized) call error_mesg ('hs_forcing','hs_forcing_init has not been called', FATAL)

!-----------------------------------------------------------------------
!     surface pressure

     if (present(kbot)) then
         do j=1,size(p_half,2)
         do i=1,size(p_half,1)
            kb = kbot(i,j)
            ps(i,j) = p_half(i,j,kb+1)
         enddo
         enddo
     else
            ps(:,:) = p_half(:,:,size(p_half,3))
     endif

!-----------------------------------------------------------------------
!     rayleigh damping of wind components near the surface

      call rayleigh_damping ( Time, ps, p_full, p_half, u, v, utnd, vtnd, mask=mask )

      if (do_conserve_energy) then
         ttnd = -((um+.5*utnd*dt)*utnd + (vm+.5*vtnd*dt)*vtnd)/CP_AIR
         tdt = tdt + ttnd
         if (id_tdt_diss > 0) used = send_data ( id_tdt_diss, ttnd, Time, is, js)
       ! vertical integral of ke dissipation
         if ( id_diss_heat > 0 ) then
          do k = 1, size(t,3)
            pmass(:,:,k) = p_half(:,:,k+1)-p_half(:,:,k)
          enddo
          diss_heat = CP_AIR/GRAV * sum( ttnd*pmass, 3)
          used = send_data ( id_diss_heat, diss_heat, Time, is, js)
         endif
      endif

      udt = udt + utnd
      vdt = vdt + vtnd

      if (id_udt > 0) used = send_data ( id_udt, utnd, Time, is, js)
      if (id_vdt > 0) used = send_data ( id_vdt, vtnd, Time, is, js)

!-----------------------------------------------------------------------
!     thermal forcing for held & suarez (1994) benchmark calculation
      call newtonian_damping ( Time, lat, ps, p_full, p_half, t, ttnd, teq, tau, mask,surface_forcing )

      tdt = tdt + ttnd
      if (id_newtonian_damping > 0) used = send_data(id_newtonian_damping, ttnd, Time, is, js)

      if(trim(local_heating_option) /= '') then
        call local_heating ( Time, is, js, lon, lat, ps, p_full, p_half, s_geo, t, ttnd )
        tdt = tdt + ttnd
        if (id_local_heating > 0) used = send_data ( id_local_heating, ttnd, Time, is, js)
      endif

      if (id_tdt > 0) used = send_data ( id_tdt, tdt, Time, is, js)
      if (id_teq > 0) used = send_data ( id_teq, teq, Time, is, js)
      if (id_tau > 0) used = send_data ( id_tau, tau, Time, is, js)

!-----------------------------------------------------------------------
!     -------- tracers -------

      call get_number_tracers(MODEL_ATMOS, num_tracers=num_tracers)
      n_hum   = get_tracer_index(MODEL_ATMOS, 'sphum') !mj
      n_pv   = get_tracer_index(MODEL_ATMOS, 'pv') !mj
      
      if(num_tracers == size(rdt,4)) then
        do n = 1, size(rdt,4)
           flux = trflux
           sink = trsink
           if (n == n_hum) then !mj
              rst = rm(:,:,:,n) + dt*rdt(:,:,:,n) !mj
              tst = tm + dt*tdt !mj
              call sphum_source_sink ( flux, sink, p_full, rst, rtnd, s_geo, tst ,dt, kbot ) !mj
              rdt(:,:,:,n) = rdt(:,:,:,n) + rtnd !mj
           elseif(n == get_tracer_index(MODEL_ATMOS,'age_of_air')) then !mj age of air tracer
              rst = rm(:,:,:,n) + dt*rdt(:,:,:,n) !mj
              call age_source_sink ( flux, sink, p_half, rst, rtnd, kbot ) !mj
              rdt(:,:,:,n) = rdt(:,:,:,n) + rtnd !mj
           elseif(n == get_tracer_index(MODEL_ATMOS,'methane')) then !mj methane tracer
              rst = rm(:,:,:,n) + dt*rdt(:,:,:,n) !mj
              call methane_source_sink ( flux, sink, p_half, rst, rtnd, kbot ) !mj
              rdt(:,:,:,n) = rdt(:,:,:,n) + rtnd !mj      
           elseif(n == n_pv) then
              rst = rm(:,:,:,n) + dt*rdt(:,:,:,n) !mj
              tst = tm + dt*tdt !mj
              call pv_tracer(vom, tst, lat, p_full, dt, rst, rtnd)
              rdt(:,:,:,n) = rdt(:,:,:,n) + rtnd !mj   
           elseif(n == get_tracer_index(MODEL_ATMOS, 'APV')) then !mj Blocking index as Schwierz2004
              rst = rm(:,:,:,n) + dt*rdt(:,:,:,n) !mj
              tst = tm + dt*tdt !mj
              if(n_pv < 1 .or. n_pv > n )then !pv tracer not present or after APV in field_table
                 call pv_tracer(vom, tst, lat, p_full, 0.0, rst, pvt) ! pvt = pv, not dpv/dt
                 call apv_tracer(pvt,p_half,dt,rst,rtnd)
              else
                 pvt = rm(:,:,:,n_pv) + dt*rdt(:,:,:,n_pv) !mj
                 call apv_tracer(pvt,p_full,dt,rst,rtnd)
              endif
              rdt(:,:,:,n) = rdt(:,:,:,n) + rtnd !mj  
           else !mj
              if (query_method('tracer_sms', MODEL_ATMOS, n, scheme, params)) then
                 if (uppercase(trim(scheme)) == 'NONE') cycle
                 if (uppercase(trim(scheme)) == 'OFF') then
                    flux = 0.; sink = 0.
                 else
                    if (parse(params,'flux',value) == 1) flux = value
                    if (parse(params,'sink',value) == 1) sink = value
                 endif
              endif
              rst = rm(:,:,:,n) + dt*rdt(:,:,:,n)
              call tracer_source_sink ( flux, sink, p_half, rst, rtnd, kbot )
              rdt(:,:,:,n) = rdt(:,:,:,n) + rtnd
           endif !mj
        enddo
      else
        call error_mesg('hs_forcing','size(rdt,4) not equal to num_tracers', FATAL)
      endif

!-----------------------------------------------------------------------

 end subroutine hs_forcing

!#######################################################################

 subroutine hs_forcing_init ( axes, Time, lonb, latb )

!-----------------------------------------------------------------------
!
!           routine for initializing the model with an
!              initial condition at rest (u & v = 0)
!
!-----------------------------------------------------------------------

           integer, intent(in) :: axes(4)
   type(time_type), intent(in) :: Time
   real, intent(in), optional, dimension(:,:) :: lonb, latb
   

!-----------------------------------------------------------------------
   integer  unit, io, ierr

!     ----- read namelist -----

#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=hs_forcing_nml, iostat=io)
     ierr = check_nml_error(io, 'hs_forcing_nml')
#else
      if (file_exist('input.nml')) then
         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=hs_forcing_nml, iostat=io, end=10)
            ierr = check_nml_error (io, 'hs_forcing_nml')
         enddo
  10     call close_file (unit)
      endif
#endif
!mj make sure we don't have conflicting choices
     if(trim(equilibrium_t_option) == 'from_file' .or. &
        &trim(equilibrium_t_option) == 'JFV' .or. &
        &trim(equilibrium_t_option) == 'strat_file' )then
        pv_sat_flag=.false.
        sat_only_flag=.false.
     endif
     if(do_seasonal_cycle) sc_flag = .true.

!     ----- write version info and namelist to log file -----

      call write_version_number (version,tagname)
      if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=hs_forcing_nml)

      if (no_forcing) return

      twopi = 2*PI

!     ----- convert local heating variables from degrees to radians -----

      xwidth  = local_heating_xwidth*PI/180.
      ywidth  = local_heating_ywidth*PI/180.
      xcenter = local_heating_xcenter*PI/180.
      ycenter = local_heating_ycenter*PI/180.

!     ----- Make sure xcenter falls in the range zero to 2*PI -----

      xcenter = xcenter - twopi*floor(xcenter/twopi)

!     ----- convert local_heating_srfamp from deg/day to deg/sec ----

      srfamp = local_heating_srfamp/SECONDS_PER_DAY
      polar_srfamp = polar_heating_srfamp/SECONDS_PER_DAY

!     ----- compute coefficients -----

! If positive, damping time units are (1/s),  value is the inverse of damping time.
! If negative, damping time units are (days), value is the damping time. It is converted to (1/s)
      
      if (ka < 0.) then
        tka = -1./(86400*ka)
      else
        tka = ka
      endif
      if (ks < 0.) then
        tks = -1./(86400*ks)
      else
        tks = ks
      endif
      if (kf < 0.) then
        vkf = -1./(86400*kf)
      else
        vkf = kf
      endif

!     ----- for tracers -----

      if (trsink < 0.) trsink = -86400.*trsink
      trdamp = 0.; if (trsink > 0.) trdamp = 1./trsink

!     ----- register diagnostic fields -----

      id_teq = register_diag_field ( mod_name, 'teq', axes(1:3), Time, &
                      'equilibrium temperature', 'deg_K'   , &
                      missing_value=missing_value, range=(/50.,400./) )

      id_tau = register_diag_field ( mod_name, 'tau', axes(1:3), Time, &
                      'Newtonian cooling damping time', 'days'   , &
                      missing_value=missing_value      )

      id_newtonian_damping = register_diag_field ( mod_name, 'tdt_ndamp', axes(1:3), Time, &
                      'Heating due to newtonian damping (deg/sec)', 'deg/sec' ,    &
                       missing_value=missing_value     )

      id_tdt = register_diag_field ( mod_name, 'tdt', axes(1:3), Time, &
                      'Total heating: newtonian damping + local heating (deg/sec)', 'deg/sec' ,    &
                       missing_value=missing_value     )

      if(trim(local_heating_option) /= '') then
        id_local_heating=register_diag_field ( mod_name, 'local_heating', axes(1:3), Time, &
                        'Local heating (deg/sec)', 'deg/sec' ,    &
                         missing_value=missing_value     )
      endif

      id_udt = register_diag_field ( mod_name, 'udt_rdamp', axes(1:3), Time, &
                      'zonal wind tendency due to rayleigh damping (m/s2)', 'm/s2',       &
                       missing_value=missing_value     )

      id_vdt = register_diag_field ( mod_name, 'vdt_rdamp', axes(1:3), Time, &
                      'meridional wind tendency due to rayleigh damping (m/s2)', 'm/s2',  &
                       missing_value=missing_value     )

      if (do_conserve_energy) then
         id_tdt_diss = register_diag_field ( mod_name, 'tdt_diss_rdamp', axes(1:3), &
                   Time, 'Dissipative heating from Rayleigh damping (deg/sec)', 'deg/sec',&
                   missing_value=missing_value     )

         id_diss_heat = register_diag_field ( mod_name, 'diss_heat_rdamp', axes(1:2), &
                   Time, 'Vertically integrated dissipative heating from Rayleigh damping (W/m2)', 'W/m2')
      endif


     if(trim(local_heating_option) == 'from_file') then
       call interpolator_init(heating_source_interp, trim(local_heating_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
     endif
     if(trim(equilibrium_t_option) == 'from_file' .or. &
          &trim(equilibrium_t_option) == 'strat_file' ) then
       call interpolator_init (temp_interp, trim(equilibrium_t_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
     endif
     if(relax_to_specified_wind) then
       call interpolator_init (u_interp,    trim(u_wind_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
       call interpolator_init (v_interp,    trim(v_wind_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
     endif

!
! mj read Newtonian time scale from file
!
     if(trim(equilibrium_tau_option) == 'from_file' .or. &
          &trim(equilibrium_tau_option) == 'strat_file' ) then
        call interpolator_init (tau_interp, trim(equilibrium_tau_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/),vert_interp=(/INTERP_LINEAR_P/))
     endif
     

     module_is_initialized  = .true.

 end subroutine hs_forcing_init

!#######################################################################

 subroutine hs_forcing_end 

!-----------------------------------------------------------------------
!
!       routine for terminating held-suarez benchmark module
!             (this routine currently does nothing)
!
!-----------------------------------------------------------------------

 if(trim(local_heating_option) == 'from_file') then
   call interpolator_end(heating_source_interp)
 endif

 if(trim(equilibrium_t_option) == 'from_file' .or. &
      &trim(equilibrium_t_option) == 'strat_file' ) then
   call interpolator_end(temp_interp)
 endif

 if(trim(equilibrium_tau_option) == 'from_file' .or. &
      &trim(equilibrium_tau_option) == 'strat_file' ) then
   call interpolator_end(tau_interp)
 endif


 if(relax_to_specified_wind) then
   call interpolator_end(u_interp)
   call interpolator_end(v_interp)
 endif
 
 module_is_initialized = .false.

 end subroutine hs_forcing_end

!#######################################################################

 subroutine newtonian_damping ( Time, lat, ps, p_full, p_half, t, tdt, teq, tau, mask, surface_forcing )
!-----------------------------------------------------------------------
!
!   routine to compute thermal forcing for held & suarez (1994)
!   benchmark calculation.
!
!-----------------------------------------------------------------------

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:)   :: lat, ps, surface_forcing
real, intent(in),  dimension(:,:,:) :: p_full, t, p_half
real, intent(out), dimension(:,:,:) :: tdt, teq, tau
real, intent(in),  dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------

          real, dimension(size(t,1),size(t,2)) :: &
     sin_lat, sin_lat_2, cos_lat_2, t_star, cos_lat_4, &
     tstr, sigma, the, tfactr, rps, p_norm, dQdistr

       real, dimension(size(t,1),size(t,2),size(t,3)) :: tdamp, dQdamp

       real, dimension(size(t,2),size(t,3)) :: tz,tauz !mj

       integer :: i, j, k, l
       real    :: tcoeff, pref, tcoeff_strat       ! tjr 06/27/03

!----------- US Standard Atmospheric Temperature 1976 ------------------ t.r.
   integer, parameter :: sat_levs = 8
   real, parameter, dimension(sat_levs) :: sat_p = (/ 1.0000000000000, &
                                                      0.2233611050922, &
                                                      0.0540329501078, &
                                                      0.0085666783593, &
                                                      0.0010945601338, &
                                                      0.0006606353133, &
                                                      0.0000390468337, &
                                                      0.0000000000001 /), &
                                           sat_t = (/ 288.150, &
                                                      216.650, &
                                                      216.650, &
                                                      228.650, &
                                                      270.650, &
                                                      270.650, &
                                                      214.650, &
                                                      186.946 /), &
                                           sat_g = (/ -6.5e-3, &
                                                       0.0e-3, &
                                                       1.0e-3, &
                                                       2.8e-3, &
                                                       0.0e-3, &
                                                      -2.8e-3, &
                                                      -2.0e-3, &
                                                       0.0e-3 /)

   real :: t_tropopause =  216.650
   real :: pif = 3.14159265358979/180.
   real :: pid = 3.14159265358979
   real, dimension(size(lat,1),size(lat,2)) :: lat_wgt, t_sat, t_pv
   real :: t0n,t0s,t_days,en,es
   integer :: days,seconds
   integer :: hits
!-------------------------mj Jucker et al (2013) troposphere  ----------
   real,dimension(size(t,1),size(t,2))  :: eps, eps_sc
!-------------------------mj Jucker et al (2014) stratosphere ----------
   real                                          :: p_t=100.e2,p_1=1.e2
   real                                          :: p_n,phipi_N,phipi_S&
        &,delta_phipi
   real                                          :: logphs,logpbd
   real,dimension(size(t,1),size(t,2))           :: P1,P2,P3,P4,Aw,As,D
   real,dimension(size(t,1),size(t,2),size(t,3)) :: Pp
   real,dimension(size(t,1),size(t,2),size(t,3)) :: teq_strat, tau_strat
!   logical,dimension(size(t,1),size(t,2),size(t,3)) :: msk
!-----------------------------------------------------------------------
!------------latitudinal constants--------------------------------------

      sin_lat  (:,:) = sin(lat(:,:))
      sin_lat_2(:,:) = sin_lat(:,:)*sin_lat(:,:)
      cos_lat_2(:,:) = 1.0-sin_lat_2(:,:)
      cos_lat_4(:,:) = cos_lat_2(:,:)*cos_lat_2(:,:)
!mj difference between southern and northern hemisphere
      where( lat(:,:) <= 0. )
         eps = epsS
      elsewhere
         eps = epsN
      end where
!mj seasonal cycle in eps
      if(sc_flag)then
         t0n=0
         t0s=days_per_year/2
!         if (.not.present(Time)) call error_mesg('newtonian_damping','sc_flag true but time not present',FATAL)
         call get_time(Time,seconds,days)
         t_days = days+seconds/86400
         es = max(0.0,sin((t_days-t0s)*2*pid/days_per_year));
         en = max(0.0,sin((t_days-t0n)*2*pid/days_per_year));
         eps_sc = eps*cos(2*pid*t_days/days_per_year)
!         eps_sc =  -eps*sin((t_days-t0s)/360*2*pif*180)
!mj
      else
         eps_sc = eps
      endif

      t_star(:,:) = t_zero - delh*sin_lat_2(:,:) - eps_sc*sin_lat(:,:)
      if ( .not. pv_sat_flag) then
         tstr  (:,:) = t_strat 
      else
         tstr  (:,:) = t_tropopause
      endif

!-----------------------------------------------------------------------
      tcoeff = (tks-tka)/(1.0-sigma_b)
      pref = P00
      rps  = 1./ps

!begin tjr 06/27/03
      tcoeff_strat = (tka-tk_strat)/(sigma_strat1 - sigma_strat2)
!end tjr 06/27/03

      do k = 1, size(t,3)

!  ----- compute equilibrium temperature (teq) ----- !mj troposphere
      if(trim(equilibrium_t_option) == 'Held_Suarez' .or. trim(equilibrium_t_option) == 'JFV' .or. trim(equilibrium_t_option) == 'strat_file' ) then
         p_norm(:,:) = p_full(:,:,k)/pref
         the   (:,:) = t_star(:,:) - delv*cos_lat_2(:,:)*log(p_norm(:,:))
         teq(:,:,k) = the(:,:)*(p_norm(:,:))**KAPPA
         teq(:,:,k) = max( teq(:,:,k), tstr(:,:) )
      else if(equilibrium_t_option == 'from_file') then !mj
         call get_temp(Time, p_half, teq)
      else if(trim(equilibrium_t_option) == 'Constant') then
      	 the(:,:) = t_zero
         teq(:,:,k) = the(:,:) - delh*abs(lat(:,:))/maxval(abs(lat(:,:)))
         teq(:,:,k) = max( teq(:,:,k), tstr(:,:) )
      else
         call error_mesg ('hs_forcing_nml', &
         '"'//trim(equilibrium_t_option)//'"  is not a valid value for equilibrium_t_option',FATAL)
      endif

!  ----- compute damping ----- !mj this is the troposphere
      sigma(:,:) = p_full(:,:,k)*rps(:,:)
      if(trim(equilibrium_tau_option) == 'from_file') then !mj
         call get_tau(Time, p_half, tdamp)
         tdamp = 1./tdamp !tdamp is damping rate, not time
      elseif(trim(equilibrium_tau_option) == 'Held_Suarez'.or. &
        &trim(equilibrium_tau_option) == 'JFV' .or. &
        &trim(equilibrium_tau_option) == 'strat_file' )then
         where (sigma(:,:) <= 1.0 .and. sigma(:,:) > sigma_b)
            tfactr(:,:) = tcoeff*(sigma(:,:)-sigma_b)
            tdamp(:,:,k) = tka + cos_lat_4(:,:)*tfactr(:,:)
         elsewhere(sigma_strat1 <= sigma(:,:) .and. sigma(:,:) < sigma_b )
            tdamp(:,:,k) = tka
         elsewhere(sigma_strat2 <= sigma(:,:) .and. sigma(:,:) < sigma_strat1)
            tfactr(:,:) = tcoeff_strat * (sigma(:,:) - sigma_strat2)
            tdamp(:,:,k) = tk_strat + tfactr(:,:)
         elsewhere
            tdamp(:,:,k) = tk_strat
         endwhere
      else
         call error_mesg ('hs_forcing_nml', &
         '"'//trim(equilibrium_tau_option)//'"  is not a valid value for equilibrium_tau_option',FATAL)         
      endif


      enddo

!  ----- add a US_SAT_1976 stratosphere and a polar vortex to teq ----- tjr
      if ( pv_sat_flag ) then

         if (sc_flag) then
            lat_wgt = 0.5*(1-tanh((lat/pif-sc_phi0n)/sc_dphin))*en + 0.5*(1-tanh((lat/pif-sc_phi0s)/sc_dphis))*es
         else
            lat_wgt = 0.5 * ( 1. - tanh((lat/pif-pv_phi0)/pv_dphi) );
         end if

         do k = 1, size(t,3)
            p_norm(:,:) = p_full(:,:,k)/pref

            do l = 1,sat_levs - 1
               where (p_norm(:,:) < p_tropopause .and. p_norm(:,:) > sat_p(l+1) .and. p_norm(:,:) <= sat_p(l)) 
                  t_sat(:,:) = sat_t(l) * (p_norm(:,:)/sat_p(l))**(-rdgas*sat_g(l)/grav);
               end where
            end do

            if(sat_only_flag)then !mj
               where ( p_norm < p_tropopause ) !mj
                  teq(:,:,k) = t_sat !mj no PV
               endwhere !mj
            else !mj
               where ( p_norm < p_tropopause )
                  t_pv = t_tropopause * (p_norm/p_tropopause)**(rdgas*pv_gamma/grav)
                  teq(:,:,k) = (1-lat_wgt) * t_sat + lat_wgt * t_pv;
               endwhere
            end if
         end do
      end if

!  ----- add Jucker et al (2014) stratosphere to teq ----- mj
      logphs  = log(p_hs)
      logpbd  = log(p_bd)
      if ( trim(equilibrium_t_option) == 'JFV' ) then
         teq_strat = 0.
         phipi_N = phi_N*pif
         phipi_S = phi_S*pif
         where ( lat .ge. phipi_S .and. lat .le. phipi_N )
            P3 = -1.96e-9*(lat/pif)**4 - 1.15e-5*(lat/pif)**2 + 1
         elsewhere ( lat .lt. phipi_S ) 
            P3 = -1.96e-9*(phi_S)**4 - 1.15e-5*(phi_S)**2 + 1
         elsewhere ( lat .gt. phipi_N )
            P3 = -1.96e-9*(phi_N)**4 - 1.15e-5*(phi_N)**2 + 1
         endwhere
         do k=1, size(t,3)
            ! tropical vertical profile
            p_norm(:,:) = log(p_full(:,:,k)/1000.e2)
            P1 = -0.537*p_norm**4 - 9.65*p_norm**3 - 60.6*p_norm**2 - 174*p_norm + 19.8
            P2 =                   0.668*p_norm**3 + 22.3*p_norm**2 + 248*p_norm + 1160.
            where ( p_full(:,:,k) .ge. p_1 )
               teq_strat(:,:,k) = max(t_strat,P1)
            elsewhere
               teq_strat(:,:,k) = max(P1,P2)
            endwhere
            ! latitudinal profile
            p_norm(:,:) = log(p_full(:,:,k)/p_t)/log(p_1/p_t)
            P4 = 1.
            where ( p_1 < p_full(:,:,k) < p_t )
               P4 = (P3 - 1.)*p_norm + 1.
            elsewhere ( p_full(:,:,k) <= p_1 )
               P4 =  P3
            endwhere 
            ! Add polar amplitudes
            p_norm(:,:) = log(p_full(:,:,k)/p_t)/(log(p_1/1000.e2) - log(p_t/1000.e2))
            p_n         = log(p_1          /p_t)/(log(p_1/1000.e2) - log(p_t/1000.e2))
            where ( p_full(:,:,k) .lt. p_1 )
               As = A_s*p_n + A_s
               where ( lat .ge. 0. )
                  Aw = (A_NH_1 - A_NH_0)*p_n + A_NH_0
               elsewhere
                  Aw = (A_SH_1 - A_SH_0)*p_n + A_SH_0
               endwhere
            elsewhere
               As = A_s*p_norm + A_s
               where ( lat .ge. 0. )
                  Aw = (A_NH_1 - A_NH_0)*p_norm + A_NH_0
               elsewhere
                  Aw = (A_SH_1 - A_SH_0)*p_norm + A_SH_0
               endwhere
            endwhere
            Aw = -abs(lat)*Aw*2/PI
            As = -abs(lat)*As*2/PI
            ! add temporal dependence (seasonal cycle)
            if ( do_seasonal_cycle ) then
               call get_time(Time,seconds,days)
               t_days = days+seconds/86400.
            else
               t_days = 0
            endif
            where ( lat .ge. 0. )
               D = cos(twopi*t_days/days_per_year)
            elsewhere
               D = cos(twopi*(t_days-0.5*days_per_year)/days_per_year)
            endwhere
            ! put everything together
            where ( D .ge. 0. )
               teq_strat(:,:,k) = teq_strat(:,:,k)*P4 + &
                    &Aw*D
            elsewhere
               teq_strat(:,:,k) = teq_strat(:,:,k)*P4 + &
                    &(As - teq_strat(:,:,k)*(1. - P4))*D
            endwhere
            ! scale Te profile by T_fact
            teq_strat(:,:,k) = (1.-T_fact)*t_strat + T_fact*teq_strat(:,:,k)
         enddo
      elseif ( trim(equilibrium_t_option) == 'strat_file' ) then
         call get_temp(Time, p_half, teq_strat)
      endif
      if ( trim(equilibrium_t_option) == 'JFV' .or. &
           &trim(equilibrium_t_option) == 'strat_file' )then
         ! merge troposphere and stratosphere
         Pp = 0.
         where ( p_full .le. p_hs .and. p_full .ge. p_bd )
            Pp  = ( logphs - log(p_full) )/( logphs - logpbd )
            teq = Pp*teq_strat + (1. - Pp)*teq
         elsewhere( p_full .lt. p_bd )
            teq = teq_strat   
         end where
      end if
      !! same for damping rate
      if ( trim(equilibrium_tau_option) == 'JFV' ) then
         tau_strat = 0.
         do k=1,size(t,3)
            ! meridonal dependence
            delta_phipi = delta_phi*pif
            where ( lat .lt. 0. .and. lat .ge. phipi_S )
               tau_strat(:,:,k) = tau_S_p + (tau_t - tau_S_p)*exp(-(lat/delta_phipi)**2)
            elsewhere ( lat .lt. phipi_S )
               tau_strat(:,:,k) = tau_S_p + (tau_t - tau_S_p)*exp(-(phi_S/delta_phi)**2)
            elsewhere ( lat .ge. 0. .and. lat .le. phi_N )
               tau_strat(:,:,k) = tau_N_p + (tau_t - tau_N_p)*exp(-(lat/delta_phipi)**2)
            elsewhere ( lat .gt. phi_N )
               tau_strat(:,:,k) = tau_S_p + (tau_t - tau_S_p)*exp(-(phi_N/delta_phi)**2)
            endwhere
            ! vertical dependence
            p_norm = log(p_full(:,:,k)/1000.e2)
            Pp(:,:,k) = 0.045*p_norm**4 + 1.38*p_norm**3 + 15.9*p_norm**2 + 81.6*p_norm + 162
            p_n = log(0.1e2/1000.e2)
            where ( p_full(:,:,k) .lt. 0.1e2 )
               Pp(:,:,k) = 0.045*p_n**4 + 1.38*p_n**3 + 15.9*p_n**2 + 81.6*p_n + 162
            endwhere
            p_n = log(p_t/1000.e2)
            P1 = 0.045*p_n**4 + 1.38*p_n**3 + 15.9*p_n**2 + 81.6*p_n + 162
         enddo
         P2 = minval(Pp,3)
         P3 = P1 - P2
         do k=1,size(t,3)
            P4 = min(1.,( Pp(:,:,k) - P2 )/P3)
            tau_strat(:,:,k) = ( tau_strat(:,:,k) - tau_m )*P4 + tau_m
         enddo
         ! scale profile 
         tau_strat(:,:,k) = (1.-tau_fact)*tau_t + tau_fact*tau_strat(:,:,k)
         ! convert days to seconds
         tau_strat = tau_strat*86400
         ! convert to damping rate
         tau_strat = 1./tau_strat
      elseif ( trim(equilibrium_tau_option) == 'strat_file' ) then
         call get_tau(Time, p_half, tau_strat)
         tau_strat = 1./tau_strat
      endif
      if ( trim(equilibrium_tau_option) == 'strat_file' .or. &
           &trim(equilibrium_tau_option) == 'JFV' ) then
         ! merge with HS troposphere
         Pp = 0.
         where ( p_full .le. p_hs .and. p_full .ge. p_bd )
            Pp  = ( logphs - log(p_full) )/( logphs - logpbd )
            tdamp = Pp*tau_strat + (1. - Pp)*tdamp
         elsewhere( p_full .lt. p_bd )
            tdamp = tau_strat  
         end where
      endif

! mj diagnostics
      tau = 1./tdamp/86400

!  ----- perform Newtonian Cooling -------------

      do k=1,size(t,3)
         tdt(:,:,k) = -tdamp(:,:,k)*(t(:,:,k)-teq(:,:,k))
      enddo

      if (present(mask)) then
         tdt = tdt * mask
         teq = teq * mask
      endif

!-----------------------------------------------------------------------

 end subroutine newtonian_damping

!#######################################################################

 subroutine rayleigh_damping ( Time, ps, p_full, p_half, u, v, udt, vdt, mask )

!-----------------------------------------------------------------------
!
!           rayleigh damping of wind components near surface
!
!-----------------------------------------------------------------------

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:)   :: ps
real, intent(in),  dimension(:,:,:) :: p_full, p_half, u, v
real, intent(out), dimension(:,:,:) :: udt, vdt
real, intent(in),  dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------

real, dimension(size(u,1),size(u,2)) :: sigma, vfactr, sfactr, rps

integer :: i,j,k
real    :: vcoeff
real, dimension(size(u,2),size(u,3)) :: uz, vz, um, vm
real :: umean, vmean

real    :: sponge_coeff !used if sponge_flag
real    :: scdamp
!-----------------------------------------------------------------------
!----------------compute damping----------------------------------------

      if(relax_to_specified_wind) then
        call get_zonal_mean_flow(Time, p_half, uz, vz)
      endif

      vcoeff = -vkf/(1.0-sigma_b)
      rps = 1./ps

      do k = 1, size(u,3)
      if (relax_to_specified_wind) then
         do j=1, size(u,2)
            umean=sum(u(:,j,k))/size(u,1)
            vmean=sum(v(:,j,k))/size(v,1)
            udt(:,j,k) = (uz(j,k)-umean)*vkf
            vdt(:,j,k) = (vz(j,k)-vmean)*vkf
         enddo
      else

         sigma(:,:) = p_full(:,:,k)*rps(:,:)

         where (sigma(:,:) <= 1.0 .and. sigma(:,:) > sigma_b)
            vfactr(:,:) = vcoeff*(sigma(:,:)-sigma_b)
            udt(:,:,k)  = vfactr(:,:)*u(:,:,k)
            vdt(:,:,k)  = vfactr(:,:)*v(:,:,k)
         elsewhere
            udt(:,:,k) = 0.0
            vdt(:,:,k) = 0.0
         endwhere

         if (sponge_flag) then   ! t.r.
            sponge_coeff = 1./sponge_tau_days/86400.
            where (p_full(:,:,k) < sponge_pbottom)
               vfactr(:,:) = -sponge_coeff*(sponge_pbottom-p_full(:,:,k))**2/(sponge_pbottom)**2
               udt(:,:,k) = udt(:,:,k) + vfactr(:,:)*u(:,:,k)
               vdt(:,:,k) = vdt(:,:,k) + vfactr(:,:)*v(:,:,k)
            endwhere
         endif

         if (scaife_flag) then
                scdamp = 1./scaife_damp/86400.
                where( (log(0.1) > log(sigma(:,:))) .and. (log(sigma(:,:)) > log(0.01)))
                   sfactr(:,:) = -scdamp * (log(0.1) - log(sigma(:,:)))
                   udt(:,:,k) = udt(:,:,k) + sfactr(:,:)*u(:,:,k)
                   vdt(:,:,k) = vdt(:,:,k) + sfactr(:,:)*v(:,:,k)
                elsewhere(log(sigma(:,:)) <= log(0.01))
                   sfactr(:,:) = -scdamp
                   udt(:,:,k) = udt(:,:,k) + sfactr(:,:)*u(:,:,k)
                   vdt(:,:,k) = vdt(:,:,k) + sfactr(:,:)*v(:,:,k)
                elsewhere
                   sfactr(:,:) = 0.
                   udt(:,:,k) = udt(:,:,k) + sfactr(:,:)*u(:,:,k)
                   vdt(:,:,k) = vdt(:,:,k) + sfactr(:,:)*v(:,:,k)
                endwhere
         endif
      endif
      enddo

      if (present(mask)) then
          udt = udt * mask
          vdt = vdt * mask
      endif

!-----------------------------------------------------------------------

 end subroutine rayleigh_damping

!#######################################################################

 subroutine sphum_source_sink ( flux, damp, p_full, r, rdt, s_geo, t, dt, kbot )
   use constants_mod !mj
!-----------------------------------------------------------------------
      real, intent(in)  :: flux, damp, p_full(:,:,:), r(:,:,:)
      real, intent(out) :: rdt(:,:,:)
   integer, intent(in), optional :: kbot(:,:)
!-----------------------------------------------------------------------
      real, dimension(size(r,1),size(r,2),size(r,3)) :: source, sink

      integer :: i, j, kb
      real    :: rdamp
!-----------------------------------------------------------------------
      real, intent(in) :: t(:,:,:), dt !mj
      real, intent(in) :: s_geo(:,:) !mj
      real, dimension(size(r,1),size(r,2),size(r,3)) :: qsat !mj
      integer, dimension(size(s_geo,1),size(s_geo,2)) :: sea_surf !mj
!-----------------------------------------------------------------------
      integer,parameter :: si=23, sl=2

      rdamp = damp
      if (rdamp < 0.) rdamp = -86400.*rdamp   ! convert days to seconds
      if (rdamp > 0.) rdamp = 1./rdamp

!------------ surface source and global sink ---------------------------
      
      source(:,:,:)=0.0
      sink(:,:,:)=0.0 !mj
      sea_surf=0 !mj

      where(s_geo < 10.*grav)sea_surf=1 !mj

      if (present(kbot)) then
         do j=1,size(r,2)
            do i=1,size(r,1)
               kb = kbot(i,j)
               qsat(i,j,:) = RDGAS*ES0*exp(-HLV*(1./t(i,j,:)-1./TFREEZE)/RVGAS)/RVGAS !mj
               qsat(i,j,2:kb) = qsat(i,j,2:kb)/p_full(i,j,2:kb) !mj
               qsat(i,j,1) = qsat(i,j,2) !mj
               if(kb.eq.size(r,3))then
                  source(i,j,kb) = max(0.,-vkf*(r(i,j,kb)-qsat(i,j,kb))) !mj specific humidity, above sea only
               endif
            enddo
         enddo
      else
        kb = size(r,3)
        qsat = RDGAS*ES0*exp(-HLV*(1./t - 1./TFREEZE)/RVGAS)/RVGAS !mj
        qsat = qsat/p_full !mj
        source(:,:,kb) = max(0.,-vkf*(r(:,:,kb)-qsat(:,:,kb))) !mj 
        source(:,:,kb) = source(:,:,kb)*sea_surf(:,:) !mj mountains
     endif

     where(r-qsat > 0.)
        sink = (r-qsat)/dt !mj direct complete removal to sat value after one time step
!!         sink = 10.*vkf*(r-qsat) !mj 10 times shorter time scale than source
!         sink = 1. + qsat*HLV*HLV/(RVGAS*CP_AIR*t*t) !mj as in Frierson 2006 JAS
!         sink = (r - qsat)/sink/dt !mj continued
     end where
      
!      sink = sink + rdamp*r !mj add global sink
      
!      sink = 0.
!      sink = r/dt
!      source = 0.

      rdt = source - sink

!-----------------------------------------------------------------------

 end subroutine sphum_source_sink

!-----------------------------------------------------------------------

 subroutine age_source_sink ( flux, damp, p_half, r, rdt, kbot )
!-----------------------------------------------------------------------
      real, intent(in)  :: flux, damp, p_half(:,:,:), r(:,:,:)
      real, intent(out) :: rdt(:,:,:)
   integer, intent(in), optional :: kbot(:,:)
!-----------------------------------------------------------------------
      real, dimension(size(r,1),size(r,2),size(r,3)) :: source,sink
      real, dimension(size(r,1),size(r,2))           :: pmass

      integer :: i, j, kb
      real    :: rdamp
!-----------------------------------------------------------------------
 
      rdamp = damp
      if (rdamp < 0.) rdamp = -86400.*rdamp   ! convert days to seconds
      if (rdamp > 0.) rdamp = 1./rdamp

!------------ simple surface source and no sink --------------------
    
      source=0.0
      sink  =0.0

      if (present(kbot)) then
         do j=1,size(r,2)
            do i=1,size(r,1)
               kb = kbot(i,j)
               pmass (i,j)    = p_half(i,j,kb+1) - p_half(i,j,kb)
               source(i,j,kb) = flux/pmass(i,j)
            enddo
         enddo
      else
         kb = size(r,3)
         pmass (:,:)    = p_half(:,:,kb+1) - p_half(:,:,kb)
         source(:,:,kb) = flux/pmass
      endif

      sink = rdamp*r

      rdt = source - sink      

!-----------------------------------------------------------------------

 end subroutine age_source_sink

!-----------------------------------------------------------------------

 subroutine methane_source_sink ( flux, damp, p_half, r, rdt, kbot )
!-----------------------------------------------------------------------
      real, intent(in)  :: flux, damp, p_half(:,:,:), r(:,:,:)
      real, intent(out) :: rdt(:,:,:)
   integer, intent(in), optional :: kbot(:,:)
!-----------------------------------------------------------------------
      real, dimension(size(r,1),size(r,2),size(r,3)) :: source, sink
      real, dimension(size(r,1),size(r,2))           :: pmass

      integer :: i, j, kb
      real    :: rdamp
!-----------------------------------------------------------------------

      rdamp = damp
      if (rdamp < 0.) rdamp = -86400.*rdamp   ! convert days to seconds
      if (rdamp > 0.) rdamp = 1./rdamp

!------------ surface source and global sink ---------------------------
     
      source(:,:,:)=0.0
      sink(:,:,:)=0.0 !mj


      if (present(kbot)) then
         do j=1,size(r,2)
            do i=1,size(r,1)
               kb = kbot(i,j)
               source(i,j,kb) = -vkf*(r(i,j,kb)-flux) !mj tracer value fixed to trflux at bottom
            enddo
         enddo
      else
        kb = size(r,3)
        source(:,:,kb) = -vkf*(r(:,:,kb)-flux) !mj tracer value fixed to trflux at bottom
     endif
      

     sink = rdamp*r
      

     rdt = source - sink

      

!-----------------------------------------------------------------------

 end subroutine methane_source_sink

!-----------------------------------------------------------------------

 subroutine pv_tracer ( vorn, temp, lat, p_full, dt, r, rdt )
 use constants_mod
!-----------------------------------------------------------------------
      real, intent(in)  :: vorn(:,:,:), temp(:,:,:)!,dZdt(:,:,:), Dtdt(:,:,:)
      real, intent(in)  :: lat(:,:), p_full(:,:,:), dt, r(:,:,:)
      real, intent(out) :: rdt(:,:,:)
!-----------------------------------------------------------------------
      integer :: i,j,k
!      real,dimension(size(rdt,1),size(rdt,2)) :: wa,dTheta
      real :: wa,dTheta
!-----------------------------------------------------------------------
     
      
      do k=2,size(rdt,3) !actual value
         do j=1,size(rdt,2)
            do i=1,size(rdt,1)
               wa = 2*OMEGA*sin(lat(i,j)) + vorn(i,j,k) !absolute vorticity
               dTheta = (1.e5/p_full(i,j,k))**KAPPA*temp(i,j,k) - (1.e5/p_full(i,j,k-1))**KAPPA*temp(i,j,k-1) ! potential temperature variation in p
               rdt(i,j,k) = -GRAV*wa*dTheta/(p_full(i,j,k)-p_full(i,j,k-1)) ! potential vorticity
            enddo
         enddo
      enddo
      rdt(:,:,1)=rdt(:,:,2)
!!$      rdt = vorn
      
      if(dt.gt.0.0)then ! normal pv tracer
         rdt = (rdt-r)/dt
      endif

!-----------------------------------------------------------------------

 end subroutine pv_tracer

!-----------------------------------------------------------------------

 subroutine apv_tracer ( pv, p_half, dt, r, rdt )
!-----------------------------------------------------------------------
      real, intent(in)  :: pv(:,:,:), p_half(:,:,:), dt, r(:,:,:)
      real, intent(out) :: rdt(:,:,:)
!-----------------------------------------------------------------------
      real, dimension(size(rdt,1),size(rdt,2)) :: del_p, apv
      integer :: i,j,k
      real    :: dp,p_low, p_high
!-----------------------------------------------------------------------
      p_low  = 15000.
      p_high = 50000.
      del_p  = 0.
      apv    = 0.
      
      do k=2,size(rdt,3) 
         do j=1,size(rdt,2)
            do i=1,size(rdt,1)
               if(p_half(i,j,k) >= p_low .and. p_half(i,j,k+1) <= p_high)then
                  dp = p_half(i,j,k+1) - p_half(i,j,k)
                  del_p(i,j) = del_p(i,j) + dp
                  apv(i,j) = apv(i,j) + pv(i,j,k)*dp
               endif
            enddo
         enddo
      enddo
      apv = apv/del_p
      do k=2,size(rdt,3)
         rdt(:,:,k) = apv(:,:)
      enddo
      rdt(:,:,1)=rdt(:,:,2)
!!$      rdt = vorn
      
      rdt = (rdt-r)/dt

!-----------------------------------------------------------------------

 end subroutine apv_tracer

!-----------------------------------------------------------------------

 subroutine tracer_source_sink ( flux, damp, p_half, r, rdt, kbot )
!   use mpp_mod, only: mpp_pe, mpp_npes !mj localized source

!-----------------------------------------------------------------------
      real, intent(in)  :: flux, damp, p_half(:,:,:), r(:,:,:)
      real, intent(out) :: rdt(:,:,:)
   integer, intent(in), optional :: kbot(:,:)
!-----------------------------------------------------------------------
      real, dimension(size(r,1),size(r,2),size(r,3)) :: source, sink
      real, dimension(size(r,1),size(r,2))           :: pmass

      integer :: i, j, kb
      real    :: rdamp
!-----------------------------------------------------------------------
!      integer,parameter :: si=23, sl=2 !mj localized source

      rdamp = damp
      if (rdamp < 0.) rdamp = -86400.*rdamp   ! convert days to seconds
      if (rdamp > 0.) rdamp = 1./rdamp

!------------ simple surface source and global sink --------------------

      source(:,:,:)=0.0

   if (present(kbot)) then
      do j=1,size(r,2)
      do i=1,size(r,1)
         kb = kbot(i,j)
         pmass (i,j)    = p_half(i,j,kb+1) - p_half(i,j,kb)
         source(i,j,kb) = flux/pmass(i,j)
      enddo
      enddo
   else
         kb = size(r,3)
         pmass (:,:)    = p_half(:,:,kb+1) - p_half(:,:,kb)
         source(:,:,kb) = flux/pmass(:,:)
   endif

!mj localized source
!     source = 0.
!     pmass (:,:)    = p_half(:,:,si+1) - p_half(:,:,si) ! mj source over tropics
!     if (mpp_pe() == mpp_npes()/2-1)source(:,sl,si) = flux/pmass(:,sl) ! mj source over tropics


     sink(:,:,:) = rdamp*r(:,:,:)
     rdt(:,:,:) = source(:,:,:)-sink(:,:,:)

!-----------------------------------------------------------------------

 end subroutine tracer_source_sink

!#######################################################################

subroutine local_heating ( Time, is, js, lon, lat, ps, p_full, p_half, surf_geopotential, tg, tdt ) !cc
use  press_and_geopot_mod, only: compute_pressures_and_heights

type(time_type), intent(in)         :: Time
integer, intent(in)                 :: is,js
real, intent(in),  dimension(:,:)   :: lon, lat, ps, surf_geopotential !cc
real, intent(in),  dimension(:,:,:) :: p_full, tg !cc
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(out), dimension(:,:,:) :: tdt

integer :: i, j, k
real :: lon_temp, x_temp, p_factor, z_factor, sig_temp !cc
real, dimension(size(lon,1),size(lon,2)) :: lon_factor
real, dimension(size(lat,1),size(lat,2)) :: lat_factor
real, dimension(size(p_full,1),size(p_full,2),size(p_full,3)) :: p_full_dummy, z_full !cc
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)) :: p_half_dummy, z_half !cc
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)) :: p_half2
do i=1,size(p_half,3)
  p_half2(:,:,i)=p_half(:,:,size(p_half,3)-i+1)
enddo

tdt(:,:,:)=0.

if(trim(local_heating_option) == 'from_file') then
   call interpolator( heating_source_interp, Time, p_half, tdt, trim(local_heating_file))
   tdt = tdt/SECONDS_PER_DAY !mj input file is in deg_K/day
else if(trim(local_heating_option) == 'Isidoro') then
   do j=1,size(lon,2)
   do i=1,size(lon,1)
     lon_temp = lon(i,j)
     ! Make sure lon_temp falls in the range zero to 2*PI
     x_temp = floor(lon_temp/twopi)
     lon_temp = lon_temp - twopi*x_temp
     lon_factor(i,j) = exp(-.5*((lon_temp-xcenter)/xwidth)**2)
     lat_factor(i,j) = exp(-.5*((lat(i,j)-ycenter)/ywidth)**2)
     do k=1,size(p_full,3)
       p_factor = exp((p_full(i,j,k)-ps(i,j))/local_heating_vert_decay)
       tdt(i,j,k) = srfamp*lon_factor(i,j)*lat_factor(i,j)*p_factor
     enddo
   enddo
   enddo
else if(trim(local_heating_option) == 'Gaussian') then
   do j=1,size(lon,2)
      do i=1,size(lon,1)
         lat_factor(i,j) = exp( -lat(i,j)**2/(2*(local_heating_latwidth)**2) )
         do k=1,size(p_full,3)
            sig_temp = p_full(i,j,k)/ps(i,j)
            p_factor = exp(-(sig_temp-local_heating_sigcenter)**2/(2*(local_heating_sigwidth)**2) )
            tdt(i,j,k) = srfamp*lat_factor(i,j)*p_factor
         enddo
      enddo
   enddo
else if(trim(local_heating_option) == 'Eichelberger') then
   call compute_pressures_and_heights(tg, ps, surf_geopotential, z_full, z_half, p_full_dummy, p_half_dummy) !cc; dummy variable because I don't want the pressure variables to be overwritten in case it screws the whole model up
   do j=1,size(lon,2)
      do i=1,size(lon,1)
         lat_factor(i,j) = exp( -lat(i,j)**2/(0.3142**2))
         do k=1,size(p_full,3)
            if(z_full(i,j,k)<=12e3) then
               z_factor = sin((pi*z_full(i,j,k))/12e3)
             tdt(i,j,k) = srfamp*lat_factor(i,j)*z_factor
          else 
             tdt(i,j,k) = 0 
          endif
       enddo
    enddo
   enddo
else
   call error_mesg ('hs_forcing_nml','"'//trim(local_heating_option)//'"  is not a valid value for local_heating_option',FATAL)
endif


if(polar_heating_option) then ! Place a half Gaussian at arctic pole to simulate polar amp
   do j=1,size(lon,2)
      do i=1,size(lon,1)
         lat_factor(i,j) = exp( -(lat(i,j)-polar_heating_latcenter)**2/(2*(polar_heating_latwidth)**2) )
         do k=1,size(p_full,3)
            sig_temp = p_full(i,j,k)/ps(i,j)
            p_factor = exp(-(sig_temp-polar_heating_sigcenter)**2/(2*(polar_heating_sigwidth)**2) )
            tdt(i,j,k) = tdt(i,j,k) + polar_srfamp*lat_factor(i,j)*p_factor
         enddo
      enddo
   enddo
endif

end subroutine local_heating

!#######################################################################


!#######################################################################

subroutine get_zonal_mean_flow ( Time, p_half, uz, vz)

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(inout), dimension(:,:) :: uz,vz

integer :: j, k
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: uf,vf
call interpolator( u_interp, p_half, uf, trim(u_wind_file))
call interpolator( v_interp, p_half, vf, trim(v_wind_file))

do j=1,size(p_half,2)
do k=1,size(p_half,3)-1
  uz(j,k)=sum(uf(:,j,k))/size(uf,1)
  vz(j,k)=sum(vf(:,j,k))/size(vf,1)
enddo
enddo
end subroutine get_zonal_mean_flow
!#######################################################################

subroutine get_zonal_mean_temp ( Time, p_half, tm)

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(inout), dimension(:,:) :: tm

integer :: j, k
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: tf
!call interpolator( temp_interp, p_half, tf, trim(equilibrium_t_file))
call interpolator( temp_interp, Time, p_half, tf, trim(equilibrium_t_file))

do j=1,size(p_half,2)
do k=1,size(p_half,3)-1
  tm(j,k)=sum(tf(:,j,k))/size(tf,1)
enddo
enddo
end subroutine get_zonal_mean_temp
!#######################################################################

subroutine get_temp ( Time, p_half, tm)

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(inout), dimension(:,:,:) :: tm

integer :: j, k
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: tf
!call interpolator( temp_interp, p_half, tf, trim(equilibrium_t_file))
call interpolator( temp_interp, Time, p_half, tm, trim(equilibrium_t_file))

end subroutine get_temp
!#######################################################################

subroutine get_zonal_mean_tau ( Time, p_half, taum)

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(inout), dimension(:,:) :: taum

integer :: j, k
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: tauf
!call interpolator( tau_interp, p_half, tauf, trim(equilibrium_tau_file))
call interpolator( tau_interp, Time, p_half, tauf, trim(equilibrium_tau_file))

do j=1,size(p_half,2)
do k=1,size(p_half,3)-1
  taum(j,k)=sum(tauf(:,j,k))/size(tauf,1)
enddo
enddo
end subroutine get_zonal_mean_tau
!#######################################################################

subroutine get_tau ( Time, p_half, taum)

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(inout), dimension(:,:,:) :: taum

integer :: j, k
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: tauf
!call interpolator( tau_interp, p_half, tauf, trim(equilibrium_tau_file))
call interpolator( tau_interp, Time, p_half, taum, trim(equilibrium_tau_file))

end subroutine get_tau
!#######################################################################

end module hs_forcing_mod

