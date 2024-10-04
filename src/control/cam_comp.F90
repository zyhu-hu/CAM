module cam_comp
!-----------------------------------------------------------------------
!
! Community Atmosphere Model (CAM) component interfaces.
!
! This interface layer is CAM specific, i.e., it deals entirely with CAM
! specific data structures.  It is the layer above this, either atm_comp_mct
! or atm_comp_esmf, which translates between CAM and either MCT or ESMF
! data structures in order to interface with the driver/coupler.
!
!-----------------------------------------------------------------------

use shr_kind_mod,      only: r8 => SHR_KIND_R8, cl=>SHR_KIND_CL, cs=>SHR_KIND_CS
use shr_sys_mod,       only: shr_sys_flush

use ESMF,              only: esmf_clock
use seq_timemgr_mod,   only: seq_timemgr_EClockGetData

use spmd_utils,        only: masterproc, mpicom
use cam_control_mod,   only: cam_ctrl_init, cam_ctrl_set_orbit, initial_run
use runtime_opts,      only: read_namelist
use time_manager,      only: timemgr_init, get_step_size, &
                             get_nstep, is_first_step, is_first_restart_step

use camsrfexch,        only: cam_out_t, cam_in_t
use ppgrid,            only: begchunk, endchunk, pver, pcols, pverp
use physics_types,     only: physics_state, physics_tend
use dyn_comp,          only: dyn_import_t, dyn_export_t

use physics_buffer,    only: physics_buffer_desc, pbuf_read_restart, pbuf_init_restart, pbuf_deallocate, pbuf_get_index, pbuf_get_field, pbuf_get_chunk, pbuf_old_tim_idx
use offline_driver,    only: offline_driver_init, offline_driver_dorun, offline_driver_run

use perf_mod
use cam_logfile,       only: iulog
use cam_abortutils,    only: endrun

implicit none
private
save

public cam_init      ! First phase of CAM initialization
public cam_run1      ! CAM run method phase 1
public cam_run2      ! CAM run method phase 2
public cam_run3      ! CAM run method phase 3
public cam_run4      ! CAM run method phase 4
public cam_final     ! CAM Finalization

type(dyn_import_t) :: dyn_in   ! Dynamics import container
type(dyn_export_t) :: dyn_out  ! Dynamics export container

type(physics_state),       pointer :: phys_state(:) => null()
type(physics_tend ),       pointer :: phys_tend(:) => null()
type(physics_buffer_desc), pointer :: pbuf2d(:,:) => null()
type(physics_buffer_desc), pointer :: pbuf(:) => null() ! added

real(r8) :: dtime_phys         ! Time step for physics tendencies.  Set by call to
                               ! stepon_run1, then passed to the phys_run*
integer :: nstep ! 

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

subroutine cam_init(EClock, &
   caseid, ctitle, start_type, dart_mode,       &
   brnch_retain_casename, aqua_planet,          &
   single_column, scmlat, scmlon,               &
   eccen, obliqr, lambm0, mvelpp,               &
   perpetual_run, perpetual_ymd, model_doi_url, &
   cam_out, cam_in)

   !-----------------------------------------------------------------------
   !
   ! CAM component initialization.
   !
   !-----------------------------------------------------------------------

   use history_defaults, only: bldfld
   use cam_initfiles,    only: cam_initfiles_open
   use dyn_grid,         only: dyn_grid_init
   use phys_grid,        only: phys_grid_init
   use physpkg,          only: phys_register, phys_init
   use chem_surfvals,    only: chem_surfvals_init
   use dyn_comp,         only: dyn_init
   use cam_restart,      only: cam_read_restart
   use stepon,           only: stepon_init
   use ionosphere_interface, only: ionosphere_init

#if (defined BFB_CAM_SCAM_IOP)
   use history_defaults, only: initialize_iop_history
#endif

   use camsrfexch,       only: hub2atm_alloc, atm2hub_alloc
   use cam_history,      only: intht
   use history_scam,     only: scm_intht
   use cam_pio_utils,    only: init_pio_subsystem
   use cam_instance,     only: inst_suffix

   ! Arguments
   type(ESMF_Clock),  intent(in) :: EClock

   character(len=cl), intent(in) :: caseid                ! case ID
   character(len=cl), intent(in) :: ctitle                ! case title
   character(len=cs), intent(in) :: start_type            ! start type: initial, restart, or branch
   logical,           intent(in) :: dart_mode             ! enables DART mode
   logical,           intent(in) :: brnch_retain_casename ! Flag to allow a branch to use the same
                                                          ! caseid as the run being branched from.
   logical,           intent(in) :: aqua_planet           ! Flag to run model in "aqua planet" mode

   logical,           intent(in) :: single_column
   real(r8),          intent(in) :: scmlat
   real(r8),          intent(in) :: scmlon

   real(r8),          intent(in) :: eccen
   real(r8),          intent(in) :: obliqr
   real(r8),          intent(in) :: lambm0
   real(r8),          intent(in) :: mvelpp

   logical,           intent(in) :: perpetual_run    ! true => perpetual mode enabled
   integer,           intent(in) :: perpetual_ymd    ! Perpetual date (YYYYMMDD)
   character(len=cl), intent(in) :: model_doi_url    ! CESM model DOI

   type(cam_out_t),   pointer    :: cam_out(:)       ! Output from CAM to surface
   type(cam_in_t) ,   pointer    :: cam_in(:)        ! Merged input state to CAM

   ! Local variables
   character(len=cs) :: filein      ! Input namelist filename
   character(len=cs) :: calendar    ! Calendar type
   integer           :: dtime       ! model timestep (sec)
   integer           :: start_ymd   ! Start date (YYYYMMDD)
   integer           :: start_tod   ! Start time of day (sec)
   integer           :: curr_ymd    ! Start date (YYYYMMDD)
   integer           :: curr_tod    ! Start time of day (sec)
   integer           :: stop_ymd    ! Stop date (YYYYMMDD)
   integer           :: stop_tod    ! Stop time of day (sec)
   integer           :: ref_ymd     ! Reference date (YYYYMMDD)
   integer           :: ref_tod     ! Reference time of day (sec)
   !-----------------------------------------------------------------------

   call init_pio_subsystem()

   ! Initializations using data passed from coupler.
   call cam_ctrl_init( &
      caseid, ctitle, start_type, dart_mode,       &
      aqua_planet, brnch_retain_casename)

   call cam_ctrl_set_orbit(eccen, obliqr, lambm0, mvelpp)

   ! Extract info from the eclock passed from coupler to initialize
   ! the local time manager
   call seq_timemgr_EClockGetData(EClock, &
      start_ymd=start_ymd, start_tod=start_tod, &
      ref_ymd=ref_ymd, ref_tod=ref_tod,         &
      stop_ymd=stop_ymd, stop_tod=stop_tod,     &
      curr_ymd=curr_ymd, curr_tod=curr_tod,     &
      dtime=dtime, calendar=calendar )

   call timemgr_init( &
      dtime, calendar, start_ymd, start_tod, ref_ymd,  &
      ref_tod, stop_ymd, stop_tod, curr_ymd, curr_tod, &
      perpetual_run, perpetual_ymd, initial_run)

   ! Read CAM namelists.
   filein = "atm_in" // trim(inst_suffix)
   call read_namelist(filein, single_column, scmlat, scmlon)

   ! Open initial or restart file, and topo file if specified.
   call cam_initfiles_open()

   ! Initialize grids and dynamics grid decomposition
   call dyn_grid_init()

   ! Initialize physics grid decomposition
   call phys_grid_init()

   ! Register advected tracers and physics buffer fields
   call phys_register ()

   ! Initialize ghg surface values before default initial distributions
   ! are set in dyn_init
   call chem_surfvals_init()

   ! initialize ionosphere
   call ionosphere_init()

   if (initial_run) then

      call dyn_init(dyn_in, dyn_out)

      ! Allocate and setup surface exchange data
      call atm2hub_alloc(cam_out)
      call hub2atm_alloc(cam_in)

   else

      call cam_read_restart(cam_in, cam_out, dyn_in, dyn_out, pbuf2d, stop_ymd, stop_tod)

#if (defined BFB_CAM_SCAM_IOP)
      call initialize_iop_history()
#endif
   end if

   call phys_init( phys_state, phys_tend, pbuf2d,  cam_out )

   call bldfld ()       ! master field list (if branch, only does hash tables)

   call stepon_init(dyn_in, dyn_out)

   call offline_driver_init()

   if (single_column) call scm_intht()
   call intht(model_doi_url)

   if (masterproc) write(iulog,*) 'cam_init complete.'

end subroutine cam_init

!
!-----------------------------------------------------------------------
!
subroutine cam_run1(cam_in, cam_out)
!-----------------------------------------------------------------------
!
! Purpose:   First phase of atmosphere model run method.
!            Runs first phase of dynamics and first phase of
!            physics (before surface model updates).
!
!-----------------------------------------------------------------------

   use physpkg,          only: phys_run1
   use stepon,           only: stepon_run1
   use ionosphere_interface,only: ionosphere_run1

   type(cam_in_t)  :: cam_in(begchunk:endchunk)
   type(cam_out_t) :: cam_out(begchunk:endchunk)

   !-----------------------------------------------------------------------

   if (offline_driver_dorun) return

   !----------------------------------------------------------
   ! First phase of dynamics (at least couple from dynamics to physics)
   ! Return time-step for physics from dynamics.
   !----------------------------------------------------------
   call t_barrierf ('sync_stepon_run1', mpicom)
   call t_startf ('stepon_run1')
   call stepon_run1( dtime_phys, phys_state, phys_tend, pbuf2d, dyn_in, dyn_out )
   call t_stopf  ('stepon_run1')

   !----------------------------------------------------------
   ! first phase of ionosphere -- write to IC file if needed
   !----------------------------------------------------------
   call ionosphere_run1(pbuf2d)

   !
   !----------------------------------------------------------
   ! PHYS_RUN Call the Physics package
   !----------------------------------------------------------
   !

   call t_barrierf ('sync_phys_run1', mpicom)
   call t_startf ('phys_run1')
   call phys_run1(phys_state, dtime_phys, phys_tend, pbuf2d,  cam_in, cam_out)
   call t_stopf  ('phys_run1')

end subroutine cam_run1

!
!-----------------------------------------------------------------------
!

subroutine cam_run2( cam_out, cam_in )
!-----------------------------------------------------------------------
!
! Purpose:   Second phase of atmosphere model run method.
!            Run the second phase physics, run methods that
!            require the surface model updates.  And run the
!            second phase of dynamics that at least couples
!            between physics to dynamics.
!
!-----------------------------------------------------------------------

   use physpkg,          only: phys_run2
   use stepon,           only: stepon_run2
   use ionosphere_interface, only: ionosphere_run2

   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
   type(cam_in_t),  intent(inout) :: cam_in(begchunk:endchunk)

   if (offline_driver_dorun) then
      call offline_driver_run( phys_state, pbuf2d, cam_out, cam_in )
      return
   endif

   !
   ! Second phase of physics (after surface model update)
   !
   call t_barrierf ('sync_phys_run2', mpicom)
   call t_startf ('phys_run2')
   call phys_run2(phys_state, dtime_phys, phys_tend, pbuf2d,  cam_out, cam_in )
   call t_stopf  ('phys_run2')


   !
   ! Second phase of dynamics (at least couple from physics to dynamics)
   !
   call t_barrierf ('sync_stepon_run2', mpicom)
   call t_startf ('stepon_run2')
   call stepon_run2( phys_state, phys_tend, dyn_in, dyn_out )
   call t_stopf  ('stepon_run2')

   !
   ! Ion transport
   !
   call t_startf('ionosphere_run2')
   call ionosphere_run2( phys_state, dyn_in, pbuf2d )
   call t_stopf ('ionosphere_run2')

   if (is_first_step() .or. is_first_restart_step()) then
      call t_startf ('cam_run2_memusage')
      call t_stopf  ('cam_run2_memusage')
   end if
end subroutine cam_run2

!
!-----------------------------------------------------------------------
!

subroutine cam_run3( cam_out )
!-----------------------------------------------------------------------
!
! Purpose:  Third phase of atmosphere model run method. This consists
!           of the third phase of the dynamics. For some dycores
!           this will be the actual dynamics run, for others the
!           dynamics happens before physics in phase 1.
!
!-----------------------------------------------------------------------
   use stepon,           only: stepon_run3

   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
!-----------------------------------------------------------------------

   if (offline_driver_dorun) return

   !
   ! Third phase of dynamics
   !
   call t_barrierf ('sync_stepon_run3', mpicom)
   call t_startf ('stepon_run3')
   call stepon_run3( dtime_phys, cam_out, phys_state, dyn_in, dyn_out )

   call t_stopf  ('stepon_run3')

   if (is_first_step() .or. is_first_restart_step()) then
      call t_startf ('cam_run3_memusage')
      call t_stopf  ('cam_run3_memusage')
   end if
end subroutine cam_run3

!
!-----------------------------------------------------------------------
!

subroutine cam_run4( cam_out, cam_in, rstwr, nlend, &
                     yr_spec, mon_spec, day_spec, sec_spec )

!-----------------------------------------------------------------------
!
! Purpose:  Final phase of atmosphere model run method. This consists
!           of all the restart output, history writes, and other
!           file output.
!
!-----------------------------------------------------------------------
   use cam_history,      only: wshist, wrapup
   use cam_restart,      only: cam_write_restart
   use qneg_module,      only: qneg_print_summary
   use time_manager,     only: is_last_step
   use phys_grid        , only: get_ncols_p ! added
   use filenames        , only: interpret_filename_spec ! added
   use rad_constituents,   only: rad_cnst_get_info ! added

   !type(cam_out_t), intent(inout)        :: cam_out(begchunk:endchunk)
   !type(cam_in_t) , intent(inout)        :: cam_in(begchunk:endchunk)
   type(cam_out_t), pointer        :: cam_out(:) ! added
   type(cam_in_t), pointer          :: cam_in(:) ! added
   logical            , intent(in)           :: rstwr           ! true => write restart file
   logical            , intent(in)           :: nlend           ! true => this is final timestep
   integer            , intent(in), optional :: yr_spec         ! Simulation year
   integer            , intent(in), optional :: mon_spec        ! Simulation month
   integer            , intent(in), optional :: day_spec        ! Simulation day
   integer            , intent(in), optional :: sec_spec        ! Seconds into current simulation day
   logical            , save                 :: do_restart=.FALSE. ! added
   logical            , save                 :: first_time=.true. ! added
   integer :: c, p, q, ncol, lchnk, i, k, m ! added
   character(len=50) :: locfn ! added
   character(len=50) :: filein  ! added
   character(len=50) :: fname_pbuf_cam    ! added
   character(len=50) :: rfilename_spec_cam = '%c.cam.r.%y-%m-%d-%s.nc' ! added
   integer  :: nmodes ! added

   ! define cam.r. buffer fields
   real(r8), pointer, dimension(:) ::  TEOUT         !(lat, lon) ; ! pbuf vars from cam.r. output - sweman
   real(r8), pointer, dimension(:,:) ::  DTCORE   !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CLDO          !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  PRER_EVAP     !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CC_T          !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CC_qv         !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CC_ql         !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CC_qi         !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CC_nl         !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CC_ni         !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CC_qlst       !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  am_evp_st     !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  evprain_st    !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  evpsnow_st    !(pbuf_00032, lat, lon)
   real(r8), pointer, dimension(:) ::  ACPRECL       !(lat, lon) ;
   real(r8), pointer, dimension(:) ::  ACGCME        !(lat, lon) ;
   real(r8), pointer, dimension(:) ::  ACNUM         !(lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  RELVAR        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  ACCRE_ENHAN   !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:) ::  pblh          !(lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  tke           !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  kvh           !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:) ::  tpert         !(lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  AST           !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  AIST          !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  ALST          !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  QIST          !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  QLST          !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CONCLD        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  CLD           !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  RAD_CLUBB     !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  WP2_nadv      !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  WP3_nadv      !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  WPTHLP_nadv   !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  WPRTP_nadv    !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  RTPTHLP_nadv  !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  RTP2_nadv     !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  THLP2_nadv    !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  UP2_nadv      !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  VP2_nadv      !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  UPWP          !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  VPWP          !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  THLM          !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  RTM           !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  UM            !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  VM            !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:,:) ::  DGNUM         !(pbuf_00128, lat, lon) ;
   real(r8), pointer, dimension(:,:,:) ::  DGNUMWET      !(pbuf_00128, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  num_c1        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  so4_c1        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  pom_c1        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  soa_c1        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  bc_c1         !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  dst_c1        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  ncl_c1        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  num_c2        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  so4_c2        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  soa_c2        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  ncl_c2        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  dst_c2        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  num_c3        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  dst_c3        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  ncl_c3        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  so4_c3        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  num_c4        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  pom_c4        !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  bc_c4         !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  DP_FLXPRC     !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  DP_FLXSNW     !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  DP_CLDLIQ     !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  DP_CLDICE     !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:) ::  cush          !(lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  QRS           !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  QRL           !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  ICIWP         !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  ICLWP         !(pbuf_00032, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  kvm           !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  turbtype      !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  smaw          !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:) ::  tauresx       !(lat, lon) ;
   real(r8), pointer, dimension(:) ::  tauresy       !(lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  qpert         !(pbuf_00033, lat, lon) ;
   real(r8), pointer, dimension(:,:) ::  T_TTEND       !(pbuf_00032, lat, lon) ;

   real(r8), pointer, dimension(:) ::  TEOUT_old         !(lat, lon) ; ! pbuf vars from cam.r. output - sweman
  real(r8), pointer, dimension(:,:) ::  DTCORE_old   !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CLDO_old          !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  PRER_EVAP_old     !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CC_T_old          !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CC_qv_old         !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CC_ql_old         !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CC_qi_old         !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CC_nl_old         !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CC_ni_old         !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CC_qlst_old       !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  am_evp_st_old     !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  evprain_st_old    !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  evpsnow_st_old    !(pbuf_00032, lat, lon)
  real(r8), pointer, dimension(:) ::  ACPRECL_old       !(lat, lon) ;
  real(r8), pointer, dimension(:) ::  ACGCME_old        !(lat, lon) ;
  real(r8), pointer, dimension(:) ::  ACNUM_old         !(lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  RELVAR_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  ACCRE_ENHAN_old   !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:) ::  pblh_old          !(lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  tke_old           !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  kvh_old           !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:) ::  tpert_old         !(lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  AST_old           !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  AIST_old          !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  ALST_old          !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  QIST_old          !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  QLST_old          !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CONCLD_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  CLD_old           !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  RAD_CLUBB_old     !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  WP2_nadv_old      !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  WP3_nadv_old      !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  WPTHLP_nadv_old   !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  WPRTP_nadv_old    !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  RTPTHLP_nadv_old  !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  RTP2_nadv_old     !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  THLP2_nadv_old    !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  UP2_nadv_old      !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  VP2_nadv_old      !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  UPWP_old          !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  VPWP_old          !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  THLM_old          !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  RTM_old           !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  UM_old            !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  VM_old            !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:,:) ::  DGNUM_old         !(pbuf_00128, lat, lon) ;
  real(r8), pointer, dimension(:,:,:) ::  DGNUMWET_old      !(pbuf_00128, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  num_c1_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  so4_c1_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  pom_c1_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  soa_c1_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  bc_c1_old         !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  dst_c1_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  ncl_c1_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  num_c2_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  so4_c2_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  soa_c2_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  ncl_c2_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  dst_c2_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  num_c3_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  dst_c3_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  ncl_c3_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  so4_c3_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  num_c4_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  pom_c4_old        !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  bc_c4_old         !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  DP_FLXPRC_old     !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  DP_FLXSNW_old     !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  DP_CLDLIQ_old     !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  DP_CLDICE_old     !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:) ::  cush_old          !(lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  QRS_old           !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  QRL_old           !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  ICIWP_old         !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  ICLWP_old         !(pbuf_00032, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  kvm_old           !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  turbtype_old      !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  smaw_old          !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:) ::  tauresx_old       !(lat, lon) ;
  real(r8), pointer, dimension(:) ::  tauresy_old       !(lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  qpert_old         !(pbuf_00033, lat, lon) ;
  real(r8), pointer, dimension(:,:) ::  T_TTEND_old       !(pbuf_00032, lat, lon) ;

  integer ::  TEOUT_idx        = 0 !(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
  integer ::  DTCORE_idx  = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CLDO_idx         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  PRER_EVAP_idx    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_T_idx         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_qv_idx        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_ql_idx        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_qi_idx        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_nl_idx        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_ni_idx        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_qlst_idx      = 0 !(pbuf_00032, lat, lon) ;
  integer ::  am_evp_st_idx    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  evprain_st_idx   = 0 !(pbuf_00032, lat, lon) ;
  integer ::  evpsnow_st_idx   = 0 !(pbuf_00032, lat, lon)
  integer ::  ACPRECL_idx      = 0 !(lat, lon) ;
  integer ::  ACGCME_idx       = 0 !(lat, lon) ;
  integer ::  ACNUM_idx        = 0 !(lat, lon) ;
  integer ::  RELVAR_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ACCRE_ENHAN_idx  = 0 !(pbuf_00032, lat, lon) ;
  integer ::  pblh_idx         = 0 !(lat, lon) ;
  integer ::  tke_idx          = 0 !(pbuf_00033, lat, lon) ;
  integer ::  kvh_idx          = 0 !(pbuf_00033, lat, lon) ;
  integer ::  tpert_idx        = 0 !(lat, lon) ;
  integer ::  AST_idx          = 0 !(pbuf_00032, lat, lon) ;
  integer ::  AIST_idx         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ALST_idx         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  QIST_idx         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  QLST_idx         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CONCLD_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CLD_idx          = 0 !(pbuf_00032, lat, lon) ;
  integer ::  RAD_CLUBB_idx    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  WP2_nadv_idx     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  WP3_nadv_idx     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  WPTHLP_nadv_idx  = 0 !(pbuf_00033, lat, lon) ;
  integer ::  WPRTP_nadv_idx   = 0 !(pbuf_00033, lat, lon) ;
  integer ::  RTPTHLP_nadv_idx = 0 !(pbuf_00033, lat, lon) ;
  integer ::  RTP2_nadv_idx    = 0 !(pbuf_00033, lat, lon) ;
  integer ::  THLP2_nadv_idx   = 0 !(pbuf_00033, lat, lon) ;
  integer ::  UP2_nadv_idx     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  VP2_nadv_idx     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  UPWP_idx         = 0 !(pbuf_00033, lat, lon) ;
  integer ::  VPWP_idx         = 0 !(pbuf_00033, lat, lon) ;
  integer ::  THLM_idx         = 0 !(pbuf_00033, lat, lon) ;
  integer ::  RTM_idx          = 0 !(pbuf_00033, lat, lon) ;
  integer ::  UM_idx           = 0 !(pbuf_00033, lat, lon) ;
  integer ::  VM_idx           = 0 !(pbuf_00033, lat, lon) ;
  integer ::  DGNUM_idx        = 0 !(pbuf_00128, lat, lon) ;
  integer ::  DGNUMWET_idx     = 0 !(pbuf_00128, lat, lon) ;
  integer ::  num_c1_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  so4_c1_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  pom_c1_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  soa_c1_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  bc_c1_idx        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  dst_c1_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ncl_c1_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  num_c2_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  so4_c2_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  soa_c2_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ncl_c2_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  dst_c2_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  num_c3_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  dst_c3_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ncl_c3_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  so4_c3_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  num_c4_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  pom_c4_idx       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  bc_c4_idx        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  DP_FLXPRC_idx    = 0 !(pbuf_00033, lat, lon) ;
  integer ::  DP_FLXSNW_idx    = 0 !(pbuf_00033, lat, lon) ;
  integer ::  DP_CLDLIQ_idx    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  DP_CLDICE_idx    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  cush_idx         = 0 !(lat, lon) ;
  integer ::  QRS_idx          = 0 !(pbuf_00032, lat, lon) ;
  integer ::  QRL_idx          = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ICIWP_idx        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ICLWP_idx        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  kvm_idx          = 0 !(pbuf_00033, lat, lon) ;
  integer ::  turbtype_idx     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  smaw_idx         = 0 !(pbuf_00033, lat, lon) ;
  integer ::  tauresx_idx      = 0 !(lat, lon) ;
  integer ::  tauresy_idx      = 0 !(lat, lon) ;
  integer ::  qpert_idx        = 0 !(pbuf_00033, lat, lon) ;
  integer ::  T_TTEND_idx      = 0 !(pbuf_00032, lat, lon) ;

  integer ::  TEOUT_oldid        = 0 !(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
  integer ::  DTCORE_oldid  = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CLDO_oldid         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  PRER_EVAP_oldid    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_T_oldid         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_qv_oldid        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_ql_oldid        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_qi_oldid        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_nl_oldid        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_ni_oldid        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CC_qlst_oldid      = 0 !(pbuf_00032, lat, lon) ;
  integer ::  am_evp_st_oldid    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  evprain_st_oldid   = 0 !(pbuf_00032, lat, lon) ;
  integer ::  evpsnow_st_oldid   = 0 !(pbuf_00032, lat, lon)
  integer ::  ACPRECL_oldid      = 0 !(lat, lon) ;
  integer ::  ACGCME_oldid       = 0 !(lat, lon) ;
  integer ::  ACNUM_oldid        = 0 !(lat, lon) ;
  integer ::  RELVAR_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ACCRE_ENHAN_oldid  = 0 !(pbuf_00032, lat, lon) ;
  integer ::  pblh_oldid         = 0 !(lat, lon) ;
  integer ::  tke_oldid          = 0 !(pbuf_00033, lat, lon) ;
  integer ::  kvh_oldid          = 0 !(pbuf_00033, lat, lon) ;
  integer ::  tpert_oldid        = 0 !(lat, lon) ;
  integer ::  AST_oldid          = 0 !(pbuf_00032, lat, lon) ;
  integer ::  AIST_oldid         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ALST_oldid         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  QIST_oldid         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  QLST_oldid         = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CONCLD_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  CLD_oldid          = 0 !(pbuf_00032, lat, lon) ;
  integer ::  RAD_CLUBB_oldid    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  WP2_nadv_oldid     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  WP3_nadv_oldid     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  WPTHLP_nadv_oldid  = 0 !(pbuf_00033, lat, lon) ;
  integer ::  WPRTP_nadv_oldid   = 0 !(pbuf_00033, lat, lon) ;
  integer ::  RTPTHLP_nadv_oldid = 0 !(pbuf_00033, lat, lon) ;
  integer ::  RTP2_nadv_oldid    = 0 !(pbuf_00033, lat, lon) ;
  integer ::  THLP2_nadv_oldid   = 0 !(pbuf_00033, lat, lon) ;
  integer ::  UP2_nadv_oldid     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  VP2_nadv_oldid     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  UPWP_oldid         = 0 !(pbuf_00033, lat, lon) ;
  integer ::  VPWP_oldid         = 0 !(pbuf_00033, lat, lon) ;
  integer ::  THLM_oldid         = 0 !(pbuf_00033, lat, lon) ;
  integer ::  RTM_oldid          = 0 !(pbuf_00033, lat, lon) ;
  integer ::  UM_oldid           = 0 !(pbuf_00033, lat, lon) ;
  integer ::  VM_oldid           = 0 !(pbuf_00033, lat, lon) ;
  integer ::  DGNUM_oldid        = 0 !(pbuf_00128, lat, lon) ;
  integer ::  DGNUMWET_oldid     = 0 !(pbuf_00128, lat, lon) ;
  integer ::  num_c1_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  so4_c1_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  pom_c1_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  soa_c1_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  bc_c1_oldid        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  dst_c1_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ncl_c1_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  num_c2_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  so4_c2_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  soa_c2_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ncl_c2_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  dst_c2_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  num_c3_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  dst_c3_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ncl_c3_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  so4_c3_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  num_c4_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  pom_c4_oldid       = 0 !(pbuf_00032, lat, lon) ;
  integer ::  bc_c4_oldid        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  DP_FLXPRC_oldid    = 0 !(pbuf_00033, lat, lon) ;
  integer ::  DP_FLXSNW_oldid    = 0 !(pbuf_00033, lat, lon) ;
  integer ::  DP_CLDLIQ_oldid    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  DP_CLDICE_oldid    = 0 !(pbuf_00032, lat, lon) ;
  integer ::  cush_oldid         = 0 !(lat, lon) ;
  integer ::  QRS_oldid          = 0 !(pbuf_00032, lat, lon) ;
  integer ::  QRL_oldid          = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ICIWP_oldid        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  ICLWP_oldid        = 0 !(pbuf_00032, lat, lon) ;
  integer ::  kvm_oldid          = 0 !(pbuf_00033, lat, lon) ;
  integer ::  turbtype_oldid     = 0 !(pbuf_00033, lat, lon) ;
  integer ::  smaw_oldid         = 0 !(pbuf_00033, lat, lon) ;
  integer ::  tauresx_oldid      = 0 !(lat, lon) ;
  integer ::  tauresy_oldid      = 0 !(lat, lon) ;
  integer ::  qpert_oldid        = 0 !(pbuf_00033, lat, lon) ;
  integer ::  T_TTEND_oldid      = 0 !(pbuf_00032, lat, lon) ;

   !----------------------------------------------------------
   ! History and restart logic: Write and/or dispose history tapes if required
   !----------------------------------------------------------
   !
   call t_barrierf ('sync_wshist', mpicom)
   call t_startf ('wshist')
   call wshist ()
   call t_stopf  ('wshist')

   !
   ! Write restart files
   !
   if (rstwr) then
      call t_startf ('cam_write_restart')
      if (present(yr_spec).and.present(mon_spec).and.present(day_spec).and.present(sec_spec)) then
         call cam_write_restart(cam_in, cam_out, dyn_out, pbuf2d, &
              yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
      else
         call cam_write_restart(cam_in, cam_out, dyn_out, pbuf2d )
      end if
      call t_stopf  ('cam_write_restart')
   end if

   !when it's time to do restart, reset buffer variables - sweidman

   if ( mod(sec_spec,21600)==0  .AND. .NOT. do_restart ) then 
      do_restart=.TRUE.

      if(masterproc) then
         print *, 'swap pbuf old to new ', sec_spec
      end if
      
      do lchnk=begchunk,endchunk 
         pbuf=> pbuf_get_chunk(pbuf2d,lchnk)
     
         ! do the same as above, get buffer values
         TEOUT_idx        = pbuf_get_index('TEOUT')!(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
         DTCORE_idx       = pbuf_get_index('DTCORE')!(pbuf_00032, lat, lon) ;
         CLDO_idx         = pbuf_get_index('CLDO')!(pbuf_00032, lat, lon) ;
         PRER_EVAP_idx    = pbuf_get_index('PRER_EVAP')!(pbuf_00032, lat, lon) ;
         CC_T_idx         = pbuf_get_index('CC_T')!(pbuf_00032, lat, lon) ;
         CC_qv_idx        = pbuf_get_index('CC_qv')!(pbuf_00032, lat, lon) ;
         CC_ql_idx        = pbuf_get_index('CC_ql')!(pbuf_00032, lat, lon) ;
         CC_qi_idx        = pbuf_get_index('CC_qi')!(pbuf_00032, lat, lon) ;
         CC_nl_idx        = pbuf_get_index('CC_nl')!(pbuf_00032, lat, lon) ;
         CC_ni_idx        = pbuf_get_index('CC_ni')!(pbuf_00032, lat, lon) ;
         CC_qlst_idx      = pbuf_get_index('CC_qlst')!(pbuf_00032, lat, lon) ;
         am_evp_st_idx    = pbuf_get_index('am_evp_st')!(pbuf_00032, lat, lon) ;
         evprain_st_idx   = pbuf_get_index('evprain_st')!(pbuf_00032, lat, lon) ;
         evpsnow_st_idx   = pbuf_get_index('evpsnow_st')!(pbuf_00032, lat, lon)
         ACPRECL_idx      = pbuf_get_index('ACPRECL')!(lat, lon) ;
         ACGCME_idx       = pbuf_get_index('ACGCME')!(lat, lon) ;
         ACNUM_idx        = pbuf_get_index('ACNUM')!(lat, lon) ;
         RELVAR_idx       = pbuf_get_index('RELVAR')!(pbuf_00032, lat, lon) ;
         ACCRE_ENHAN_idx  = pbuf_get_index('ACCRE_ENHAN')!(pbuf_00032, lat, lon) ;
         pblh_idx         = pbuf_get_index('pblh')!(lat, lon) ;
         tke_idx          = pbuf_get_index('tke')!(pbuf_00033, lat, lon) ;
         kvh_idx          = pbuf_get_index('kvh')!(pbuf_00033, lat, lon) ;
         tpert_idx        = pbuf_get_index('tpert')!(lat, lon) ;
         AST_idx          = pbuf_get_index('AST')!(pbuf_00032, lat, lon) ;
         AIST_idx         = pbuf_get_index('AIST')!(pbuf_00032, lat, lon) ;
         ALST_idx         = pbuf_get_index('ALST')!(pbuf_00032, lat, lon) ;
         QIST_idx         = pbuf_get_index('QIST')!(pbuf_00032, lat, lon) ;
         QLST_idx         = pbuf_get_index('QLST')!(pbuf_00032, lat, lon) ;
         CONCLD_idx       = pbuf_get_index('CONCLD')!(pbuf_00032, lat, lon) ;
         CLD_idx          = pbuf_get_index('CLD')!(pbuf_00032, lat, lon) ;
         RAD_CLUBB_idx    = pbuf_get_index('RAD_CLUBB')!(pbuf_00032, lat, lon) ;
         WP2_nadv_idx     = pbuf_get_index('WP2_nadv')!(pbuf_00033, lat, lon) ;
         WP3_nadv_idx     = pbuf_get_index('WP3_nadv')!(pbuf_00033, lat, lon) ;
         WPTHLP_nadv_idx  = pbuf_get_index('WPTHLP_nadv')!(pbuf_00033, lat, lon) ;
         WPRTP_nadv_idx   = pbuf_get_index('WPRTP_nadv')!(pbuf_00033, lat, lon) ;
         RTPTHLP_nadv_idx = pbuf_get_index('RTPTHLP_nadv')!(pbuf_00033, lat, lon) ;
         RTP2_nadv_idx    = pbuf_get_index('RTP2_nadv')!(pbuf_00033, lat, lon) ;
         THLP2_nadv_idx   = pbuf_get_index('THLP2_nadv')!(pbuf_00033, lat, lon) ;
         UP2_nadv_idx     = pbuf_get_index('UP2_nadv')!(pbuf_00033, lat, lon) ;
         VP2_nadv_idx     = pbuf_get_index('VP2_nadv')!(pbuf_00033, lat, lon) ;
         UPWP_idx         = pbuf_get_index('UPWP')!(pbuf_00033, lat, lon) ;
         VPWP_idx         = pbuf_get_index('VPWP')!(pbuf_00033, lat, lon) ;
         THLM_idx         = pbuf_get_index('THLM')!(pbuf_00033, lat, lon) ;
         RTM_idx          = pbuf_get_index('RTM')!(pbuf_00033, lat, lon) ;
         UM_idx           = pbuf_get_index('UM')!(pbuf_00033, lat, lon) ;
         VM_idx           = pbuf_get_index('VM')!(pbuf_00033, lat, lon) ;
         DGNUM_idx        = pbuf_get_index('DGNUM')!(pbuf_00128, lat, lon) ;
         DGNUMWET_idx     = pbuf_get_index('DGNUMWET')!(pbuf_00128, lat, lon) ;
         num_c1_idx       = pbuf_get_index('num_c1')!(pbuf_00032, lat, lon) ;
         so4_c1_idx       = pbuf_get_index('so4_c1')!(pbuf_00032, lat, lon) ;
         pom_c1_idx       = pbuf_get_index('pom_c1')!(pbuf_00032, lat, lon) ;
         soa_c1_idx       = pbuf_get_index('soa_c1')!(pbuf_00032, lat, lon) ;
         bc_c1_idx        = pbuf_get_index('bc_c1')!(pbuf_00032, lat, lon) ;
         dst_c1_idx       = pbuf_get_index('dst_c1')!(pbuf_00032, lat, lon) ;
         ncl_c1_idx       = pbuf_get_index('ncl_c1')!(pbuf_00032, lat, lon) ;
         num_c2_idx       = pbuf_get_index('num_c2')!(pbuf_00032, lat, lon) ;
         so4_c2_idx       = pbuf_get_index('so4_c2')!(pbuf_00032, lat, lon) ;
         soa_c2_idx       = pbuf_get_index('soa_c2')!(pbuf_00032, lat, lon) ;
         ncl_c2_idx       = pbuf_get_index('ncl_c2')!(pbuf_00032, lat, lon) ;
         dst_c2_idx       = pbuf_get_index('dst_c2')!(pbuf_00032, lat, lon) ;
         num_c3_idx       = pbuf_get_index('num_c3')!(pbuf_00032, lat, lon) ;
         dst_c3_idx       = pbuf_get_index('dst_c3')!(pbuf_00032, lat, lon) ;
         ncl_c3_idx       = pbuf_get_index('ncl_c3')!(pbuf_00032, lat, lon) ;
         so4_c3_idx       = pbuf_get_index('so4_c3')!(pbuf_00032, lat, lon) ;
         num_c4_idx       = pbuf_get_index('num_c4')!(pbuf_00032, lat, lon) ;
         pom_c4_idx       = pbuf_get_index('pom_c4')!(pbuf_00032, lat, lon) ;
         bc_c4_idx        = pbuf_get_index('bc_c4')!(pbuf_00032, lat, lon) ;
         DP_FLXPRC_idx    = pbuf_get_index('DP_FLXPRC')!(pbuf_00033, lat, lon) ;
         DP_FLXSNW_idx    = pbuf_get_index('DP_FLXSNW')!(pbuf_00033, lat, lon) ;
         DP_CLDLIQ_idx    = pbuf_get_index('DP_CLDLIQ')!(pbuf_00032, lat, lon) ;
         DP_CLDICE_idx    = pbuf_get_index('DP_CLDICE')!(pbuf_00032, lat, lon) ;
         cush_idx         = pbuf_get_index('cush')!(lat, lon) ;
         QRS_idx          = pbuf_get_index('QRS')!(pbuf_00032, lat, lon) ;
         QRL_idx          = pbuf_get_index('QRL')!(pbuf_00032, lat, lon) ;
         ICIWP_idx        = pbuf_get_index('ICIWP')!(pbuf_00032, lat, lon) ;
         ICLWP_idx        = pbuf_get_index('ICLWP')!(pbuf_00032, lat, lon) ;
         kvm_idx          = pbuf_get_index('kvm')!(pbuf_00033, lat, lon) ;
         turbtype_idx     = pbuf_get_index('turbtype')!(pbuf_00033, lat, lon) ;
         smaw_idx         = pbuf_get_index('smaw')!(pbuf_00033, lat, lon) ;
         tauresx_idx      = pbuf_get_index('tauresx')!(lat, lon) ;
         tauresy_idx      = pbuf_get_index('tauresy')!(lat, lon) ;
         qpert_idx        = pbuf_get_index('qpert')!(pbuf_00033, lat, lon) ;
         T_TTEND_idx      = pbuf_get_index('T_TTEND')!(pbuf_00032, lat, lon) ;

         TEOUT_oldid        = pbuf_get_index('TEOUT_OLD')!(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
         DTCORE_oldid       = pbuf_get_index('DTCORE_OLD')!(pbuf_00032, lat, lon) ;
         CLDO_oldid         = pbuf_get_index('CLDO_OLD')!(pbuf_00032, lat, lon) ;
         PRER_EVAP_oldid    = pbuf_get_index('PRER_EVAP_OLD')!(pbuf_00032, lat, lon) ;
         CC_T_oldid         = pbuf_get_index('CC_T_OLD')!(pbuf_00032, lat, lon) ;
         CC_qv_oldid        = pbuf_get_index('CC_qv_OLD')!(pbuf_00032, lat, lon) ;
         CC_ql_oldid        = pbuf_get_index('CC_ql_OLD')!(pbuf_00032, lat, lon) ;
         CC_qi_oldid        = pbuf_get_index('CC_qi_OLD')!(pbuf_00032, lat, lon) ;
         CC_nl_oldid        = pbuf_get_index('CC_nl_OLD')!(pbuf_00032, lat, lon) ;
         CC_ni_oldid        = pbuf_get_index('CC_ni_OLD')!(pbuf_00032, lat, lon) ;
         CC_qlst_oldid      = pbuf_get_index('CC_qlst_OLD')!(pbuf_00032, lat, lon) ;
         am_evp_st_oldid    = pbuf_get_index('am_evp_st_OLD')!(pbuf_00032, lat, lon) ;
         evprain_st_oldid   = pbuf_get_index('evprain_st_OLD')!(pbuf_00032, lat, lon) ;
         evpsnow_st_oldid   = pbuf_get_index('evpsnow_st_OLD')!(pbuf_00032, lat, lon)
         ACPRECL_oldid      = pbuf_get_index('ACPRECL_OLD')!(lat, lon) ;
         ACGCME_oldid       = pbuf_get_index('ACGCME_OLD')!(lat, lon) ;
         ACNUM_oldid        = pbuf_get_index('ACNUM_OLD')!(lat, lon) ;
         RELVAR_oldid       = pbuf_get_index('RELVAR_OLD')!(pbuf_00032, lat, lon) ;
         ACCRE_ENHAN_oldid  = pbuf_get_index('ACCRE_ENHAN_OLD')!(pbuf_00032, lat, lon) ;
         pblh_oldid         = pbuf_get_index('pblh_OLD')!(lat, lon) ;
         tke_oldid          = pbuf_get_index('tke_OLD')!(pbuf_00033, lat, lon) ;
         kvh_oldid          = pbuf_get_index('kvh_OLD')!(pbuf_00033, lat, lon) ;
         tpert_oldid        = pbuf_get_index('tpert_OLD')!(lat, lon) ;
         AST_oldid          = pbuf_get_index('AST_OLD')!(pbuf_00032, lat, lon) ;
         AIST_oldid         = pbuf_get_index('AIST_OLD')!(pbuf_00032, lat, lon) ;
         ALST_oldid         = pbuf_get_index('ALST_OLD')!(pbuf_00032, lat, lon) ;
         QIST_oldid         = pbuf_get_index('QIST_OLD')!(pbuf_00032, lat, lon) ;
         QLST_oldid         = pbuf_get_index('QLST_OLD')!(pbuf_00032, lat, lon) ;
         CONCLD_oldid       = pbuf_get_index('CONCLD_OLD')!(pbuf_00032, lat, lon) ;
         CLD_oldid          = pbuf_get_index('CLD_OLD')!(pbuf_00032, lat, lon) ;
         RAD_CLUBB_oldid    = pbuf_get_index('RAD_CLUBB_OLD')!(pbuf_00032, lat, lon) ;
         WP2_nadv_oldid     = pbuf_get_index('WP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
         WP3_nadv_oldid     = pbuf_get_index('WP3_nadv_OLD')!(pbuf_00033, lat, lon) ;
         WPTHLP_nadv_oldid  = pbuf_get_index('WPTHLP_nadv_OLD')!(pbuf_00033, lat, lon) ;
         WPRTP_nadv_oldid   = pbuf_get_index('WPRTP_nadv_OLD')!(pbuf_00033, lat, lon) ;
         RTPTHLP_nadv_oldid = pbuf_get_index('RTPTHLP_nadv_OLD')!(pbuf_00033, lat, lon) ;
         RTP2_nadv_oldid    = pbuf_get_index('RTP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
         THLP2_nadv_oldid   = pbuf_get_index('THLP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
         UP2_nadv_oldid     = pbuf_get_index('UP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
         VP2_nadv_oldid     = pbuf_get_index('VP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
         UPWP_oldid         = pbuf_get_index('UPWP_OLD')!(pbuf_00033, lat, lon) ;
         VPWP_oldid         = pbuf_get_index('VPWP_OLD')!(pbuf_00033, lat, lon) ;
         THLM_oldid         = pbuf_get_index('THLM_OLD')!(pbuf_00033, lat, lon) ;
         RTM_oldid          = pbuf_get_index('RTM_OLD')!(pbuf_00033, lat, lon) ;
         UM_oldid           = pbuf_get_index('UM_OLD')!(pbuf_00033, lat, lon) ;
         VM_oldid           = pbuf_get_index('VM_OLD')!(pbuf_00033, lat, lon) ;
         DGNUM_oldid        = pbuf_get_index('DGNUM_OLD')!(pbuf_00128, lat, lon) ;
         DGNUMWET_oldid     = pbuf_get_index('DGNUMWET_OLD')!(pbuf_00128, lat, lon) ;
         num_c1_oldid       = pbuf_get_index('num_c1_OLD')!(pbuf_00032, lat, lon) ;
         so4_c1_oldid       = pbuf_get_index('so4_c1_OLD')!(pbuf_00032, lat, lon) ;
         pom_c1_oldid       = pbuf_get_index('pom_c1_OLD')!(pbuf_00032, lat, lon) ;
         soa_c1_oldid       = pbuf_get_index('soa_c1_OLD')!(pbuf_00032, lat, lon) ;
         bc_c1_oldid        = pbuf_get_index('bc_c1_OLD')!(pbuf_00032, lat, lon) ;
         dst_c1_oldid       = pbuf_get_index('dst_c1_OLD')!(pbuf_00032, lat, lon) ;
         ncl_c1_oldid       = pbuf_get_index('ncl_c1_OLD')!(pbuf_00032, lat, lon) ;
         num_c2_oldid       = pbuf_get_index('num_c2_OLD')!(pbuf_00032, lat, lon) ;
         so4_c2_oldid       = pbuf_get_index('so4_c2_OLD')!(pbuf_00032, lat, lon) ;
         soa_c2_oldid       = pbuf_get_index('soa_c2_OLD')!(pbuf_00032, lat, lon) ;
         ncl_c2_oldid       = pbuf_get_index('ncl_c2_OLD')!(pbuf_00032, lat, lon) ;
         dst_c2_oldid       = pbuf_get_index('dst_c2_OLD')!(pbuf_00032, lat, lon) ;
         num_c3_oldid       = pbuf_get_index('num_c3_OLD')!(pbuf_00032, lat, lon) ;
         dst_c3_oldid       = pbuf_get_index('dst_c3_OLD')!(pbuf_00032, lat, lon) ;
         ncl_c3_oldid       = pbuf_get_index('ncl_c3_OLD')!(pbuf_00032, lat, lon) ;
         so4_c3_oldid       = pbuf_get_index('so4_c3_OLD')!(pbuf_00032, lat, lon) ;
         num_c4_oldid       = pbuf_get_index('num_c4_OLD')!(pbuf_00032, lat, lon) ;
         pom_c4_oldid       = pbuf_get_index('pom_c4_OLD')!(pbuf_00032, lat, lon) ;
         bc_c4_oldid        = pbuf_get_index('bc_c4_OLD')!(pbuf_00032, lat, lon) ;
         DP_FLXPRC_oldid    = pbuf_get_index('DP_FLXPRC_OLD')!(pbuf_00033, lat, lon) ;
         DP_FLXSNW_oldid    = pbuf_get_index('DP_FLXSNW_OLD')!(pbuf_00033, lat, lon) ;
         DP_CLDLIQ_oldid    = pbuf_get_index('DP_CLDLIQ_OLD')!(pbuf_00032, lat, lon) ;
         DP_CLDICE_oldid    = pbuf_get_index('DP_CLDICE_OLD')!(pbuf_00032, lat, lon) ;
         cush_oldid         = pbuf_get_index('cush_OLD')!(lat, lon) ;
         QRS_oldid          = pbuf_get_index('QRS_OLD')!(pbuf_00032, lat, lon) ;
         QRL_oldid          = pbuf_get_index('QRL_OLD')!(pbuf_00032, lat, lon) ;
         ICIWP_oldid        = pbuf_get_index('ICIWP_OLD')!(pbuf_00032, lat, lon) ;
         ICLWP_oldid        = pbuf_get_index('ICLWP_OLD')!(pbuf_00032, lat, lon) ;
         kvm_oldid          = pbuf_get_index('kvm_OLD')!(pbuf_00033, lat, lon) ;
         turbtype_oldid     = pbuf_get_index('turbtype_OLD')!(pbuf_00033, lat, lon) ;
         smaw_oldid         = pbuf_get_index('smaw_OLD')!(pbuf_00033, lat, lon) ;
         tauresx_oldid      = pbuf_get_index('tauresx_OLD')!(lat, lon) ;
         tauresy_oldid      = pbuf_get_index('tauresy_OLD')!(lat, lon) ;
         qpert_oldid        = pbuf_get_index('qpert_OLD')!(pbuf_00033, lat, lon) ;
         T_TTEND_oldid      = pbuf_get_index('T_TTEND_OLD')!(pbuf_00032, lat, lon) ;

         ! call buffer
         call pbuf_get_field(pbuf,TEOUT_idx        ,TEOUT)!(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
         call pbuf_get_field(pbuf,DTCORE_idx       ,DTCORE)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CLDO_idx         ,CLDO)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,PRER_EVAP_idx    ,PRER_EVAP)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_T_idx         ,CC_T)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_qv_idx        ,CC_qv)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_ql_idx        ,CC_ql)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_qi_idx        ,CC_qi)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_nl_idx        ,CC_nl)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_ni_idx        ,CC_ni)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_qlst_idx      ,CC_qlst)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,am_evp_st_idx    ,am_evp_st)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,evprain_st_idx   ,evprain_st)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,evpsnow_st_idx   ,evpsnow_st)!(pbuf_00032, lat, lon)
         call pbuf_get_field(pbuf,ACPRECL_idx      ,ACPRECL)!(lat, lon) ;
         call pbuf_get_field(pbuf,ACGCME_idx       ,ACGCME)!(lat, lon) ;
         call pbuf_get_field(pbuf,ACNUM_idx        ,ACNUM)!(lat, lon) ;
         call pbuf_get_field(pbuf,RELVAR_idx       ,RELVAR)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ACCRE_ENHAN_idx  ,ACCRE_ENHAN)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,pblh_idx         ,pblh)!(lat, lon) ;
         call pbuf_get_field(pbuf,tke_idx          ,tke)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,kvh_idx          ,kvh)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,tpert_idx        ,tpert)!(lat, lon) ;
         call pbuf_get_field(pbuf,AST_idx          ,AST)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,AIST_idx         ,AIST)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ALST_idx         ,ALST)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,QIST_idx         ,QIST)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,QLST_idx         ,QLST)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CONCLD_idx       ,CONCLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CLD_idx          ,CLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,RAD_CLUBB_idx    ,RAD_CLUBB)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,WP2_nadv_idx     ,WP2_nadv)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,WP3_nadv_idx     ,WP3_nadv)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,WPTHLP_nadv_idx  ,WPTHLP_nadv)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,WPRTP_nadv_idx   ,WPRTP_nadv)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,RTPTHLP_nadv_idx ,RTPTHLP_nadv)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,RTP2_nadv_idx    ,RTP2_nadv)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,THLP2_nadv_idx   ,THLP2_nadv)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,UP2_nadv_idx     ,UP2_nadv)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,VP2_nadv_idx     ,VP2_nadv)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,UPWP_idx         ,UPWP)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,VPWP_idx         ,VPWP)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,THLM_idx         ,THLM)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,RTM_idx          ,RTM)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,UM_idx           ,UM)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,VM_idx           ,VM)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,DGNUM_idx        ,DGNUM)!(pbuf_00128, lat, lon) ;
         call pbuf_get_field(pbuf,DGNUMWET_idx     ,DGNUMWET)!(pbuf_00128, lat, lon) ;
         call pbuf_get_field(pbuf,num_c1_idx       ,num_c1)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,so4_c1_idx       ,so4_c1)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,pom_c1_idx       ,pom_c1)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,soa_c1_idx       ,soa_c1)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,bc_c1_idx        ,bc_c1)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,dst_c1_idx       ,dst_c1)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ncl_c1_idx       ,ncl_c1)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,num_c2_idx       ,num_c2)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,so4_c2_idx       ,so4_c2)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,soa_c2_idx       ,soa_c2)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ncl_c2_idx       ,ncl_c2)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,dst_c2_idx       ,dst_c2)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,num_c3_idx       ,num_c3)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,dst_c3_idx       ,dst_c3)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ncl_c3_idx       ,ncl_c3)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,so4_c3_idx       ,so4_c3)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,num_c4_idx       ,num_c4)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,pom_c4_idx       ,pom_c4)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,bc_c4_idx        ,bc_c4)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,DP_FLXPRC_idx    ,DP_FLXPRC)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,DP_FLXSNW_idx    ,DP_FLXSNW)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,DP_CLDLIQ_idx    ,DP_CLDLIQ)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,DP_CLDICE_idx    ,DP_CLDICE)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,cush_idx         ,cush)!(lat, lon) ;
         call pbuf_get_field(pbuf,QRS_idx          ,QRS)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,QRL_idx          ,QRL)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ICIWP_idx        ,ICIWP)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ICLWP_idx        ,ICLWP)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,kvm_idx          ,kvm)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,turbtype_idx     ,turbtype)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,smaw_idx         ,smaw)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,tauresx_idx      ,tauresx)!(lat, lon) ;
         call pbuf_get_field(pbuf,tauresy_idx      ,tauresy)!(lat, lon) ;
         call pbuf_get_field(pbuf,qpert_idx        ,qpert)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,T_TTEND_idx      ,T_TTEND)!(pbuf_00032, lat, lon) ;

         call pbuf_get_field(pbuf,TEOUT_oldid        ,TEOUT_OLD)!(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
         call pbuf_get_field(pbuf,DTCORE_oldid       ,DTCORE_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CLDO_oldid         ,CLDO_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,PRER_EVAP_oldid    ,PRER_EVAP_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_T_oldid         ,CC_T_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_qv_oldid        ,CC_qv_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_ql_oldid        ,CC_ql_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_qi_oldid        ,CC_qi_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_nl_oldid        ,CC_nl_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_ni_oldid        ,CC_ni_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CC_qlst_oldid      ,CC_qlst_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,am_evp_st_oldid    ,am_evp_st_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,evprain_st_oldid   ,evprain_st_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,evpsnow_st_oldid   ,evpsnow_st_OLD)!(pbuf_00032, lat, lon)
         call pbuf_get_field(pbuf,ACPRECL_oldid      ,ACPRECL_OLD)!(lat, lon) ;
         call pbuf_get_field(pbuf,ACGCME_oldid       ,ACGCME_OLD)!(lat, lon) ;
         call pbuf_get_field(pbuf,ACNUM_oldid        ,ACNUM_OLD)!(lat, lon) ;
         call pbuf_get_field(pbuf,RELVAR_oldid       ,RELVAR_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ACCRE_ENHAN_oldid  ,ACCRE_ENHAN_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,pblh_oldid         ,pblh_OLD)!(lat, lon) ;
         call pbuf_get_field(pbuf,tke_oldid          ,tke_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,kvh_oldid          ,kvh_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,tpert_oldid        ,tpert_OLD)!(lat, lon) ;
         call pbuf_get_field(pbuf,AST_oldid          ,AST_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,AIST_oldid         ,AIST_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ALST_oldid         ,ALST_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,QIST_oldid         ,QIST_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,QLST_oldid         ,QLST_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CONCLD_oldid       ,CONCLD_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,CLD_oldid          ,CLD_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,RAD_CLUBB_oldid    ,RAD_CLUBB_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,WP2_nadv_oldid     ,WP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,WP3_nadv_oldid     ,WP3_nadv_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,WPTHLP_nadv_oldid  ,WPTHLP_nadv_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,WPRTP_nadv_oldid   ,WPRTP_nadv_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,RTPTHLP_nadv_oldid ,RTPTHLP_nadv_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,RTP2_nadv_oldid    ,RTP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,THLP2_nadv_oldid   ,THLP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,UP2_nadv_oldid     ,UP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,VP2_nadv_oldid     ,VP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,UPWP_oldid         ,UPWP_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,VPWP_oldid         ,VPWP_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,THLM_oldid         ,THLM_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,RTM_oldid          ,RTM_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,UM_oldid           ,UM_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,VM_oldid           ,VM_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,DGNUM_oldid        ,DGNUM_OLD)!(pbuf_00128, lat, lon) ;
         call pbuf_get_field(pbuf,DGNUMWET_oldid     ,DGNUMWET_OLD)!(pbuf_00128, lat, lon) ;
         call pbuf_get_field(pbuf,num_c1_oldid       ,num_c1_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,so4_c1_oldid       ,so4_c1_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,pom_c1_oldid       ,pom_c1_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,soa_c1_oldid       ,soa_c1_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,bc_c1_oldid        ,bc_c1_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,dst_c1_oldid       ,dst_c1_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ncl_c1_oldid       ,ncl_c1_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,num_c2_oldid       ,num_c2_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,so4_c2_oldid       ,so4_c2_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,soa_c2_oldid       ,soa_c2_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ncl_c2_oldid       ,ncl_c2_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,dst_c2_oldid       ,dst_c2_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,num_c3_oldid       ,num_c3_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,dst_c3_oldid       ,dst_c3_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ncl_c3_oldid       ,ncl_c3_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,so4_c3_oldid       ,so4_c3_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,num_c4_oldid       ,num_c4_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,pom_c4_oldid       ,pom_c4_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,bc_c4_oldid        ,bc_c4_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,DP_FLXPRC_oldid    ,DP_FLXPRC_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,DP_FLXSNW_oldid    ,DP_FLXSNW_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,DP_CLDLIQ_oldid    ,DP_CLDLIQ_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,DP_CLDICE_oldid    ,DP_CLDICE_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,cush_oldid         ,cush_OLD)!(lat, lon) ;
         call pbuf_get_field(pbuf,QRS_oldid          ,QRS_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,QRL_oldid          ,QRL_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ICIWP_oldid        ,ICIWP_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,ICLWP_oldid        ,ICLWP_OLD)!(pbuf_00032, lat, lon) ;
         call pbuf_get_field(pbuf,kvm_oldid          ,kvm_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,turbtype_oldid     ,turbtype_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,smaw_oldid         ,smaw_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,tauresx_oldid      ,tauresx_OLD)!(lat, lon) ;
         call pbuf_get_field(pbuf,tauresy_oldid      ,tauresy_OLD)!(lat, lon) ;
         call pbuf_get_field(pbuf,qpert_oldid        ,qpert_OLD)!(pbuf_00033, lat, lon) ;
         call pbuf_get_field(pbuf,T_TTEND_oldid      ,T_TTEND_OLD)!(pbuf_00032, lat, lon) ;

         ! swap old to current
      ncol = get_ncols_p(lchnk)
      call rad_cnst_get_info(0, nmodes=nmodes)
      do i = 1, ncol 
      do k = 1, pver ! 2d

         do m = 1, nmodes ! 3d
            DGNUM_old(i,k,m) = DGNUM(i,k,m) 
            DGNUMWET_old(i,k,m) = DGNUMWET(i,k,m)   
         end do ! nmodes
         
         DTCORE_old(i,k) = DTCORE(i,k)
         CLDO_old(i,k) = CLDO(i,k)
         PRER_EVAP_old(i,k) = PRER_EVAP(i,k)
         CC_T_old(i,k) = CC_T(i,k)
         CC_qv_old(i,k) = CC_qv(i,k)
         CC_ql_old(i,k) = CC_ql(i,k)
         CC_qi_old(i,k) = CC_qi(i,k)
         CC_nl_old(i,k) = CC_nl(i,k)
         CC_ni_old(i,k) = CC_ni(i,k)
         CC_qlst_old(i,k) = CC_qlst(i,k)
         am_evp_st_old(i,k) = am_evp_st(i,k)
         evprain_st_old(i,k) = evprain_st(i,k)
         evpsnow_st_old(i,k) = evpsnow_st(i,k)
         RELVAR_old(i,k) = RELVAR(i,k)
         ACCRE_ENHAN_old(i,k) = ACCRE_ENHAN(i,k)
         AST_old(i,k) = AST(i,k)
         AIST_old(i,k) = AIST(i,k)
         ALST_old(i,k) = ALST(i,k)
         QIST_old(i,k) = QIST(i,k)
         QLST_old(i,k) = QLST(i,k)
         CONCLD_old(i,k) = CONCLD(i,k)
         CLD_old(i,k) = CLD(i,k)
         RAD_CLUBB_old(i,k) = RAD_CLUBB(i,k)
         num_c1_old(i,k) = num_c1(i,k)
         so4_c1_old(i,k) = so4_c1(i,k)
         pom_c1_old(i,k) = pom_c1(i,k)
         soa_c1_old(i,k) = soa_c1(i,k)
         bc_c1_old(i,k) = bc_c1(i,k)
         dst_c1_old(i,k) = dst_c1(i,k)
         ncl_c1_old(i,k) = ncl_c1(i,k)
         num_c2_old(i,k) = num_c2(i,k)
         so4_c2_old(i,k) = so4_c2(i,k)
         soa_c2_old(i,k) = soa_c2(i,k)
         ncl_c2_old(i,k) = ncl_c2(i,k)
         dst_c2_old(i,k) = dst_c2(i,k)
         num_c3_old(i,k) = num_c3(i,k)
         dst_c3_old(i,k) = dst_c3(i,k)
         ncl_c3_old(i,k) = ncl_c3(i,k)
         so4_c3_old(i,k) = so4_c3(i,k)
         num_c4_old(i,k) = num_c4(i,k)
         pom_c4_old(i,k) = pom_c4(i,k)
         bc_c4_old(i,k) = bc_c4(i,k)
         DP_CLDLIQ_old(i,k) = DP_CLDLIQ(i,k)
         DP_CLDICE_old(i,k) = DP_CLDICE(i,k)
         QRS_old(i,k) = QRS(i,k)
         QRL_old(i,k) = QRL(i,k)
         ICIWP_old(i,k) = ICIWP(i,k)
         ICLWP_old(i,k) = ICLWP(i,k)
         T_TTEND_old(i,k) = T_TTEND(i,k)
      end do ! k=1,pver

      do k = 1, pverp ! 2d
         tke_old(i,k) = tke(i,k)
         kvh_old(i,k) = kvh(i,k)
         WP2_nadv_old(i,k) = WP2_nadv(i,k)
         WP3_nadv_old(i,k) = WP3_nadv(i,k)
         WPTHLP_nadv_old(i,k) = WPTHLP_nadv(i,k)
         WPRTP_nadv_old(i,k) = WPRTP_nadv(i,k)
         RTPTHLP_nadv_old(i,k) = RTPTHLP_nadv(i,k)
         RTP2_nadv_old(i,k) = RTP2_nadv(i,k)
         THLP2_nadv_old(i,k) = THLP2_nadv(i,k)
         UP2_nadv_old(i,k) = UP2_nadv(i,k)
         VP2_nadv_old(i,k) = VP2_nadv(i,k)
         UPWP_old(i,k) = UPWP(i,k)
         VPWP_old(i,k) = VPWP(i,k)
         THLM_old(i,k) = THLM(i,k)
         RTM_old(i,k) = RTM(i,k)
         UM_old(i,k) = UM(i,k)
         VM_old(i,k) = VM(i,k)
         DP_FLXPRC_old(i,k) = DP_FLXPRC(i,k)
         DP_FLXSNW_old(i,k) = DP_FLXSNW(i,k)
         kvm_old(i,k) = kvm(i,k)
         turbtype_old(i,k) = turbtype(i,k)
         smaw_old(i,k) = smaw(i,k)
         qpert_old(i,k) = qpert(i,k)
      end do ! k=1,pverp
      
      TEOUT_old(i) = TEOUT(i)
      ACPRECL_old(i) = ACPRECL(i)
      ACGCME_old(i) = ACGCME(i)
      ACNUM_old(i) = ACNUM(i)
      pblh_old(i) = pblh(i)
      tpert_old(i) = tpert(i)
      cush_old(i) = cush(i)
      tauresx_old(i) = tauresx(i)
      tauresy_old(i) = tauresy(i)
      end do ! i=1,ncol
      end do ! lchunk

   end if ! mod(sec_spec,21600)==0  .AND. .NOT. do_restart 

   if ( mod(sec_spec,21600)==10800 .AND. do_restart ) then

      do_restart=.FALSE.

      !filein = "atm_in"
      !fname_pbuf_cam = interpret_filename_spec( rfilename_spec_cam, yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec=sec_spec-10800 )
      !locfn= './' // trim(fname_pbuf_cam)
      do lchnk=begchunk,endchunk
      pbuf=> pbuf_get_chunk(pbuf2d,lchnk)

      if(masterproc) then
         print *, 'swap pbuf new to old ', sec_spec
      end if

      ! get index of buffer vars
      TEOUT_idx        = pbuf_get_index('TEOUT')!(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
      DTCORE_idx       = pbuf_get_index('DTCORE')!(pbuf_00032, lat, lon) ;
      CLDO_idx         = pbuf_get_index('CLDO')!(pbuf_00032, lat, lon) ;
      PRER_EVAP_idx    = pbuf_get_index('PRER_EVAP')!(pbuf_00032, lat, lon) ;
      CC_T_idx         = pbuf_get_index('CC_T')!(pbuf_00032, lat, lon) ;
      CC_qv_idx        = pbuf_get_index('CC_qv')!(pbuf_00032, lat, lon) ;
      CC_ql_idx        = pbuf_get_index('CC_ql')!(pbuf_00032, lat, lon) ;
      CC_qi_idx        = pbuf_get_index('CC_qi')!(pbuf_00032, lat, lon) ;
      CC_nl_idx        = pbuf_get_index('CC_nl')!(pbuf_00032, lat, lon) ;
      CC_ni_idx        = pbuf_get_index('CC_ni')!(pbuf_00032, lat, lon) ;
      CC_qlst_idx      = pbuf_get_index('CC_qlst')!(pbuf_00032, lat, lon) ;
      am_evp_st_idx    = pbuf_get_index('am_evp_st')!(pbuf_00032, lat, lon) ;
      evprain_st_idx   = pbuf_get_index('evprain_st')!(pbuf_00032, lat, lon) ;
      evpsnow_st_idx   = pbuf_get_index('evpsnow_st')!(pbuf_00032, lat, lon)
      ACPRECL_idx      = pbuf_get_index('ACPRECL')!(lat, lon) ;
      ACGCME_idx       = pbuf_get_index('ACGCME')!(lat, lon) ;
      ACNUM_idx        = pbuf_get_index('ACNUM')!(lat, lon) ;
      RELVAR_idx       = pbuf_get_index('RELVAR')!(pbuf_00032, lat, lon) ;
      ACCRE_ENHAN_idx  = pbuf_get_index('ACCRE_ENHAN')!(pbuf_00032, lat, lon) ;
      pblh_idx         = pbuf_get_index('pblh')!(lat, lon) ;
      tke_idx          = pbuf_get_index('tke')!(pbuf_00033, lat, lon) ;
      kvh_idx          = pbuf_get_index('kvh')!(pbuf_00033, lat, lon) ;
      tpert_idx        = pbuf_get_index('tpert')!(lat, lon) ;
      AST_idx          = pbuf_get_index('AST')!(pbuf_00032, lat, lon) ;
      AIST_idx         = pbuf_get_index('AIST')!(pbuf_00032, lat, lon) ;
      ALST_idx         = pbuf_get_index('ALST')!(pbuf_00032, lat, lon) ;
      QIST_idx         = pbuf_get_index('QIST')!(pbuf_00032, lat, lon) ;
      QLST_idx         = pbuf_get_index('QLST')!(pbuf_00032, lat, lon) ;
      CONCLD_idx       = pbuf_get_index('CONCLD')!(pbuf_00032, lat, lon) ;
      CLD_idx          = pbuf_get_index('CLD')!(pbuf_00032, lat, lon) ;
      RAD_CLUBB_idx    = pbuf_get_index('RAD_CLUBB')!(pbuf_00032, lat, lon) ;
      WP2_nadv_idx     = pbuf_get_index('WP2_nadv')!(pbuf_00033, lat, lon) ;
      WP3_nadv_idx     = pbuf_get_index('WP3_nadv')!(pbuf_00033, lat, lon) ;
      WPTHLP_nadv_idx  = pbuf_get_index('WPTHLP_nadv')!(pbuf_00033, lat, lon) ;
      WPRTP_nadv_idx   = pbuf_get_index('WPRTP_nadv')!(pbuf_00033, lat, lon) ;
      RTPTHLP_nadv_idx = pbuf_get_index('RTPTHLP_nadv')!(pbuf_00033, lat, lon) ;
      RTP2_nadv_idx    = pbuf_get_index('RTP2_nadv')!(pbuf_00033, lat, lon) ;
      THLP2_nadv_idx   = pbuf_get_index('THLP2_nadv')!(pbuf_00033, lat, lon) ;
      UP2_nadv_idx     = pbuf_get_index('UP2_nadv')!(pbuf_00033, lat, lon) ;
      VP2_nadv_idx     = pbuf_get_index('VP2_nadv')!(pbuf_00033, lat, lon) ;
      UPWP_idx         = pbuf_get_index('UPWP')!(pbuf_00033, lat, lon) ;
      VPWP_idx         = pbuf_get_index('VPWP')!(pbuf_00033, lat, lon) ;
      THLM_idx         = pbuf_get_index('THLM')!(pbuf_00033, lat, lon) ;
      RTM_idx          = pbuf_get_index('RTM')!(pbuf_00033, lat, lon) ;
      UM_idx           = pbuf_get_index('UM')!(pbuf_00033, lat, lon) ;
      VM_idx           = pbuf_get_index('VM')!(pbuf_00033, lat, lon) ;
      DGNUM_idx        = pbuf_get_index('DGNUM')!(pbuf_00128, lat, lon) ;
      DGNUMWET_idx     = pbuf_get_index('DGNUMWET')!(pbuf_00128, lat, lon) ;
      num_c1_idx       = pbuf_get_index('num_c1')!(pbuf_00032, lat, lon) ;
      so4_c1_idx       = pbuf_get_index('so4_c1')!(pbuf_00032, lat, lon) ;
      pom_c1_idx       = pbuf_get_index('pom_c1')!(pbuf_00032, lat, lon) ;
      soa_c1_idx       = pbuf_get_index('soa_c1')!(pbuf_00032, lat, lon) ;
      bc_c1_idx        = pbuf_get_index('bc_c1')!(pbuf_00032, lat, lon) ;
      dst_c1_idx       = pbuf_get_index('dst_c1')!(pbuf_00032, lat, lon) ;
      ncl_c1_idx       = pbuf_get_index('ncl_c1')!(pbuf_00032, lat, lon) ;
      num_c2_idx       = pbuf_get_index('num_c2')!(pbuf_00032, lat, lon) ;
      so4_c2_idx       = pbuf_get_index('so4_c2')!(pbuf_00032, lat, lon) ;
      soa_c2_idx       = pbuf_get_index('soa_c2')!(pbuf_00032, lat, lon) ;
      ncl_c2_idx       = pbuf_get_index('ncl_c2')!(pbuf_00032, lat, lon) ;
      dst_c2_idx       = pbuf_get_index('dst_c2')!(pbuf_00032, lat, lon) ;
      num_c3_idx       = pbuf_get_index('num_c3')!(pbuf_00032, lat, lon) ;
      dst_c3_idx       = pbuf_get_index('dst_c3')!(pbuf_00032, lat, lon) ;
      ncl_c3_idx       = pbuf_get_index('ncl_c3')!(pbuf_00032, lat, lon) ;
      so4_c3_idx       = pbuf_get_index('so4_c3')!(pbuf_00032, lat, lon) ;
      num_c4_idx       = pbuf_get_index('num_c4')!(pbuf_00032, lat, lon) ;
      pom_c4_idx       = pbuf_get_index('pom_c4')!(pbuf_00032, lat, lon) ;
      bc_c4_idx        = pbuf_get_index('bc_c4')!(pbuf_00032, lat, lon) ;
      DP_FLXPRC_idx    = pbuf_get_index('DP_FLXPRC')!(pbuf_00033, lat, lon) ;
      DP_FLXSNW_idx    = pbuf_get_index('DP_FLXSNW')!(pbuf_00033, lat, lon) ;
      DP_CLDLIQ_idx    = pbuf_get_index('DP_CLDLIQ')!(pbuf_00032, lat, lon) ;
      DP_CLDICE_idx    = pbuf_get_index('DP_CLDICE')!(pbuf_00032, lat, lon) ;
      cush_idx         = pbuf_get_index('cush')!(lat, lon) ;
      QRS_idx          = pbuf_get_index('QRS')!(pbuf_00032, lat, lon) ;
      QRL_idx          = pbuf_get_index('QRL')!(pbuf_00032, lat, lon) ;
      ICIWP_idx        = pbuf_get_index('ICIWP')!(pbuf_00032, lat, lon) ;
      ICLWP_idx        = pbuf_get_index('ICLWP')!(pbuf_00032, lat, lon) ;
      kvm_idx          = pbuf_get_index('kvm')!(pbuf_00033, lat, lon) ;
      turbtype_idx     = pbuf_get_index('turbtype')!(pbuf_00033, lat, lon) ;
      smaw_idx         = pbuf_get_index('smaw')!(pbuf_00033, lat, lon) ;
      tauresx_idx      = pbuf_get_index('tauresx')!(lat, lon) ;
      tauresy_idx      = pbuf_get_index('tauresy')!(lat, lon) ;
      qpert_idx        = pbuf_get_index('qpert')!(pbuf_00033, lat, lon) ;
      T_TTEND_idx      = pbuf_get_index('T_TTEND')!(pbuf_00032, lat, lon) ;

      TEOUT_oldid        = pbuf_get_index('TEOUT_OLD')!(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
      DTCORE_oldid       = pbuf_get_index('DTCORE_OLD')!(pbuf_00032, lat, lon) ;
      CLDO_oldid         = pbuf_get_index('CLDO_OLD')!(pbuf_00032, lat, lon) ;
      PRER_EVAP_oldid    = pbuf_get_index('PRER_EVAP_OLD')!(pbuf_00032, lat, lon) ;
      CC_T_oldid         = pbuf_get_index('CC_T_OLD')!(pbuf_00032, lat, lon) ;
      CC_qv_oldid        = pbuf_get_index('CC_qv_OLD')!(pbuf_00032, lat, lon) ;
      CC_ql_oldid        = pbuf_get_index('CC_ql_OLD')!(pbuf_00032, lat, lon) ;
      CC_qi_oldid        = pbuf_get_index('CC_qi_OLD')!(pbuf_00032, lat, lon) ;
      CC_nl_oldid        = pbuf_get_index('CC_nl_OLD')!(pbuf_00032, lat, lon) ;
      CC_ni_oldid        = pbuf_get_index('CC_ni_OLD')!(pbuf_00032, lat, lon) ;
      CC_qlst_oldid      = pbuf_get_index('CC_qlst_OLD')!(pbuf_00032, lat, lon) ;
      am_evp_st_oldid    = pbuf_get_index('am_evp_st_OLD')!(pbuf_00032, lat, lon) ;
      evprain_st_oldid   = pbuf_get_index('evprain_st_OLD')!(pbuf_00032, lat, lon) ;
      evpsnow_st_oldid   = pbuf_get_index('evpsnow_st_OLD')!(pbuf_00032, lat, lon)
      ACPRECL_oldid      = pbuf_get_index('ACPRECL_OLD')!(lat, lon) ;
      ACGCME_oldid       = pbuf_get_index('ACGCME_OLD')!(lat, lon) ;
      ACNUM_oldid        = pbuf_get_index('ACNUM_OLD')!(lat, lon) ;
      RELVAR_oldid       = pbuf_get_index('RELVAR_OLD')!(pbuf_00032, lat, lon) ;
      ACCRE_ENHAN_oldid  = pbuf_get_index('ACCRE_ENHAN_OLD')!(pbuf_00032, lat, lon) ;
      pblh_oldid         = pbuf_get_index('pblh_OLD')!(lat, lon) ;
      tke_oldid          = pbuf_get_index('tke_OLD')!(pbuf_00033, lat, lon) ;
      kvh_oldid          = pbuf_get_index('kvh_OLD')!(pbuf_00033, lat, lon) ;
      tpert_oldid        = pbuf_get_index('tpert_OLD')!(lat, lon) ;
      AST_oldid          = pbuf_get_index('AST_OLD')!(pbuf_00032, lat, lon) ;
      AIST_oldid         = pbuf_get_index('AIST_OLD')!(pbuf_00032, lat, lon) ;
      ALST_oldid         = pbuf_get_index('ALST_OLD')!(pbuf_00032, lat, lon) ;
      QIST_oldid         = pbuf_get_index('QIST_OLD')!(pbuf_00032, lat, lon) ;
      QLST_oldid         = pbuf_get_index('QLST_OLD')!(pbuf_00032, lat, lon) ;
      CONCLD_oldid       = pbuf_get_index('CONCLD_OLD')!(pbuf_00032, lat, lon) ;
      CLD_oldid          = pbuf_get_index('CLD_OLD')!(pbuf_00032, lat, lon) ;
      RAD_CLUBB_oldid    = pbuf_get_index('RAD_CLUBB_OLD')!(pbuf_00032, lat, lon) ;
      WP2_nadv_oldid     = pbuf_get_index('WP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
      WP3_nadv_oldid     = pbuf_get_index('WP3_nadv_OLD')!(pbuf_00033, lat, lon) ;
      WPTHLP_nadv_oldid  = pbuf_get_index('WPTHLP_nadv_OLD')!(pbuf_00033, lat, lon) ;
      WPRTP_nadv_oldid   = pbuf_get_index('WPRTP_nadv_OLD')!(pbuf_00033, lat, lon) ;
      RTPTHLP_nadv_oldid = pbuf_get_index('RTPTHLP_nadv_OLD')!(pbuf_00033, lat, lon) ;
      RTP2_nadv_oldid    = pbuf_get_index('RTP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
      THLP2_nadv_oldid   = pbuf_get_index('THLP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
      UP2_nadv_oldid     = pbuf_get_index('UP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
      VP2_nadv_oldid     = pbuf_get_index('VP2_nadv_OLD')!(pbuf_00033, lat, lon) ;
      UPWP_oldid         = pbuf_get_index('UPWP_OLD')!(pbuf_00033, lat, lon) ;
      VPWP_oldid         = pbuf_get_index('VPWP_OLD')!(pbuf_00033, lat, lon) ;
      THLM_oldid         = pbuf_get_index('THLM_OLD')!(pbuf_00033, lat, lon) ;
      RTM_oldid          = pbuf_get_index('RTM_OLD')!(pbuf_00033, lat, lon) ;
      UM_oldid           = pbuf_get_index('UM_OLD')!(pbuf_00033, lat, lon) ;
      VM_oldid           = pbuf_get_index('VM_OLD')!(pbuf_00033, lat, lon) ;
      DGNUM_oldid        = pbuf_get_index('DGNUM_OLD')!(pbuf_00128, lat, lon) ;
      DGNUMWET_oldid     = pbuf_get_index('DGNUMWET_OLD')!(pbuf_00128, lat, lon) ;
      num_c1_oldid       = pbuf_get_index('num_c1_OLD')!(pbuf_00032, lat, lon) ;
      so4_c1_oldid       = pbuf_get_index('so4_c1_OLD')!(pbuf_00032, lat, lon) ;
      pom_c1_oldid       = pbuf_get_index('pom_c1_OLD')!(pbuf_00032, lat, lon) ;
      soa_c1_oldid       = pbuf_get_index('soa_c1_OLD')!(pbuf_00032, lat, lon) ;
      bc_c1_oldid        = pbuf_get_index('bc_c1_OLD')!(pbuf_00032, lat, lon) ;
      dst_c1_oldid       = pbuf_get_index('dst_c1_OLD')!(pbuf_00032, lat, lon) ;
      ncl_c1_oldid       = pbuf_get_index('ncl_c1_OLD')!(pbuf_00032, lat, lon) ;
      num_c2_oldid       = pbuf_get_index('num_c2_OLD')!(pbuf_00032, lat, lon) ;
      so4_c2_oldid       = pbuf_get_index('so4_c2_OLD')!(pbuf_00032, lat, lon) ;
      soa_c2_oldid       = pbuf_get_index('soa_c2_OLD')!(pbuf_00032, lat, lon) ;
      ncl_c2_oldid       = pbuf_get_index('ncl_c2_OLD')!(pbuf_00032, lat, lon) ;
      dst_c2_oldid       = pbuf_get_index('dst_c2_OLD')!(pbuf_00032, lat, lon) ;
      num_c3_oldid       = pbuf_get_index('num_c3_OLD')!(pbuf_00032, lat, lon) ;
      dst_c3_oldid       = pbuf_get_index('dst_c3_OLD')!(pbuf_00032, lat, lon) ;
      ncl_c3_oldid       = pbuf_get_index('ncl_c3_OLD')!(pbuf_00032, lat, lon) ;
      so4_c3_oldid       = pbuf_get_index('so4_c3_OLD')!(pbuf_00032, lat, lon) ;
      num_c4_oldid       = pbuf_get_index('num_c4_OLD')!(pbuf_00032, lat, lon) ;
      pom_c4_oldid       = pbuf_get_index('pom_c4_OLD')!(pbuf_00032, lat, lon) ;
      bc_c4_oldid        = pbuf_get_index('bc_c4_OLD')!(pbuf_00032, lat, lon) ;
      DP_FLXPRC_oldid    = pbuf_get_index('DP_FLXPRC_OLD')!(pbuf_00033, lat, lon) ;
      DP_FLXSNW_oldid    = pbuf_get_index('DP_FLXSNW_OLD')!(pbuf_00033, lat, lon) ;
      DP_CLDLIQ_oldid    = pbuf_get_index('DP_CLDLIQ_OLD')!(pbuf_00032, lat, lon) ;
      DP_CLDICE_oldid    = pbuf_get_index('DP_CLDICE_OLD')!(pbuf_00032, lat, lon) ;
      cush_oldid         = pbuf_get_index('cush_OLD')!(lat, lon) ;
      QRS_oldid          = pbuf_get_index('QRS_OLD')!(pbuf_00032, lat, lon) ;
      QRL_oldid          = pbuf_get_index('QRL_OLD')!(pbuf_00032, lat, lon) ;
      ICIWP_oldid        = pbuf_get_index('ICIWP_OLD')!(pbuf_00032, lat, lon) ;
      ICLWP_oldid        = pbuf_get_index('ICLWP_OLD')!(pbuf_00032, lat, lon) ;
      kvm_oldid          = pbuf_get_index('kvm_OLD')!(pbuf_00033, lat, lon) ;
      turbtype_oldid     = pbuf_get_index('turbtype_OLD')!(pbuf_00033, lat, lon) ;
      smaw_oldid         = pbuf_get_index('smaw_OLD')!(pbuf_00033, lat, lon) ;
      tauresx_oldid      = pbuf_get_index('tauresx_OLD')!(lat, lon) ;
      tauresy_oldid      = pbuf_get_index('tauresy_OLD')!(lat, lon) ;
      qpert_oldid        = pbuf_get_index('qpert_OLD')!(pbuf_00033, lat, lon) ;
      T_TTEND_oldid      = pbuf_get_index('T_TTEND_OLD')!(pbuf_00032, lat, lon) ;

      ! call buffer
      call pbuf_get_field(pbuf,TEOUT_idx        ,TEOUT)!(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
      call pbuf_get_field(pbuf,DTCORE_idx       ,DTCORE)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CLDO_idx         ,CLDO)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,PRER_EVAP_idx    ,PRER_EVAP)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_T_idx         ,CC_T)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_qv_idx        ,CC_qv)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_ql_idx        ,CC_ql)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_qi_idx        ,CC_qi)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_nl_idx        ,CC_nl)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_ni_idx        ,CC_ni)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_qlst_idx      ,CC_qlst)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,am_evp_st_idx    ,am_evp_st)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,evprain_st_idx   ,evprain_st)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,evpsnow_st_idx   ,evpsnow_st)!(pbuf_00032, lat, lon)
      call pbuf_get_field(pbuf,ACPRECL_idx      ,ACPRECL)!(lat, lon) ;
      call pbuf_get_field(pbuf,ACGCME_idx       ,ACGCME)!(lat, lon) ;
      call pbuf_get_field(pbuf,ACNUM_idx        ,ACNUM)!(lat, lon) ;
      call pbuf_get_field(pbuf,RELVAR_idx       ,RELVAR)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ACCRE_ENHAN_idx  ,ACCRE_ENHAN)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,pblh_idx         ,pblh)!(lat, lon) ;
      call pbuf_get_field(pbuf,tke_idx          ,tke)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,kvh_idx          ,kvh)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,tpert_idx        ,tpert)!(lat, lon) ;
      call pbuf_get_field(pbuf,AST_idx          ,AST)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,AIST_idx         ,AIST)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ALST_idx         ,ALST)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,QIST_idx         ,QIST)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,QLST_idx         ,QLST)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CONCLD_idx       ,CONCLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CLD_idx          ,CLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,RAD_CLUBB_idx    ,RAD_CLUBB)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,WP2_nadv_idx     ,WP2_nadv)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,WP3_nadv_idx     ,WP3_nadv)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,WPTHLP_nadv_idx  ,WPTHLP_nadv)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,WPRTP_nadv_idx   ,WPRTP_nadv)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,RTPTHLP_nadv_idx ,RTPTHLP_nadv)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,RTP2_nadv_idx    ,RTP2_nadv)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,THLP2_nadv_idx   ,THLP2_nadv)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,UP2_nadv_idx     ,UP2_nadv)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,VP2_nadv_idx     ,VP2_nadv)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,UPWP_idx         ,UPWP)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,VPWP_idx         ,VPWP)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,THLM_idx         ,THLM)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,RTM_idx          ,RTM)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,UM_idx           ,UM)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,VM_idx           ,VM)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,DGNUM_idx        ,DGNUM)!(pbuf_00128, lat, lon) ;
      call pbuf_get_field(pbuf,DGNUMWET_idx     ,DGNUMWET)!(pbuf_00128, lat, lon) ;
      call pbuf_get_field(pbuf,num_c1_idx       ,num_c1)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,so4_c1_idx       ,so4_c1)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,pom_c1_idx       ,pom_c1)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,soa_c1_idx       ,soa_c1)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,bc_c1_idx        ,bc_c1)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,dst_c1_idx       ,dst_c1)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ncl_c1_idx       ,ncl_c1)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,num_c2_idx       ,num_c2)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,so4_c2_idx       ,so4_c2)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,soa_c2_idx       ,soa_c2)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ncl_c2_idx       ,ncl_c2)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,dst_c2_idx       ,dst_c2)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,num_c3_idx       ,num_c3)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,dst_c3_idx       ,dst_c3)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ncl_c3_idx       ,ncl_c3)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,so4_c3_idx       ,so4_c3)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,num_c4_idx       ,num_c4)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,pom_c4_idx       ,pom_c4)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,bc_c4_idx        ,bc_c4)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,DP_FLXPRC_idx    ,DP_FLXPRC)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,DP_FLXSNW_idx    ,DP_FLXSNW)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,DP_CLDLIQ_idx    ,DP_CLDLIQ)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,DP_CLDICE_idx    ,DP_CLDICE)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,cush_idx         ,cush)!(lat, lon) ;
      call pbuf_get_field(pbuf,QRS_idx          ,QRS)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,QRL_idx          ,QRL)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ICIWP_idx        ,ICIWP)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ICLWP_idx        ,ICLWP)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,kvm_idx          ,kvm)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,turbtype_idx     ,turbtype)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,smaw_idx         ,smaw)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,tauresx_idx      ,tauresx)!(lat, lon) ;
      call pbuf_get_field(pbuf,tauresy_idx      ,tauresy)!(lat, lon) ;
      call pbuf_get_field(pbuf,qpert_idx        ,qpert)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,T_TTEND_idx      ,T_TTEND)!(pbuf_00032, lat, lon) ;

      call pbuf_get_field(pbuf,TEOUT_oldid        ,TEOUT_OLD)!(lat, lon) ; ! pbuf vars from cam.r. output - sweidman
      call pbuf_get_field(pbuf,DTCORE_oldid       ,DTCORE_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CLDO_oldid         ,CLDO_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,PRER_EVAP_oldid    ,PRER_EVAP_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_T_oldid         ,CC_T_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_qv_oldid        ,CC_qv_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_ql_oldid        ,CC_ql_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_qi_oldid        ,CC_qi_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_nl_oldid        ,CC_nl_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_ni_oldid        ,CC_ni_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CC_qlst_oldid      ,CC_qlst_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,am_evp_st_oldid    ,am_evp_st_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,evprain_st_oldid   ,evprain_st_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,evpsnow_st_oldid   ,evpsnow_st_OLD)!(pbuf_00032, lat, lon)
      call pbuf_get_field(pbuf,ACPRECL_oldid      ,ACPRECL_OLD)!(lat, lon) ;
      call pbuf_get_field(pbuf,ACGCME_oldid       ,ACGCME_OLD)!(lat, lon) ;
      call pbuf_get_field(pbuf,ACNUM_oldid        ,ACNUM_OLD)!(lat, lon) ;
      call pbuf_get_field(pbuf,RELVAR_oldid       ,RELVAR_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ACCRE_ENHAN_oldid  ,ACCRE_ENHAN_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,pblh_oldid         ,pblh_OLD)!(lat, lon) ;
      call pbuf_get_field(pbuf,tke_oldid          ,tke_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,kvh_oldid          ,kvh_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,tpert_oldid        ,tpert_OLD)!(lat, lon) ;
      call pbuf_get_field(pbuf,AST_oldid          ,AST_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,AIST_oldid         ,AIST_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ALST_oldid         ,ALST_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,QIST_oldid         ,QIST_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,QLST_oldid         ,QLST_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CONCLD_oldid       ,CONCLD_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,CLD_oldid          ,CLD_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,RAD_CLUBB_oldid    ,RAD_CLUBB_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,WP2_nadv_oldid     ,WP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,WP3_nadv_oldid     ,WP3_nadv_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,WPTHLP_nadv_oldid  ,WPTHLP_nadv_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,WPRTP_nadv_oldid   ,WPRTP_nadv_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,RTPTHLP_nadv_oldid ,RTPTHLP_nadv_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,RTP2_nadv_oldid    ,RTP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,THLP2_nadv_oldid   ,THLP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,UP2_nadv_oldid     ,UP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,VP2_nadv_oldid     ,VP2_nadv_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,UPWP_oldid         ,UPWP_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,VPWP_oldid         ,VPWP_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,THLM_oldid         ,THLM_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,RTM_oldid          ,RTM_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,UM_oldid           ,UM_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,VM_oldid           ,VM_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,DGNUM_oldid        ,DGNUM_OLD)!(pbuf_00128, lat, lon) ;
      call pbuf_get_field(pbuf,DGNUMWET_oldid     ,DGNUMWET_OLD)!(pbuf_00128, lat, lon) ;
      call pbuf_get_field(pbuf,num_c1_oldid       ,num_c1_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,so4_c1_oldid       ,so4_c1_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,pom_c1_oldid       ,pom_c1_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,soa_c1_oldid       ,soa_c1_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,bc_c1_oldid        ,bc_c1_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,dst_c1_oldid       ,dst_c1_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ncl_c1_oldid       ,ncl_c1_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,num_c2_oldid       ,num_c2_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,so4_c2_oldid       ,so4_c2_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,soa_c2_oldid       ,soa_c2_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ncl_c2_oldid       ,ncl_c2_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,dst_c2_oldid       ,dst_c2_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,num_c3_oldid       ,num_c3_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,dst_c3_oldid       ,dst_c3_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ncl_c3_oldid       ,ncl_c3_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,so4_c3_oldid       ,so4_c3_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,num_c4_oldid       ,num_c4_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,pom_c4_oldid       ,pom_c4_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,bc_c4_oldid        ,bc_c4_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,DP_FLXPRC_oldid    ,DP_FLXPRC_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,DP_FLXSNW_oldid    ,DP_FLXSNW_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,DP_CLDLIQ_oldid    ,DP_CLDLIQ_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,DP_CLDICE_oldid    ,DP_CLDICE_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,cush_oldid         ,cush_OLD)!(lat, lon) ;
      call pbuf_get_field(pbuf,QRS_oldid          ,QRS_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,QRL_oldid          ,QRL_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ICIWP_oldid        ,ICIWP_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,ICLWP_oldid        ,ICLWP_OLD)!(pbuf_00032, lat, lon) ;
      call pbuf_get_field(pbuf,kvm_oldid          ,kvm_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,turbtype_oldid     ,turbtype_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,smaw_oldid         ,smaw_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,tauresx_oldid      ,tauresx_OLD)!(lat, lon) ;
      call pbuf_get_field(pbuf,tauresy_oldid      ,tauresy_OLD)!(lat, lon) ;
      call pbuf_get_field(pbuf,qpert_oldid        ,qpert_OLD)!(pbuf_00033, lat, lon) ;
      call pbuf_get_field(pbuf,T_TTEND_oldid      ,T_TTEND_OLD)!(pbuf_00032, lat, lon) ;

      ! swap normal to old
      ncol = get_ncols_p(lchnk)
      call rad_cnst_get_info(0, nmodes=nmodes) 
      do i = 1, ncol 
      do k = 1, pver ! 2d
         do m = 1, nmodes
            DGNUM(i,k,m) = DGNUM_old(i,k,m)       !(pbuf_00128, lat, lon) ;
            DGNUMWET(i,k,m) = DGNUMWET_old(i,k,m)   
         end do ! nmodes
         
         DTCORE(i,k) = DTCORE_old(i,k)
         CLDO(i,k) = CLDO_old(i,k)
         PRER_EVAP(i,k) = PRER_EVAP_old(i,k)
         CC_T(i,k) = CC_T_old(i,k)
         CC_qv(i,k) = CC_qv_old(i,k)
         CC_ql(i,k) = CC_ql_old(i,k)
         CC_qi(i,k) = CC_qi_old(i,k)
         CC_nl(i,k) = CC_nl_old(i,k)
         CC_ni(i,k) = CC_ni_old(i,k)
         CC_qlst(i,k) = CC_qlst_old(i,k)
         am_evp_st(i,k) = am_evp_st_old(i,k)
         evprain_st(i,k) = evprain_st_old(i,k)
         evpsnow_st(i,k) = evpsnow_st_old(i,k)
         RELVAR(i,k) = RELVAR_old(i,k)
         ACCRE_ENHAN(i,k) = ACCRE_ENHAN_old(i,k)
         AST(i,k) = AST_old(i,k)
         AIST(i,k) = AIST_old(i,k)
         ALST(i,k) = ALST_old(i,k)
         QIST(i,k) = QIST_old(i,k)
         QLST(i,k) = QLST_old(i,k)
         CONCLD(i,k) = CONCLD_old(i,k)
         CLD(i,k) = CLD_old(i,k)
         RAD_CLUBB(i,k) = RAD_CLUBB_old(i,k)
         num_c1(i,k) = num_c1_old(i,k)
         so4_c1(i,k) = so4_c1_old(i,k)
         pom_c1(i,k) = pom_c1_old(i,k)
         soa_c1(i,k) = soa_c1_old(i,k)
         bc_c1(i,k) = bc_c1_old(i,k)
         dst_c1(i,k) = dst_c1_old(i,k)
         ncl_c1(i,k) = ncl_c1_old(i,k)
         num_c2(i,k) = num_c2_old(i,k)
         so4_c2(i,k) = so4_c2_old(i,k)
         soa_c2(i,k) = soa_c2_old(i,k)
         ncl_c2(i,k) = ncl_c2_old(i,k)
         dst_c2(i,k) = dst_c2_old(i,k)
         num_c3(i,k) = num_c3_old(i,k)
         dst_c3(i,k) = dst_c3_old(i,k)
         ncl_c3(i,k) = ncl_c3_old(i,k)
         so4_c3(i,k) = so4_c3_old(i,k)
         num_c4(i,k) = num_c4_old(i,k)
         pom_c4(i,k) = pom_c4_old(i,k)
         bc_c4(i,k) = bc_c4_old(i,k)
         DP_CLDLIQ(i,k) = DP_CLDLIQ_old(i,k)
         DP_CLDICE(i,k) = DP_CLDICE_old(i,k)
         QRS(i,k) = QRS_old(i,k)
         QRL(i,k) = QRL_old(i,k)
         ICIWP(i,k) = ICIWP_old(i,k)
         ICLWP(i,k) = ICLWP_old(i,k)
         T_TTEND(i,k) = T_TTEND_old(i,k)

      end do ! k=1,pver

      do k = 1, pverp
         tke(i,k) = tke_old(i,k)
         kvh(i,k) = kvh_old(i,k)
         WP2_nadv(i,k) = WP2_nadv_old(i,k)
         WP3_nadv(i,k) = WP3_nadv_old(i,k)
         WPTHLP_nadv(i,k) = WPTHLP_nadv_old(i,k)
         WPRTP_nadv(i,k) = WPRTP_nadv_old(i,k)
         RTPTHLP_nadv(i,k) = RTPTHLP_nadv_old(i,k)
         RTP2_nadv(i,k) = RTP2_nadv_old(i,k)
         THLP2_nadv(i,k) = THLP2_nadv_old(i,k)
         UP2_nadv(i,k) = UP2_nadv_old(i,k)
         VP2_nadv(i,k) = VP2_nadv_old(i,k)
         UPWP(i,k) = UPWP_old(i,k)
         VPWP(i,k) = VPWP_old(i,k)
         THLM(i,k) = THLM_old(i,k)
         RTM(i,k) = RTM_old(i,k)
         UM(i,k) = UM_old(i,k)
         VM(i,k) = VM_old(i,k)
         DP_FLXPRC(i,k) = DP_FLXPRC_old(i,k)
         DP_FLXSNW(i,k) = DP_FLXSNW_old(i,k)
         kvm(i,k) = kvm_old(i,k)
         turbtype(i,k) = turbtype_old(i,k)
         smaw(i,k) = smaw_old(i,k)
         qpert(i,k) = qpert_old(i,k)
      end do ! k=1,pverp

      TEOUT(i) = TEOUT_old(i)
      ACPRECL(i) = ACPRECL_old(i)
      ACGCME(i) = ACGCME_old(i)
      ACNUM(i) = ACNUM_old(i)
      pblh(i) = pblh_old(i)
      tpert(i) = tpert_old(i)
      cush(i) = cush_old(i)
      tauresx(i) = tauresx_old(i)
      tauresy(i) = tauresy_old(i)
      end do ! i=1,ncol
      end do ! lchunk

      ! zero out tendencies
      do lchnk = begchunk,endchunk
         ncol = get_ncols_p(lchnk)
         do k = 1, pver
            do i = 1, ncol
                phys_tend(lchnk)%dudt (i,k)=0
                phys_tend(lchnk)%dvdt (i,k)=0
            end do
         end do
      end do
      do lchnk = begchunk,endchunk
         ncol = get_ncols_p(lchnk)
         do i = 1, ncol
                phys_tend(lchnk)%flx_net (i)=0
                phys_tend(lchnk)%te_tnd (i)=0
                phys_tend(lchnk)%tw_tnd (i)=0
         end do
      end do
   end if !  mod(sec_spec,21600)==10800 .AND. do_restart

   call t_startf ('cam_run4_wrapup')
   call wrapup(rstwr, nlend)
   call t_stopf  ('cam_run4_wrapup')

   call qneg_print_summary(is_last_step())

   call shr_sys_flush(iulog)

end subroutine cam_run4

!
!-----------------------------------------------------------------------
!

subroutine cam_final( cam_out, cam_in )
!-----------------------------------------------------------------------
!
! Purpose:  CAM finalization.
!
!-----------------------------------------------------------------------
   use stepon,           only: stepon_final
   use physpkg,          only: phys_final
   use cam_initfiles,    only: cam_initfiles_close
   use camsrfexch,       only: atm2hub_deallocate, hub2atm_deallocate
   use ionosphere_interface, only: ionosphere_final

   !
   ! Arguments
   !
   type(cam_out_t), pointer :: cam_out(:) ! Output from CAM to surface
   type(cam_in_t),  pointer :: cam_in(:)   ! Input from merged surface to CAM

   ! Local variables
   integer :: nstep           ! Current timestep number.
   !-----------------------------------------------------------------------

   call phys_final( phys_state, phys_tend , pbuf2d)
   call stepon_final(dyn_in, dyn_out)
   call ionosphere_final()

   if (initial_run) then
      call cam_initfiles_close()
   end if

   call hub2atm_deallocate(cam_in)
   call atm2hub_deallocate(cam_out)

   ! This flush attempts to ensure that asynchronous diagnostic prints from all
   ! processes do not get mixed up with the "END OF MODEL RUN" message printed
   ! by masterproc below.  The test-model script searches for this message in the
   ! output log to figure out if CAM completed successfully.
   call shr_sys_flush( 0 )       ! Flush all output to standard error
   call shr_sys_flush( iulog )   ! Flush all output to the CAM log file

   if (masterproc) then
      nstep = get_nstep()
      write(iulog,9300) nstep-1,nstep
9300  format (//'Number of completed timesteps:',i6,/,'Time step ',i6, &
                ' partially done to provide convectively adjusted and ', &
                'time filtered values for history tape.')
      write(iulog,*)' '
      write(iulog,*)'******* END OF MODEL RUN *******'
   end if

end subroutine cam_final

!-----------------------------------------------------------------------

end module cam_comp
