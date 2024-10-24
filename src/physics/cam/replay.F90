module replay
!=====================================================================
!
! Purpose: Implement replay: force U,V,T,Q towards reanalysis 
!
! Author: Sarah Weidman
!
! Description:
!         
!=====================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,get_curr_date,get_step_size
  use phys_grid   ,   only:scatter_field_to_chunk, get_ncols_p, get_lat_all_p,get_lon_all_p, &
        get_rlat_all_p,get_lat_p,get_lon_p
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc
  use cam_logfile ,   only:iulog
  use physics_types,    only: physics_state, physics_tend, physics_ptend, physics_update, &
        physics_ptend_init
  use ppgrid,           only: begchunk, endchunk, pcols, pver, pverp, psubcols
  use phys_control, only: phys_getopts
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: replay_readnl
  public:: replay_register
  public:: replay_correction
  private:: read_netcdf_replay
  private:: interpret_filename_replay
  public:: Replay_Model ! for if statements

  ! Replay parameters
  logical          :: Replay_Model       =.false.
  character(len=cl):: Replay_Path
  character(len=cs):: Replay_File_Template
  integer          :: Replay_Beg_Year
  real(r8)         :: Replay_coef

  ! replay observation arrays
  real(r8),allocatable::Ufield3d (:,:,:) !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Vfield3d (:,:,:) !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Tfield3d (:,:,:) !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Qfield3d (:,:,:) !(pcols,pver,begchunk:endchunk)

contains
  !================================================================
  subroutine replay_readnl(nlfile)
   ! 
   ! REPLAY_READNL: Initialize default values controlling the Replay 
   !                 process. Then read namelist values to override 
   !                 them.
   !===============================================================
   use namelist_utils,only:find_group_name
   use units         ,only:getunit,freeunit
   !
   ! Arguments
   !-------------
   character(len=*),intent(in)::nlfile
   !
   ! Local Values
   !---------------
   integer ierr,unitn

   namelist /replay_nl/ Replay_Model,Replay_Path,                       &
                         Replay_File_Template, Replay_Beg_Year,         &
                         Replay_coef 
                         


   ! Set Default Namelist values
   !-----------------------------
   Replay_Model         = .false.
   Replay_Path          = '/n/holylfs04/LABS/kuang_lab/Lab/sweidman/MERRA2_OG/MERRA2_f19/'
   Replay_File_Template = 'MERRA2_%y%m%d_%h.nc'
   Replay_Beg_Year      = 1980
   Replay_coef          = 1._r8

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'replay_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,replay_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('replay_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Replay_Path         ,len(Replay_Path)         ,mpichar,0,mpicom)
   call mpibcast(Replay_File_Template,len(Replay_File_Template),mpichar,0,mpicom)
   call mpibcast(Replay_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(Replay_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Replay_coef         , 1, mpir8 , 0, mpicom)
#endif

if(masterproc) then
  write(iulog,*) ' '
  write(iulog,*) '---------------------------------------------------------'
  write(iulog,*) '  MODEL REPLAY INITIALIZED WITH THE FOLLOWING SETTINGS: '
  write(iulog,*) '---------------------------------------------------------'
  write(iulog,*) 'REPLAY: Replay_Model=',Replay_Model
  write(iulog,*) 'REPLAY: Replay_Path=',Replay_Path
  write(iulog,*) 'REPLAY: Replay_File_Template =',Replay_File_Template
  write(iulog,*) 'REPLAY: Replay_Beg_Year =',Replay_Beg_Year
  write(iulog,*) 'REPLAY: Replay_coef  =',Replay_coef
endif

   ! End Routine
   !------------
   return
  end subroutine ! replay_readnl
  !================================================================

  subroutine replay_register

    use physics_buffer,     only: pbuf_add_field, dtype_r8
    use rad_constituents,   only: rad_cnst_get_info
    use constituents,     only : pcnst
#ifdef CRM
    use crmdims,            only: crm_nx, crm_ny, crm_nz 
#endif
    !---------------------------Local variables-----------------------------
    !
    ! ids for buffer vars
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
    integer ::  crm_u_oldid        = 0 
    integer ::  crm_v_oldid        = 0 
    integer ::  crm_w_oldid        = 0 
    integer ::  crm_t_oldid        = 0
    integer ::  crm_qrad_oldid        = 0 
    integer ::  crm_qt_oldid        = 0 
    integer ::  crm_qp_oldid        = 0 
    integer ::  crm_qn_oldid        = 0  

    logical :: use_spcam ! added for spcam
    integer  :: nmodes
    !-----------------------------------------------------------------------

    call rad_cnst_get_info(0, nmodes=nmodes)
    call phys_getopts( use_spcam_out = use_SPCAM)

    ! add pbuf vars from cam.r. to buffer - sweidman
    call pbuf_add_field('TEOUT_OLD', 'global', dtype_r8, (/pcols/), TEOUT_oldid) !(lat, lon) ; 
    call pbuf_add_field('DTCORE_OLD', 'global', dtype_r8, (/pcols,pver/), DTCORE_oldid  ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CLDO_OLD', 'global', dtype_r8, (/pcols,pver/),  CLDO_oldid         ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('PRER_EVAP_OLD', 'global', dtype_r8, (/pcols,pver/),  PRER_EVAP_oldid    ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CC_T_OLD', 'global', dtype_r8, (/pcols,pver/),  CC_T_oldid         ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CC_qv_OLD', 'global', dtype_r8, (/pcols,pver/),  CC_qv_oldid        ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CC_ql_OLD', 'global', dtype_r8, (/pcols,pver/),  CC_ql_oldid        ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CC_qi_OLD', 'global', dtype_r8, (/pcols,pver/),  CC_qi_oldid        ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CC_nl_OLD', 'global', dtype_r8, (/pcols,pver/),  CC_nl_oldid        ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CC_ni_OLD', 'global', dtype_r8, (/pcols,pver/),  CC_ni_oldid        ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CC_qlst_OLD', 'global', dtype_r8, (/pcols,pver/),  CC_qlst_oldid      ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('am_evp_st_OLD', 'global', dtype_r8, (/pcols,pver/),  am_evp_st_oldid    ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('evprain_st_OLD', 'global', dtype_r8, (/pcols,pver/),  evprain_st_oldid   ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('evpsnow_st_OLD', 'global', dtype_r8, (/pcols,pver/),  evpsnow_st_oldid   ) !(pbuf_00032, lat, lon)
    call pbuf_add_field('ACPRECL_OLD', 'global', dtype_r8, (/pcols/), ACPRECL_oldid      ) !(lat, lon) ;
    call pbuf_add_field('ACGCME_OLD', 'global', dtype_r8, (/pcols/), ACGCME_oldid       ) !(lat, lon) ;
    call pbuf_add_field('ACNUM_OLD', 'global', dtype_r8, (/pcols/), ACNUM_oldid        ) !(lat, lon) ;
    call pbuf_add_field('RELVAR_OLD', 'global', dtype_r8, (/pcols,pver/),  RELVAR_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('ACCRE_ENHAN_OLD', 'global', dtype_r8, (/pcols,pver/),  ACCRE_ENHAN_oldid  ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('pblh_OLD', 'global', dtype_r8, (/pcols/), pblh_oldid         ) !(lat, lon) ;
    call pbuf_add_field('tke_OLD', 'global', dtype_r8, (/pcols,pverp/),  tke_oldid          ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('kvh_OLD', 'global', dtype_r8, (/pcols,pverp/),  kvh_oldid          ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('tpert_OLD', 'global', dtype_r8, (/pcols/), tpert_oldid        ) !(lat, lon) ;
    call pbuf_add_field('AST_OLD', 'global', dtype_r8, (/pcols,pver/),  AST_oldid          ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('AIST_OLD', 'global', dtype_r8, (/pcols,pver/),  AIST_oldid         ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('ALST_OLD', 'global', dtype_r8, (/pcols,pver/),  ALST_oldid         ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('QIST_OLD', 'global', dtype_r8, (/pcols,pver/),  QIST_oldid         ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('QLST_OLD', 'global', dtype_r8, (/pcols,pver/),  QLST_oldid         ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CONCLD_OLD', 'global', dtype_r8, (/pcols,pver/),  CONCLD_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('CLD_OLD', 'global', dtype_r8, (/pcols,pver/),  CLD_oldid          ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('RAD_CLUBB_OLD', 'global', dtype_r8, (/pcols,pver/),  RAD_CLUBB_oldid    ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('WP2_nadv_OLD', 'global', dtype_r8, (/pcols,pverp/),  WP2_nadv_oldid     ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('WP3_nadv_OLD', 'global', dtype_r8, (/pcols,pverp/),  WP3_nadv_oldid     ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('WPTHLP_nadv_OLD', 'global', dtype_r8, (/pcols,pverp/),  WPTHLP_nadv_oldid  ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('WPRTP_nadv_OLD', 'global', dtype_r8, (/pcols,pverp/),  WPRTP_nadv_oldid   ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('RTPTHLP_nadv_OLD', 'global', dtype_r8, (/pcols,pverp/),  RTPTHLP_nadv_oldid ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('RTP2_nadv_OLD', 'global', dtype_r8, (/pcols,pverp/),  RTP2_nadv_oldid    ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('THLP2_nadv_OLD', 'global', dtype_r8, (/pcols,pverp/),  THLP2_nadv_oldid   ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('UP2_nadv_OLD', 'global', dtype_r8, (/pcols,pverp/),  UP2_nadv_oldid     ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('VP2_nadv_OLD', 'global', dtype_r8, (/pcols,pverp/),  VP2_nadv_oldid     ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('UPWP_OLD', 'global', dtype_r8, (/pcols,pverp/),  UPWP_oldid         ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('VPWP_OLD', 'global', dtype_r8, (/pcols,pverp/),  VPWP_oldid         ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('THLM_OLD', 'global', dtype_r8, (/pcols,pverp/),  THLM_oldid         ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('RTM_OLD', 'global', dtype_r8, (/pcols,pverp/),  RTM_oldid          ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('UM_OLD', 'global', dtype_r8, (/pcols,pverp/),  UM_oldid           ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('VM_OLD', 'global', dtype_r8, (/pcols,pverp/),  VM_oldid           ) !(pbuf_00033, lat, lon) ;
    
    call pbuf_add_field('DGNUM_OLD', 'global', dtype_r8, (/pcols,pverp,nmodes/),  DGNUM_oldid        ) !(pbuf_00128, lat, lon) ;
    call pbuf_add_field('DGNUMWET_OLD', 'global', dtype_r8, (/pcols,pverp,nmodes/),  DGNUMWET_oldid     ) !(pbuf_00128, lat, lon) ;
    call pbuf_add_field('num_c1_OLD', 'global', dtype_r8, (/pcols,pver/),  num_c1_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('so4_c1_OLD', 'global', dtype_r8, (/pcols,pver/),  so4_c1_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('pom_c1_OLD', 'global', dtype_r8, (/pcols,pver/),  pom_c1_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('soa_c1_OLD', 'global', dtype_r8, (/pcols,pver/),  soa_c1_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('bc_c1_OLD', 'global', dtype_r8, (/pcols,pver/),  bc_c1_oldid        ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('dst_c1_OLD', 'global', dtype_r8, (/pcols,pver/),  dst_c1_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('ncl_c1_OLD', 'global', dtype_r8, (/pcols,pver/),  ncl_c1_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('num_c2_OLD', 'global', dtype_r8, (/pcols,pver/),  num_c2_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('so4_c2_OLD', 'global', dtype_r8, (/pcols,pver/),  so4_c2_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('soa_c2_OLD', 'global', dtype_r8, (/pcols,pver/),  soa_c2_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('ncl_c2_OLD', 'global', dtype_r8, (/pcols,pver/),  ncl_c2_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('dst_c2_OLD', 'global', dtype_r8, (/pcols,pver/),  dst_c2_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('num_c3_OLD', 'global', dtype_r8, (/pcols,pver/),  num_c3_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('dst_c3_OLD', 'global', dtype_r8, (/pcols,pver/),  dst_c3_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('ncl_c3_OLD', 'global', dtype_r8, (/pcols,pver/),  ncl_c3_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('so4_c3_OLD', 'global', dtype_r8, (/pcols,pver/),  so4_c3_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('num_c4_OLD', 'global', dtype_r8, (/pcols,pver/),  num_c4_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('pom_c4_OLD', 'global', dtype_r8, (/pcols,pver/),  pom_c4_oldid       ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('bc_c4_OLD', 'global', dtype_r8, (/pcols,pver/),  bc_c4_oldid        ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('DP_FLXPRC_OLD', 'global', dtype_r8, (/pcols,pverp/),  DP_FLXPRC_oldid    ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('DP_FLXSNW_OLD', 'global', dtype_r8, (/pcols,pverp/),  DP_FLXSNW_oldid    ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('DP_CLDLIQ_OLD', 'global', dtype_r8, (/pcols,pver/),  DP_CLDLIQ_oldid    ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('DP_CLDICE_OLD', 'global', dtype_r8, (/pcols,pver/),  DP_CLDICE_oldid    ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('cush_OLD', 'global', dtype_r8, (/pcols/), cush_oldid         ) !(lat, lon) ;
    call pbuf_add_field('QRS_OLD', 'global', dtype_r8, (/pcols,pver/),  QRS_oldid          ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('QRL_OLD', 'global', dtype_r8, (/pcols,pver/),  QRL_oldid          ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('ICIWP_OLD', 'global', dtype_r8, (/pcols,pver/),  ICIWP_oldid        ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('ICLWP_OLD', 'global', dtype_r8, (/pcols,pver/),  ICLWP_oldid        ) !(pbuf_00032, lat, lon) ;
    call pbuf_add_field('kvm_OLD', 'global', dtype_r8, (/pcols,pverp/),  kvm_oldid          ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('turbtype_OLD', 'global', dtype_r8, (/pcols,pverp/),  turbtype_oldid     ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('smaw_OLD', 'global', dtype_r8, (/pcols,pverp/),  smaw_oldid         ) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('tauresx_OLD', 'global', dtype_r8, (/pcols/), tauresx_oldid      ) !(lat, lon) ;
    call pbuf_add_field('tauresy_OLD', 'global', dtype_r8, (/pcols/), tauresy_oldid      ) !(lat, lon) ;
    call pbuf_add_field('qpert_OLD', 'global', dtype_r8, (/pcols,pcnst/), qpert_oldid) !(pbuf_00033, lat, lon) ;
    call pbuf_add_field('T_TTEND_OLD', 'global', dtype_r8, (/pcols,pver/),  T_TTEND_oldid      ) !(pbuf_00032, lat, lon) ;

    if (use_SPCAM) then
#ifdef CRM
  if(masterproc) write(iulog,*) 'registering CRM buffer vars'
      call pbuf_add_field('CRM_U_OLD', 'global', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/), crm_u_oldid)
      call pbuf_add_field('CRM_V_OLD', 'global', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/), crm_v_oldid)
      call pbuf_add_field('CRM_W_OLD', 'global', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/), crm_w_oldid)
      call pbuf_add_field('CRM_T_OLD', 'global', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/), crm_t_oldid)
      call pbuf_add_field('CRM_QRAD_OLD', 'global',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/), crm_qrad_oldid)
      call pbuf_add_field('CRM_QT_OLD', 'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), crm_qt_oldid)
      call pbuf_add_field('CRM_QP_OLD', 'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), crm_qp_oldid)
      call pbuf_add_field('CRM_QN_OLD', 'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), crm_qn_oldid)
#endif
    end if

  end subroutine replay_register 

  !================================================================

  subroutine read_netcdf_replay(analysis_file, Replay_nlon, Replay_nlat)
    use netcdf
    use ppgrid,    only: pver,pcols,begchunk,endchunk
    use phys_grid,   only: scatter_field_to_chunk
    use dyn_grid      ,only: get_horiz_grid_dim_d
    use error_messages,only: alloc_err
    use cam_abortutils, only:endrun
    !implicit none
 
    ! Arguments
    !-------------
    character(len=*),intent(in):: analysis_file
 
    ! Local values
    !-------------
    integer lev
    integer nlon,nlat,plev,istat
    integer Replay_nlon,Replay_nlat,Replay_nlev
    integer ncid,varid
    integer ilat,ilon,ilev
    real(r8) Xanal(Replay_nlon,Replay_nlat,pver) ! TODO: need to define Replay_lat et al
    real(r8) Lat_anal(Replay_nlat)
    real(r8) Lon_anal(Replay_nlon)
    real(r8) Xtrans(Replay_nlon,pver,Replay_nlat)
    !integer  hdim1_d,hdim2_d
 
    Replay_nlev=pver
 
    if(masterproc) then
    ! Open the given file
      !-----------------------
    istat=nf90_open(trim(analysis_file),NF90_NOWRITE,ncid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*)'NF90_OPEN: failed for file ',trim(analysis_file)
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    ! Read in Dimensions
      !--------------------
    istat=nf90_inq_dimid(ncid,'lon',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_inquire_dimension(ncid,varid,len=nlon)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    istat=nf90_inq_dimid(ncid,'lat',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_inquire_dimension(ncid,varid,len=nlat)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    istat=nf90_inq_dimid(ncid,'lev',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_inquire_dimension(ncid,varid,len=plev)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    istat=nf90_inq_varid(ncid,'lon',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Lon_anal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    istat=nf90_inq_varid(ncid,'lat',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Lat_anal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    if((Replay_nlon.ne.nlon).or.(Replay_nlat.ne.nlat).or.(plev.ne.pver)) then
     write(iulog,*) 'ERROR: replay_update_analyses_fv: nlon=',nlon,' Replay_nlon=',Replay_nlon
     write(iulog,*) 'ERROR: replay_update_analyses_fv: nlat=',nlat,' Replay_nlat=',Replay_nlat
     write(iulog,*) 'ERROR: replay_update_analyses_fv: plev=',plev,' pver=',pver
     call endrun('replay_update_analyses_fv: analyses dimension mismatch')
    endif
 
    ! Read in, transpose lat/lev indices, 
      ! and scatter data arrays
      !----------------------------------
    istat=nf90_inq_varid(ncid,'U',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Xanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    do ilat=1,nlat
    do ilev=1,plev
    do ilon=1,nlon
      Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    end do
    end do
    end do
  endif ! (masterproc) then
  call scatter_field_to_chunk(1,Replay_nlev,1,Replay_nlon,Xtrans,   &
                              Ufield3d(1,1,begchunk))
 
  if(masterproc) then
    istat=nf90_inq_varid(ncid,'V',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Xanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    do ilat=1,nlat
    do ilev=1,plev
    do ilon=1,nlon
      Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    end do
    end do
    end do
  endif ! (masterproc) then
  call scatter_field_to_chunk(1,Replay_nlev,1,Replay_nlon,Xtrans,   &
                              Vfield3d(1,1,begchunk))
 
  if(masterproc) then
    istat=nf90_inq_varid(ncid,'T',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Xanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    do ilat=1,nlat
    do ilev=1,plev
    do ilon=1,nlon
      Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    end do
    end do
    end do
  endif ! (masterproc) then
  call scatter_field_to_chunk(1,Replay_nlev,1,Replay_nlon,Xtrans,   &
                              Tfield3d(1,1,begchunk))
 
  if(masterproc) then
    istat=nf90_inq_varid(ncid,'Q',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Xanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    do ilat=1,nlat
    do ilev=1,plev
    do ilon=1,nlon
      Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    end do
    end do
    end do
 
  ! Close analysis file
  istat=nf90_close(ncid)
  if(istat.ne.NF90_NOERR) then
    write(iulog,*) nf90_strerror(istat)
    call endrun ('UPDATE_ANALYSES_EUL')
  endif
 
  endif ! (masterproc) then
  call scatter_field_to_chunk(1,Replay_nlev,1,Replay_nlon,Xtrans,   &
                              Qfield3d(1,1,begchunk))
 
    ! End Routine
    !------------
    return
 
   end subroutine read_netcdf_replay

!================================================================

   character(len=cl) function interpret_filename_replay( filename_spec, case, &
   yr_spec, mon_spec, day_spec, hr_spec, sec_spec )

! Create a filename from a filename specifier. The 
! filename specifyer includes codes for setting things such as the
! year, month, day, seconds in day, caseid, and tape number. 
!
! Interpret filename specifyer string with: 
!
!      %c for case, 
!      %y for year
!      %m for month
!      %d for day
!      %h for modstep
!      %% for the "%" character
!
! If the filename specifyer has spaces " ", they will be trimmed out
! of the resulting filename.

   ! arguments
   character(len=*), intent(in)           :: filename_spec   ! Filename specifier to use
   character(len=*), intent(in), optional :: case            ! Optional casename
   integer         , intent(in), optional :: yr_spec         ! Simulation year
   integer         , intent(in), optional :: mon_spec        ! Simulation month
   integer         , intent(in), optional :: day_spec        ! Simulation day
   integer         , intent(in), optional :: hr_spec         ! Modstep from replay correction
   integer         , intent(in), optional :: sec_spec        ! Simulation seconds of day

   ! Local variables
   integer :: year  ! Simulation year
   integer :: month ! Simulation month
   integer :: day   ! Simulation day
   integer :: ncsec   ! Seconds into current simulation day
   integer :: modstep ! modstep into current simulation day
   character(len=cl) :: string    ! Temporary character string 
   character(len=cl) :: format    ! Format character string 
   integer :: i, n  ! Loop variables
   logical :: done
   !-----------------------------------------------------------------------------


   if ( len_trim(filename_spec) == 0 )then
      call endrun ('INTERPRET_FILENAME_SPEC: filename specifier is empty')
   end if
   if ( index(trim(filename_spec)," ") /= 0 )then
      call endrun ('INTERPRET_FILENAME_SPEC: filename specifier can not contain a space:'//trim(filename_spec))
   end if
   !
   ! Determine year, month, day and sec to put in filename
   !
   if (present(yr_spec) .and. present(mon_spec) .and. present(day_spec) .and. present(hr_spec) .and. present(sec_spec)) then
      year  = yr_spec
      month = mon_spec
      day   = day_spec
      modstep = hr_spec
      ncsec = sec_spec
   end if
   !
   ! Go through each character in the filename specifyer and interpret if special string
   !
   i = 1
   interpret_filename_replay = ''
   do while ( i <= len_trim(filename_spec) )
      !
      ! If following is an expansion string
      !
      if ( filename_spec(i:i) == "%" )then
         i = i + 1
         select case( filename_spec(i:i) )
         case( 'y' )   ! year
            if ( year > 99999   ) then
               format = '(i6.6)'
            else if ( year > 9999    ) then
               format = '(i5.5)'
            else
               format = '(i4.4)'
            end if
            write(string,format) year
         case( 'm' )   ! month
            write(string,'(i2.2)') month
         case( 'd' )   ! day
            write(string,'(i2.2)') day
         case( 'h' )   ! 3-hour period
            write(string,'(i1.1)') modstep
         case( 's' )   ! second
            write(string,'(i5.5)') ncsec
         case( '%' )   ! percent character
            string = "%"
         case default
            call endrun ('INTERPRET_FILENAME_SPEC: Invalid expansion character: '//filename_spec(i:i))
         end select
         !
         ! Otherwise take normal text up to the next "%" character
         !
      else
         n = index( filename_spec(i:), "%" )
         if ( n == 0 ) n = len_trim( filename_spec(i:) ) + 1
         if ( n == 0 ) exit 
         string = filename_spec(i:n+i-2)
         i = n + i - 2
      end if
      if ( len_trim(interpret_filename_replay) == 0 )then
        interpret_filename_replay = trim(string)
      else
         if ( (len_trim(interpret_filename_replay)+len_trim(string)) >= cl )then
            call endrun ('INTERPRET_FILENAME_SPEC: Resultant filename too long')
         end if
         interpret_filename_replay = trim(interpret_filename_replay) // trim(string)
      end if
      i = i + 1

   end do
   if ( len_trim(interpret_filename_replay) == 0 )then
      call endrun ('INTERPRET_FILENAME_SPEC: Resulting filename is empty')
   end if

end function interpret_filename_replay

   subroutine replay_correction (state,tend,ztodt)

    !----------------------------------------------------------------------- 
    ! Purpose: 
    !-----------------------------------------------------------------------
        use physics_buffer, only : pbuf_get_index, pbuf_get_field,physics_buffer_desc, pbuf_set_field, pbuf_add_field, dtype_r8, pbuf_get_index, pbuf_old_tim_idx
        use phys_grid,    only:  gather_chunk_to_field, scatter_field_to_chunk
        use dyn_grid,     only: get_horiz_grid_dim_d
        use time_manager, only: get_nstep, get_curr_date, get_start_date
        use geopotential, only: geopotential_dse
        use physconst,    only: zvir, gravit, cpairv, rair,cpair
        use cam_pio_utils,    only: cam_pio_openfile, cam_pio_get_decomp ! sweid - replace with cam_pio_get_decomp ?
        use cam_grid_support,   only: cam_grid_get_decomp, cam_grid_id, cam_grid_dimensions ! trying this? - sweid
        use pio,          only: pio_write_darray, pio_read_darray, file_desc_t, var_desc_t, io_desc_t, pio_offset, pio_setframe, pio_double, pio_write, pio_nowrite, pio_inq_varid, pio_def_var, pio_closefile
        use pio_types, only : file_desc_t, var_desc_t, io_desc_t
        use ioFileMod,     only: getfil
        use cam_history,    only: addfld, outfld
        use shr_mem_mod,       only: shr_mem_init, shr_mem_getusage
        use pmgrid,          only: plon, plat
        use constituents,     only: cnst_get_ind, pcnst
        use check_energy,    only: check_energy_chng
        use error_messages, only: alloc_err 
    
       real(r8) :: qtarget_c(pcols,begchunk:endchunk,pver)
       real(r8) :: utarget_c(pcols,begchunk:endchunk,pver)
       real(r8) :: vtarget_c(pcols,begchunk:endchunk,pver)
       real(r8) :: starget_c(pcols,begchunk:endchunk,pver)
       real(r8) :: ttarget_c(pcols,begchunk:endchunk,pver)
    
       real(r8) :: qforcing(pcols,begchunk:endchunk,pver)
       real(r8) :: uforcing(pcols,begchunk:endchunk,pver)
       real(r8) :: vforcing(pcols,begchunk:endchunk,pver)
       real(r8) :: sforcing(pcols,begchunk:endchunk,pver)
    
       integer, save :: nstep_count
    
    ! Arguments
        type(physics_state), intent(inout) :: state(begchunk:endchunk)
        type(physics_tend), intent(inout) :: tend(begchunk:endchunk)
        real(r8) , intent(in) :: ztodt
    ! Local workspace
        type(physics_ptend)   :: ptend                  ! indivdualparameterization tendencies
        real(r8) ::     arrq(pcols,begchunk:endchunk,pver),arrt(pcols,begchunk:endchunk,pver),arru(pcols,begchunk:endchunk,pver),arrv(pcols,begchunk:endchunk,pver)       ! Input array,chunked
        real(r8) :: zmean                        ! temporary zonal mean value
        integer :: i, j, qconst, ifld,n ,ilat,ilon,istat                 ! longitude, latitude,field, and global column indices
        integer :: hdim1, hdim2, c, ncols, k, istep, modstep
        integer :: hdim1_d, hdim2_d, Replay_nlon, Replay_nlat
        real(r8) ::rlat(pcols),damping_coef,wrk,forcingtime,dampingtime!,tmprand(128,64)
        real :: tmp_zero
        real(r8) :: zero(pcols)                    ! array of zeros
        integer :: ilat_all(pcols)
        integer :: ndays, day, mon, yr, ncsec
        integer :: yr_start, mon_start, day_start, tod_start
        integer :: modstep6hr, modstep3hr
        logical :: fileexists
        logical,save :: corrector_step, started   
    !----- 
        integer :: ierr,csize,indw                          !!Added  
        integer jerr
        character(len=256) :: ncdata_loc,filen,ncwrite_loc,outn
        integer :: dims2d(2)
        integer :: resul,lchnk
        real(r8), pointer :: londeg(:,:)
        type(file_desc_t) :: File
    
       !---------------------------flamraoui 2------------------------------
      ! variable for reading netcdf 
        character(len=256)        :: fileName
    
        integer, dimension(11) :: monarray
    
        real(r8), allocatable :: tmpfield(:)
        type(io_desc_t), pointer :: iodesc
        integer                  :: dims(3), gdims(3), nhdims ! added - sweid
        integer                  :: physgrid ! added - sweid
        type(var_desc_t) :: vardesc
        real(r8) :: msize,mrss
 
        logical  :: lq(pcnst)
    
      !---------------------------flamraoui------------------------------
      
       istep=get_nstep()
       nstep_count=istep
    
    !On first time step, make sure we're starting a "clean" run
    
    if (nstep_count==0 .AND. started .ne. .TRUE.) then
        started=.TRUE.
        corrector_step=.FALSE.
        
        #if ( defined SPMD )
            do c = begchunk, endchunk
                call get_rlat_all_p(c,pcols,rlat)
                call get_lat_all_p(c,pcols,ilat_all)
                rlat=rlat*180._r8/3.14159
                ncols = get_ncols_p(c)
                do i = 1, ncols
                    do k=1,pver
                        state(c)%qforce(i,k) = 0.0
                        state(c)%uforce(i,k) = 0.0
                        state(c)%vforce(i,k) = 0.0
                        state(c)%sforce(i,k) = 0.0
                    end do
                end do
            end do
        #endif
    endif
    
    !then figure out what day it is
    
      !---------------------------flamraoui------------------------------
         ndays = 1+ istep/48     !  number of days based on istep
    !------------------------------------------------     
    
    call get_curr_date(yr, mon, day, ncsec)
    call get_start_date(yr_start, mon_start, day_start, tod_start)
    
    yr = yr + (Replay_Beg_Year - yr_start)
    
    if (masterproc) then
       write(iulog,*)  'iyear = ', yr
       write(iulog,*)  'imonth = ', mon
       write(iulog,*)  'iday = ', day
       write(iulog,*) "step count from the beginning", nstep_count
    endif
 
    !define constants
        forcingtime=1._r8*21600._r8       ! six hours in seconds
    
 
        call get_horiz_grid_dim_d(hdim1, hdim2)
    
    
    modstep=int((mod(istep+6, 48)) / 6)
    modstep6hr=  (mod(istep, 12)) 
    modstep3hr=  (mod(istep, 6)) 
    
    
    
      !------------------------------------------
       
       
    !goals at end of "corrector" run
    !a) reset forcing to zero
    !b) set flag to false 
       
       if (modstep6hr == 11 ) then
    !       if (corrector_step) then
               #if ( defined SPMD )
               do c = begchunk, endchunk
                   call get_rlat_all_p(c,pcols,rlat)
                   call get_lat_all_p(c,pcols,ilat_all)
                   rlat=rlat*180._r8/3.14159
                   ncols = get_ncols_p(c)
                   do i = 1, ncols
                       do k=1,pver
                            state(c)%qforce(i,k) = 0.0
                            state(c)%uforce(i,k) = 0.0
                            state(c)%vforce(i,k) = 0.0
                            state(c)%sforce(i,k) = 0.0 
                       end do
                   end do
               end do
               #endif
    !       a) wipe the forcing to zero
               corrector_step=.FALSE.
    
        endif
      !      set corrector_step is false

    
    if (masterproc) then
    print*, "modstep", modstep    ! timestep during the day 
    print*, "modstep6hr", modstep6hr  ! reset to zero every 6 hrs 
    print*, "modstep3hr", modstep3hr  ! reset to zero every 3 hrs
    print*, "corrector_step", corrector_step
    print*, state(begchunk)%sforce(1,1)
    endif
    
    !
    !if in "corrector" step: divide difference by 6h to get tendency and apply it to physics
    !
    if (corrector_step) then
       if(masterproc) then
          print *, 'applying corrector ptend'
          !write(iulog,*) 'replay: Reading analyses:',trim(filename)
       endif
    #if ( defined SPMD )
    do c = begchunk, endchunk
           ! reallocate ptend
           call cnst_get_ind('Q',indw)
           lq(:)   =.false.
           lq(indw)=.true.
           call physics_ptend_init(ptend, state(c)%psetcols, "none", ls=.true.,  lu=.true., lv=.true., lq=lq) 
 
           ncols = get_ncols_p(c)
    !  print *, "=========doloop ====>  ncols= ", ncols
    do i = 1, ncols
    do k=1,pver   
    
        ptend%q(i,k,indw) = ptend%q(i,k,indw) + state(c)%qforce(i,k)/forcingtime*Replay_coef 
        ptend%u(i,k) = ptend%u(i,k) + state(c)%uforce(i,k)/forcingtime*Replay_coef 
        ptend%v(i,k) = ptend%v(i,k) + state(c)%vforce(i,k)/forcingtime*Replay_coef 
        ptend%s(i,k) = ptend%s(i,k) + state(c)%sforce(i,k)/forcingtime*Replay_coef 
    
    !
    end do
    end do
    !print *, "ptendusize",size(ptend%u)
    !apply tendencies to model

        call physics_update (state(c), ptend, ztodt, tend(c)) ! this calls ptend deallocate
        call check_energy_chng(state(c), tend(c), "replay", istep, ztodt, zero, zero, zero, zero)
    end do
    #endif
    endif
    
    
    !goals at end of "clean" run
    ! a) reading merra target
    ! b) define forcing to be difference with the state divided by 6hrs 
    ! c) write (save) it to disk 
    ! d) set corrector_step to true 
    ! e) reset clock and restart states (not here)
    
        if  (modstep6hr==5 .AND. .NOT. corrector_step ) then

          fileexists=.FALSE.
    
    do while (.NOT. fileexists )

      if (masterproc) write(iulog,*) "modstep, ncsec + 1800", modstep, ncsec+1800
    
      filename=interpret_filename_replay(Replay_File_Template      , &
          yr_spec=yr , &
          mon_spec=mon, &
          day_spec=day  , &
          hr_spec=modstep, &
          sec_spec=ncsec+1800    )

      if(masterproc) then
        write(iulog,*) trim(Replay_Path)//trim(filename)
      endif
      
      INQUIRE(FILE=trim(Replay_Path)//trim(filename), EXIST=fileexists)
      
      if (.not. fileexists) print*, 'file missing', filename
      
      day=day-1
      if (day==0) then
      day=31
      mon=mon-1
      if (mon==0) then
      mon=12
      yr=yr-1
      end if
      end if
      
    end do ! checking if file exists
   
    !-------------------------------------------------------------------------------------
    if (masterproc) write(iulog,*) "Reanalysis filename = ", trim(Replay_Path)//trim(filename)   ! print filename used 
    !-----------------------------finish calling filename--------------------------------------  
    
         
          call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
          Replay_nlon=hdim1_d
          Replay_nlat=hdim2_d
          
          allocate(Ufield3d(pcols,pver,begchunk:endchunk),stat=istat)
          call alloc_err(istat,'replay_init','Ufield3d',pcols*pver*((endchunk-begchunk)+1))
          allocate(Vfield3d(pcols,pver,begchunk:endchunk),stat=istat)
          call alloc_err(istat,'replay_init','Vfield3d',pcols*pver*((endchunk-begchunk)+1))
          allocate(Tfield3d(pcols,pver,begchunk:endchunk),stat=istat)
          call alloc_err(istat,'replay_init','Tfield3d',pcols*pver*((endchunk-begchunk)+1))
          allocate(Qfield3d(pcols,pver,begchunk:endchunk),stat=istat)
          call alloc_err(istat,'replay_init','Qfield3d',pcols*pver*((endchunk-begchunk)+1))
 
          Ufield3d(:pcols,:pver,begchunk:endchunk)=0._r8
          Vfield3d(:pcols,:pver,begchunk:endchunk)=0._r8
          Tfield3d(:pcols,:pver,begchunk:endchunk)=0._r8
          Qfield3d(:pcols,:pver,begchunk:endchunk)=0._r8
 
          call read_netcdf_replay(trim(Replay_Path)//trim(filename), Replay_nlon, Replay_nlat)
    
       !call pio_closefile(File)
       if(masterproc) then 
       write(iulog,*) "done read in reanalysis"
       write(iulog,*) "state(c)%sforce(1,1): ", state(begchunk)%sforce(1,1)
       write(iulog,*) "begchunk: ", begchunk
       write(iulog,*) "anal_field T(1,1,1): ", Tfield3d(1,1,begchunk)
       endif
    
            do c = begchunk, endchunk
                ncols = get_ncols_p(c)
                do i = 1, ncols
                    do k=1,pver
                        state(c)%qforce(i,k)=((Qfield3d(i,k,c)-state(c)%q(i,k,1)))
                        state(c)%uforce(i,k)=((Ufield3d(i,k,c)-state(c)%u(i,k)))
                        state(c)%vforce(i,k)=((Vfield3d(i,k,c)-state(c)%v(i,k)))
                        state(c)%sforce(i,k)=((Tfield3d(i,k,c)-state(c)%t(i,k)))*cpair
                    end do
                end do
            end do
 
            if(masterproc) then
            write(iulog,*) "done update state"
            write(iulog,*) "state(c)%sforce(1,1): ", state(begchunk)%sforce(1,1)
            write(iulog,*) "state(c)%sforce(2,2): ", state(begchunk)%sforce(2,2)
            endif
    
            !deallocate(tmpfield)
            deallocate(Tfield3d)
            deallocate(Ufield3d)
            deallocate(Vfield3d)
            deallocate(Qfield3d)
            !deallocate(Zfield3d)
    !writing happens in cam_diagnostics
    
            corrector_step=.TRUE.
    
    end if
    
    end subroutine replay_correction


  !================================================================

end module replay
