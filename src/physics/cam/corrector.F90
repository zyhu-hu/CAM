module corrector 
!=====================================================================
!
! Purpose: Implement nudging of the model state of U,V,T,Q, and/or PS
!          using constant state independent forcing. 
!
! Author: Sarah Weidman
!
! Description:
!    
!    This module assumes that the user has {U,V,S,Q} forcing values 
!    which have been preprocessed onto the current model grid and adjusted 
!    for differences in topography. It is also assumed that these resulting 
!    values and are stored in individual files which are indexed with respect 
!    to year, month, day, and second of the day. When the model is inbetween 
!    the given begining and ending times, a relaxation forcing is added to 
!    nudge the model with these forcing values. 
!
!    The nudging of the model with the forcing data is controlled by 
!    the 'corrector_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which nudging is applied, the strength of the nudging
!    tendencies, and its spatial distribution. 
!
!    FORCING:
!    --------
!    Nudging tendencies are applied as a relaxation force from the given
!    forcing in the input files, multiplied by a timescale / strength.
!    The nudging strength Alpha=[0.,1.] for each 
!    variable is specified by the 'Force_Xcoef' values. Where X={U,V,T,Q,PS}
!
!           F_nudge = Alpha*((Forcing(t_curr))/TimeScale
!
!
!    WINDOWING:
!    ----------
!    The region of applied nudging can be limited using Horizontal/Vertical 
!    window functions that are constructed using a parameterization of the 
!    Heaviside step function. 
!
!    The Heaviside window function is the product of separate horizonal and vertical 
!    windows that are controled via 12 parameters:
!
!        Force_Hwin_lat0:     Specify the horizontal center of the window in degrees. 
!        Force_Hwin_lon0:     The longitude must be in the range [0,360] and the 
!                             latitude should be [-90,+90].
!        Force_Hwin_latWidth: Specify the lat and lon widths of the window as positive 
!        Force_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999) 
!                             renders the window a constant in that direction.
!        Force_Hwin_latDelta: Controls the sharpness of the window transition with a 
!        Force_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step 
!                             function while a large value yeilds a smoother transition.
!        Force_Hwin_Invert  : A logical flag used to invert the horizontal window function 
!                             to get its compliment.(e.g. to nudge outside a given window).
!
!        Force_Vwin_Lindex:   In the vertical, the window is specified in terms of model 
!        Force_Vwin_Ldelta:   level indcies. The High and Low transition levels should 
!        Force_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also 
!        Force_Vwin_Hdelta:   specified in terms of model indices. For a window function 
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition 
!                             lengths should be set to 0.001 
!        Force_Vwin_Invert  : A logical flag used to invert the vertical window function 
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent 
!                 of the vertical (30 levels):
!                        Force_Hwin_lat0     = 0.         Force_Vwin_Lindex = 0.
!                        Force_Hwin_latWidth = 30.        Force_Vwin_Ldelta = 0.001
!                        Force_Hwin_latDelta = 5.0        Force_Vwin_Hindex = 31.
!                        Force_Hwin_lon0     = 180.       Force_Vwin_Hdelta = 0.001 
!                        Force_Hwin_lonWidth = 999.       Force_Vwin_Invert = .false.
!                        Force_Hwin_lonDelta = 1.0
!                        Force_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Force_Hwin_Invert = .true.
!
!    A user can preview the window resulting from a given set of namelist values before 
!    running the model. Lookat_NudgeWindow.ncl is a script avalable in the tools directory 
!    which will read in the values for a given namelist and display the resulting window.
!
!    The module is currently configured for only 1 window function. It can readily be 
!    extended for multiple windows if the need arises.
!
!
! Input/Output Values:
!
!    &nudging_nl
!      Force_Model         - LOGICAL toggle to activate nudging.
!                              TRUE  -> Nudging is on.
!                              FALSE -> Nudging is off.                            [DEFAULT]
!
!      Force_Path          - CHAR path to the analyses files.
!                              (e.g. '/glade/scratch/USER/inputdata/nudging/ERAI-Data/')
!
!      Force_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!                              (e.g. '%y/ERAI_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc')
!
!      Force_Times_Per_Day - INT Number of analyses files available per day.
!                              1 --> daily analyses.
!                              4 --> 6 hourly analyses.                            [DEFAULT]
!                              8 --> 3 hourly.
!
!
!      Force_Beg_Year      - INT nudging begining year.  [1979- ]
!      Force_Beg_Month     - INT nudging begining month. [1-12]
!      Force_Beg_Day       - INT nudging begining day.   [1-31]
!      Force_End_Year      - INT nudging ending year.    [1979-]
!      Force_End_Month     - INT nudging ending month.   [1-12]
!      Force_End_Day       - INT nudging ending day.     [1-31]
!
!
!      Force_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Force_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Force_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Force_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Force_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No Nudging of this variable)
!                                         1 == CONSTANT (Spatially Uniform Nudging)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      Force_Ucoef         - REAL fractional nudging coeffcient for U. 
!      Force_Vcoef         - REAL fractional nudging coeffcient for V. 
!      Force_Tcoef         - REAL fractional nudging coeffcient for T. 
!      Force_Qcoef         - REAL fractional nudging coeffcient for Q. 
!      Force_PScoef        - REAL fractional nudging coeffcient for PS. 
!
!                                 The strength of the nudging is specified as a fractional 
!                                 coeffcient between [0,1].
!           
!      Force_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Force_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Force_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Force_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Force_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Force_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Force_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!      Force_Vwin_Lindex   - REAL LO model index of transition
!      Force_Vwin_Hindex   - REAL HI model index of transition
!      Force_Vwin_Ldelta   - REAL LO transition length 
!      Force_Vwin_Hdelta   - REAL HI transition length 
!      Force_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!    /
!
!================
!
! TO DO:
! -----------
!    ** Implement Ps Nudging????
!          
!=====================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,get_curr_date,get_step_size
  use phys_grid   ,   only:scatter_field_to_chunk,gather_chunk_to_field
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc
  use cam_logfile ,   only:iulog
#ifdef SPMD
  use mpishorthand
#endif
  use torch_ftn
  use iso_fortran_env

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: Force_Model,Force_ON
  public:: corrector_readnl
  public:: corrector_init
  public:: corrector_timestep_init
  public:: nncorrector_timestep_init
  public:: corrector_timestep_tend
  public:: init_neural_net
  private::corrector_update_analyses_fv
  private::corrector_set_PSprofile
  private::corrector_set_profile

  ! corrector Parameters
  !--------------------
  logical          :: Force_Model       =.false.
  logical          :: Force_ON          =.false.
  logical          :: Force_Initialized =.false.
  character(len=cl):: Force_Path
  character(len=cs):: Force_File,Force_File_Template
  integer          :: Force_Times_Per_Day
  real(r8)         :: Force_Ucoef,Force_Vcoef
  integer          :: Force_Uprof,Force_Vprof
  real(r8)         :: Force_Qcoef,Force_Tcoef
  integer          :: Force_Qprof,Force_Tprof
  real(r8)         :: Force_PScoef
  integer          :: Force_PSprof
  integer          :: Force_Beg_Year ,Force_Beg_Month
  integer          :: Force_Beg_Day  ,Force_Beg_Sec
  integer          :: Force_End_Year ,Force_End_Month
  integer          :: Force_End_Day  ,Force_End_Sec
  integer          :: Force_Curr_Year,Force_Curr_Month
  integer          :: Force_Curr_Day ,Force_Curr_Sec
  integer          :: Force_Next_Year,Force_Next_Month
  integer          :: Force_Next_Day ,Force_Next_Sec
  integer          :: Force_Step
  real(r8)         :: Force_Hwin_lat0
  real(r8)         :: Force_Hwin_latWidth
  real(r8)         :: Force_Hwin_latDelta
  real(r8)         :: Force_Hwin_lon0
  real(r8)         :: Force_Hwin_lonWidth
  real(r8)         :: Force_Hwin_lonDelta
  logical          :: Force_Hwin_Invert = .false.
  real(r8)         :: Force_Hwin_lo
  real(r8)         :: Force_Hwin_hi
  real(r8)         :: Force_Vwin_Hindex
  real(r8)         :: Force_Vwin_Hdelta
  real(r8)         :: Force_Vwin_Lindex
  real(r8)         :: Force_Vwin_Ldelta
  logical          :: Force_Vwin_Invert =.false.
  real(r8)         :: Force_Vwin_lo
  real(r8)         :: Force_Vwin_hi
  real(r8)         :: Force_Hwin_latWidthH
  real(r8)         :: Force_Hwin_lonWidthH
  real(r8)         :: Force_Hwin_max
  real(r8)         :: Force_Hwin_min

  ! corrector State Arrays
  !-----------------------
  integer Force_nlon,Force_nlat,Force_ncol,Force_nlev
  real(r8),allocatable::Target_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_PS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Force_Utau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Force_Vtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Force_Stau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Force_Qtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Force_PStau (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Force_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Force_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Force_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Force_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Force_PSstep(:,:)    !(pcols,begchunk:endchunk)

  real(r8),allocatable::Model_state_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_state_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_state_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_state_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_state_QLIQ     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_state_QICE     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_state_OMEGA     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_state_PS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_SOLIN    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_LHFLX    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_SHFLX    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_SNOWHLND    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_PHIS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_TAUX    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_TAUY    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_TS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_ICEFRAC    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_LANDFRAC    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_lat    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_lon    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_tod    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_state_toy    (:,:)    !(pcols,begchunk:endchunk)

  ! corrector Observation Arrays
  !-----------------------------
  logical :: Force_File_Present

  ! NN related variables
  !---------------------
  integer :: nn_inputlength  = 197     ! length of NN input vector
  integer :: nn_outputlength = 104     ! length of NN output vector
  character(len=256)    :: torch_model='/n/holylfs04/LABS/kuang_lab/Lab/kuanglfs/zeyuanhu/climcorr/swin_test_dim1024_depth8_v2_2nodes_r4.pt'
  type(torch_module), allocatable :: torch_mod(:)

contains
  !================================================================
  subroutine corrector_readnl(nlfile)
   ! 
   ! corrector_READNL: Initialize default values controlling the corrector 
   !                 process. Then read namelist values to override 
   !                 them.
   !===============================================================
   use ppgrid        ,only: pver
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

   namelist /corrector_nl/ Force_Model,Force_Path,                       &
                         Force_File_Template,                          &
                         Force_Times_Per_Day,                          &
                         Force_Ucoef ,Force_Uprof,                     &
                         Force_Vcoef ,Force_Vprof,                     &
                         Force_Qcoef ,Force_Qprof,                     &
                         Force_Tcoef ,Force_Tprof,                     &
                         Force_PScoef,Force_PSprof,                    &
                         Force_Beg_Year,Force_Beg_Month,Force_Beg_Day, &
                         Force_End_Year,Force_End_Month,Force_End_Day, &
                         Force_Hwin_lat0,Force_Hwin_lon0,              &
                         Force_Hwin_latWidth,Force_Hwin_lonWidth,      &
                         Force_Hwin_latDelta,Force_Hwin_lonDelta,      &
                         Force_Hwin_Invert,                            &
                         Force_Vwin_Lindex,Force_Vwin_Hindex,          &
                         Force_Vwin_Ldelta,Force_Vwin_Hdelta,          &
                         Force_Vwin_Invert                            

   ! corrector is NOT initialized yet, For now
   ! corrector will always begin/end at midnight.
   !--------------------------------------------
   Force_Initialized =.false.
   Force_ON          =.false.
   Force_Beg_Sec=0
   Force_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   Force_Model         = .false.
   Force_Path          = '/n/home04/sweidman/holylfs04/IC_CESM2/'
   Force_File_Template = 'spcam_replay.%m-%d-%s.nc'
   Force_Times_Per_Day = 4
   Force_Ucoef         = 1._r8
   Force_Vcoef         = 1._r8
   Force_Qcoef         = 1._r8
   Force_Tcoef         = 1._r8
   Force_PScoef        = 0._r8
   Force_Uprof         = 1 
   Force_Vprof         = 1
   Force_Qprof         = 1
   Force_Tprof         = 1
   Force_PSprof        = 0
   Force_Beg_Year      = 1980
   Force_Beg_Month     = 1
   Force_Beg_Day       = 1
   Force_End_Year      = 2020
   Force_End_Month     = 12
   Force_End_Day       = 31
   Force_Hwin_lat0     = 0._r8
   Force_Hwin_latWidth = 9999._r8
   Force_Hwin_latDelta = 1.0_r8
   Force_Hwin_lon0     = 180._r8
   Force_Hwin_lonWidth = 9999._r8
   Force_Hwin_lonDelta = 1.0_r8
   Force_Hwin_Invert   = .false.
   Force_Hwin_lo       = 0.0_r8
   Force_Hwin_hi       = 1.0_r8
   Force_Vwin_Hindex   = float(pver+1)
   Force_Vwin_Hdelta   = 0.001_r8
   Force_Vwin_Lindex   = 0.0_r8
   Force_Vwin_Ldelta   = 0.001_r8
   Force_Vwin_Invert   = .false.
   Force_Vwin_lo       = 0.0_r8
   Force_Vwin_hi       = 1.0_r8

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'corrector_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,corrector_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('corrector_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(Force_Hwin_Invert) then
     Force_Hwin_lo = 1.0_r8
     Force_Hwin_hi = 0.0_r8
   else
     Force_Hwin_lo = 0.0_r8
     Force_Hwin_hi = 1.0_r8
   endif

   if(Force_Vwin_Invert) then
     Force_Vwin_lo = 1.0_r8
     Force_Vwin_hi = 0.0_r8
   else
     Force_Vwin_lo = 0.0_r8
     Force_Vwin_hi = 1.0_r8
   endif

   ! Check for valid namelist values 
   !----------------------------------
   if((Force_Hwin_lat0.lt.-90._r8).or.(Force_Hwin_lat0.gt.+90._r8)) then
     write(iulog,*) 'corrector: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'corrector:  Force_Hwin_lat0=',Force_Hwin_lat0
     call endrun('corrector_readnl:: ERROR in namelist')
   endif

   if((Force_Hwin_lon0.lt.0._r8).or.(Force_Hwin_lon0.ge.360._r8)) then
     write(iulog,*) 'corrector: Window lon0 must be in [0,+360)'
     write(iulog,*) 'corrector:  Force_Hwin_lon0=',Force_Hwin_lon0
     call endrun('corrector_readnl:: ERROR in namelist')
   endif

   if((Force_Vwin_Lindex.gt.Force_Vwin_Hindex)                         .or. &
      (Force_Vwin_Hindex.gt.float(pver+1)).or.(Force_Vwin_Hindex.lt.0._r8).or. &
      (Force_Vwin_Lindex.gt.float(pver+1)).or.(Force_Vwin_Lindex.lt.0._r8)   ) then
     write(iulog,*) 'corrector: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'corrector: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'corrector: Lindex must be LE than Hindex'
     write(iulog,*) 'corrector:  Force_Vwin_Lindex=',Force_Vwin_Lindex
     write(iulog,*) 'corrector:  Force_Vwin_Hindex=',Force_Vwin_Hindex
     call endrun('corrector_readnl:: ERROR in namelist')
   endif

   if((Force_Hwin_latDelta.le.0._r8).or.(Force_Hwin_lonDelta.le.0._r8).or. &
      (Force_Vwin_Hdelta  .le.0._r8).or.(Force_Vwin_Ldelta  .le.0._r8)    ) then
     write(iulog,*) 'corrector: Window Deltas must be positive'
     write(iulog,*) 'corrector:  Force_Hwin_latDelta=',Force_Hwin_latDelta
     write(iulog,*) 'corrector:  Force_Hwin_lonDelta=',Force_Hwin_lonDelta
     write(iulog,*) 'corrector:  Force_Vwin_Hdelta=',Force_Vwin_Hdelta
     write(iulog,*) 'corrector:  Force_Vwin_Ldelta=',Force_Vwin_Ldelta
     call endrun('corrector_readnl:: ERROR in namelist')

   endif

   if((Force_Hwin_latWidth.le.0._r8).or.(Force_Hwin_lonWidth.le.0._r8)) then
     write(iulog,*) 'corrector: Window widths must be positive'
     write(iulog,*) 'corrector:  Force_Hwin_latWidth=',Force_Hwin_latWidth
     write(iulog,*) 'corrector:  Force_Hwin_lonWidth=',Force_Hwin_lonWidth
     call endrun('corrector_readnl:: ERROR in namelist')
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Force_Path         ,len(Force_Path)         ,mpichar,0,mpicom)
   call mpibcast(Force_File_Template,len(Force_File_Template),mpichar,0,mpicom)
   call mpibcast(Force_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(Force_Initialized  , 1, mpilog, 0, mpicom)
   call mpibcast(Force_ON           , 1, mpilog, 0, mpicom)
   call mpibcast(Force_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Force_Ucoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Vcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Tcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Qcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_PScoef       , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Uprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Force_Vprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Force_Tprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Force_Qprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Force_PSprof       , 1, mpiint, 0, mpicom)
   call mpibcast(Force_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Force_Beg_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Force_Beg_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Force_Beg_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Force_End_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Force_End_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Force_End_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Force_End_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Force_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Force_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Force_Vwin_Invert,   1, mpilog, 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! corrector_readnl
  !================================================================


  !================================================================
  subroutine corrector_init
   ! 
   ! corrector_INIT: Allocate space and initialize corrector values
   !===============================================================
   use ppgrid        ,only: pver,pcols,begchunk,endchunk
   use error_messages,only: alloc_err
   use dycore        ,only: dycore_is
   use dyn_grid      ,only: get_horiz_grid_dim_d
   use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
   use cam_history   ,only: addfld,add_default
   use shr_const_mod ,only: SHR_CONST_PI
   use filenames     ,only: interpret_filename_spec

   ! Local values
   !----------------
   integer  Year,Month,Day,Sec
   integer  YMD1,YMD
   logical  After_Beg,Before_End
   integer  istat,lchnk,ncol,icol,ilev
   integer  hdim1_d,hdim2_d
   integer  dtime
   real(r8) rlat,rlon
   real(r8) Wprof(pver)
   real(r8) lonp,lon0,lonn,latp,lat0,latn
   real(r8) Val1_p,Val2_p,Val3_p,Val4_p
   real(r8) Val1_0,Val2_0,Val3_0,Val4_0
   real(r8) Val1_n,Val2_n,Val3_n,Val4_n
   integer               nn

   ! Get the time step size
   !------------------------
   dtime = get_step_size()

   ! Allocate Space for corrector data arrays
   !-----------------------------------------
   allocate(Target_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Target_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Target_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Target_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Target_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Target_PS',pcols*((endchunk-begchunk)+1))

   allocate(Model_state_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Model_state_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_state_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Model_state_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_state_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Model_state_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_state_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Model_state_Q',pcols*pver*((endchunk-begchunk)+1))
    allocate(Model_state_QLIQ(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_QLIQ',pcols*pver*((endchunk-begchunk)+1))
    allocate(Model_state_QICE(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_QICE',pcols*pver*((endchunk-begchunk)+1))
    allocate(Model_state_OMEGA(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_OMEGA',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_state_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Model_state_PS',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_SOLIN(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_SOLIN',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_LHFLX(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_LHFLX',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_SHFLX(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_SHFLX',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_SNOWHLND(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_SNOWHLND',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_PHIS(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_PHIS',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_TAUX(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_TAUX',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_TAUY(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_TAUY',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_TS(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_TS',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_ICEFRAC(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_ICEFRAC',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_LANDFRAC(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_LANDFRAC',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_lat(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_lat',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_lon(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_lon',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_tod(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_tod',pcols*((endchunk-begchunk)+1))
    allocate(Model_state_toy(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Model_state_toy',pcols*((endchunk-begchunk)+1))


   ! Allocate Space for spatial dependence of 
   ! corrector Coefs and corrector Forcing.
   !-------------------------------------------
   allocate(Force_Utau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_Utau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Force_Vtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_Vtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Force_Stau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_Stau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Force_Qtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_Qtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Force_PStau(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_PStau',pcols*((endchunk-begchunk)+1))

   allocate(Force_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_Ustep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Force_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_Vstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Force_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_Sstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Force_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_Qstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Force_PSstep(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'corrector_init','Force_PSstep',pcols*((endchunk-begchunk)+1))

   ! Register output fields with the cam history module
   !-----------------------------------------------------
   call addfld( 'Force_U',(/ 'lev' /),'I','m/s/s'  ,'U corrector Tendency')
   call addfld( 'Force_V',(/ 'lev' /),'I','m/s/s'  ,'V corrector Tendency')
   call addfld( 'Force_T',(/ 'lev' /),'I','K/s'    ,'T corrector Tendency')
   call addfld( 'Force_Q',(/ 'lev' /),'I','kg/kg/s','Q corrector Tendency')

   ! Values initialized only by masterproc
   !-----------------------------------------
   if(masterproc) then

     ! Set the Stepping intervals for Model and corrector values
     ! Ensure that the Model_Step is not smaller then one timestep
     !  and not larger then the Force_Step.
     !--------------------------------------------------------
     Force_Step=86400/Force_Times_Per_Day

     ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     Force_nlon=hdim1_d
     Force_nlat=hdim2_d
     Force_ncol=hdim1_d*hdim2_d
     Force_nlev=pver

     ! Check the time relative to the corrector window
     !------------------------------------------------
     call get_curr_date(Year,Month,Day,Sec)
     YMD=(Year*10000) + (Month*100) + Day
     YMD1=(Force_Beg_Year*10000) + (Force_Beg_Month*100) + Force_Beg_Day
     call timemgr_time_ge(YMD1,Force_Beg_Sec,         &
                          YMD ,Sec          ,After_Beg)
     YMD1=(Force_End_Year*10000) + (Force_End_Month*100) + Force_End_Day
     call timemgr_time_ge(YMD ,Sec          ,          &
                          YMD1,Force_End_Sec,Before_End)
  
     if((After_Beg).and.(Before_End)) then
       ! Set Time indicies so that the next call to 
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Force_Next_Year =Year
       Force_Next_Month=Month
       Force_Next_Day  =Day
       Force_Next_Sec  =(Sec/Force_Step)*Force_Step
     elseif(.not.After_Beg) then
       ! Set Time indicies to corrector start,
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Force_Next_Year =Force_Beg_Year
       Force_Next_Month=Force_Beg_Month
       Force_Next_Day  =Force_Beg_Day
       Force_Next_Sec  =Force_Beg_Sec
     elseif(.not.Before_End) then
       ! corrector will never occur, so switch it off
       !--------------------------------------------
       Force_Model=.false.
       Force_ON   =.false.
       write(iulog,*) ' '
       write(iulog,*) 'corrector: WARNING - corrector has been requested by it will'
       write(iulog,*) 'corrector:           never occur for the given time values'
       write(iulog,*) ' '
     endif

     ! Initialize values for window function  
     !----------------------------------------
     lonp= 180._r8
     lon0=   0._r8
     lonn=-180._r8
     latp=  90._r8-Force_Hwin_lat0
     lat0=   0._r8
     latn= -90._r8-Force_Hwin_lat0
    
     Force_Hwin_lonWidthH=Force_Hwin_lonWidth/2._r8
     Force_Hwin_latWidthH=Force_Hwin_latWidth/2._r8

     Val1_p=(1._r8+tanh((Force_Hwin_lonWidthH+lonp)/Force_Hwin_lonDelta))/2._r8
     Val2_p=(1._r8+tanh((Force_Hwin_lonWidthH-lonp)/Force_Hwin_lonDelta))/2._r8
     Val3_p=(1._r8+tanh((Force_Hwin_latWidthH+latp)/Force_Hwin_latDelta))/2._r8
     Val4_p=(1._r8+tanh((Force_Hwin_latWidthH-latp)/Force_Hwin_latDelta))/2_r8
     Val1_0=(1._r8+tanh((Force_Hwin_lonWidthH+lon0)/Force_Hwin_lonDelta))/2._r8
     Val2_0=(1._r8+tanh((Force_Hwin_lonWidthH-lon0)/Force_Hwin_lonDelta))/2._r8
     Val3_0=(1._r8+tanh((Force_Hwin_latWidthH+lat0)/Force_Hwin_latDelta))/2._r8
     Val4_0=(1._r8+tanh((Force_Hwin_latWidthH-lat0)/Force_Hwin_latDelta))/2._r8

     Val1_n=(1._r8+tanh((Force_Hwin_lonWidthH+lonn)/Force_Hwin_lonDelta))/2._r8
     Val2_n=(1._r8+tanh((Force_Hwin_lonWidthH-lonn)/Force_Hwin_lonDelta))/2._r8
     Val3_n=(1._r8+tanh((Force_Hwin_latWidthH+latn)/Force_Hwin_latDelta))/2._r8
     Val4_n=(1._r8+tanh((Force_Hwin_latWidthH-latn)/Force_Hwin_latDelta))/2._r8

     Force_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
     Force_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

     Force_File_Present=.false.

     ! Initialization is done, 
     !--------------------------
     Force_Initialized=.true.

     ! Check that this is a valid DYCORE model
     !------------------------------------------
     if(.not.dycore_is('LR')) then
       call endrun('corrector IS CURRENTLY ONLY CONFIGURED FOR FV')
     endif

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) '  MODEL corrector INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) 'corrector: Force_Model=',Force_Model
     write(iulog,*) 'corrector: Force_Path=',Force_Path
     write(iulog,*) 'corrector: Force_File_Template =',Force_File_Template
     write(iulog,*) 'corrector: Force_Times_Per_Day=',Force_Times_Per_Day
     write(iulog,*) 'corrector: Force_Step=',Force_Step
     write(iulog,*) 'corrector: Force_Ucoef  =',Force_Ucoef
     write(iulog,*) 'corrector: Force_Vcoef  =',Force_Vcoef
     write(iulog,*) 'corrector: Force_Qcoef  =',Force_Qcoef
     write(iulog,*) 'corrector: Force_Tcoef  =',Force_Tcoef
     write(iulog,*) 'corrector: Force_PScoef =',Force_PScoef
     write(iulog,*) 'corrector: Force_Uprof  =',Force_Uprof
     write(iulog,*) 'corrector: Force_Vprof  =',Force_Vprof
     write(iulog,*) 'corrector: Force_Qprof  =',Force_Qprof
     write(iulog,*) 'corrector: Force_Tprof  =',Force_Tprof
     write(iulog,*) 'corrector: Force_PSprof =',Force_PSprof
     write(iulog,*) 'corrector: Force_Beg_Year =',Force_Beg_Year
     write(iulog,*) 'corrector: Force_Beg_Month=',Force_Beg_Month
     write(iulog,*) 'corrector: Force_Beg_Day  =',Force_Beg_Day
     write(iulog,*) 'corrector: Force_End_Year =',Force_End_Year
     write(iulog,*) 'corrector: Force_End_Month=',Force_End_Month
     write(iulog,*) 'corrector: Force_End_Day  =',Force_End_Day
     write(iulog,*) 'corrector: Force_Hwin_lat0     =',Force_Hwin_lat0
     write(iulog,*) 'corrector: Force_Hwin_latWidth =',Force_Hwin_latWidth
     write(iulog,*) 'corrector: Force_Hwin_latDelta =',Force_Hwin_latDelta
     write(iulog,*) 'corrector: Force_Hwin_lon0     =',Force_Hwin_lon0
     write(iulog,*) 'corrector: Force_Hwin_lonWidth =',Force_Hwin_lonWidth
     write(iulog,*) 'corrector: Force_Hwin_lonDelta =',Force_Hwin_lonDelta
     write(iulog,*) 'corrector: Force_Hwin_Invert   =',Force_Hwin_Invert  
     write(iulog,*) 'corrector: Force_Hwin_lo       =',Force_Hwin_lo
     write(iulog,*) 'corrector: Force_Hwin_hi       =',Force_Hwin_hi
     write(iulog,*) 'corrector: Force_Vwin_Hindex   =',Force_Vwin_Hindex
     write(iulog,*) 'corrector: Force_Vwin_Hdelta   =',Force_Vwin_Hdelta
     write(iulog,*) 'corrector: Force_Vwin_Lindex   =',Force_Vwin_Lindex
     write(iulog,*) 'corrector: Force_Vwin_Ldelta   =',Force_Vwin_Ldelta
     write(iulog,*) 'corrector: Force_Vwin_Invert   =',Force_Vwin_Invert  
     write(iulog,*) 'corrector: Force_Vwin_lo       =',Force_Vwin_lo
     write(iulog,*) 'corrector: Force_Vwin_hi       =',Force_Vwin_hi
     write(iulog,*) 'corrector: Force_Hwin_latWidthH=',Force_Hwin_latWidthH
     write(iulog,*) 'corrector: Force_Hwin_lonWidthH=',Force_Hwin_lonWidthH
     write(iulog,*) 'corrector: Force_Hwin_max      =',Force_Hwin_max
     write(iulog,*) 'corrector: Force_Hwin_min      =',Force_Hwin_min
     write(iulog,*) 'corrector: Force_Initialized   =',Force_Initialized

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Force_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Force_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Force_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Force_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Force_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Force_Model         ,            1, mpilog, 0, mpicom)
   call mpibcast(Force_ON            ,            1, mpilog, 0, mpicom)
   call mpibcast(Force_Initialized   ,            1, mpilog, 0, mpicom)
   call mpibcast(Force_ncol          ,            1, mpiint, 0, mpicom)
   call mpibcast(Force_nlev          ,            1, mpiint, 0, mpicom)
   call mpibcast(Force_nlon          ,            1, mpiint, 0, mpicom)
   call mpibcast(Force_nlat          ,            1, mpiint, 0, mpicom)
   call mpibcast(Force_Hwin_max      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_min      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_lonWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Force_Hwin_latWidthH,            1, mpir8 , 0, mpicom)
#endif

!!DIAG
   if(masterproc) then
     write(iulog,*) 'corrector: corrector_init() SIZE#',(9*pcols*pver*((endchunk-begchunk)+1)*1)
     write(iulog,*) 'corrector: corrector_init() MB:',float(8*9*pcols*pver*((endchunk-begchunk)+1)*1)/(1024._r8*1024._r8)
     write(iulog,*) 'corrector: corrector_init() pcols=',pcols,' pver=',pver
     write(iulog,*) 'corrector: corrector_init() begchunk:',begchunk,' endchunk=',endchunk
     write(iulog,*) 'corrector: corrector_init() Force_File_Present=',Force_File_Present
   endif
!!DIAG

   ! Initialize the analysis filename at the NEXT time for startup.
   !---------------------------------------------------------------
   Force_File=interpret_filename_spec(Force_File_Template      , &
                                       yr_spec=Force_Beg_Year , &
                                      mon_spec=Force_Beg_Month, &
                                      day_spec=Force_Beg_Day  , &
                                      sec_spec=Force_Beg_Sec    )

   if(masterproc) then
    write(iulog,*) 'corrector: Reading forcing:',trim(Force_Path)//trim(Force_File)
   endif

  call corrector_update_analyses_fv (trim(Force_Path)//trim(Force_File))

   ! Initialize corrector Coeffcient profiles in local arrays
   ! Load zeros into corrector arrays
   !------------------------------------------------------
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call corrector_set_profile(rlat,rlon,Force_Uprof,Wprof,pver)
       Force_Utau(icol,:,lchnk)=Wprof(:)
       call corrector_set_profile(rlat,rlon,Force_Vprof,Wprof,pver)
       Force_Vtau(icol,:,lchnk)=Wprof(:)
       call corrector_set_profile(rlat,rlon,Force_Tprof,Wprof,pver)
       Force_Stau(icol,:,lchnk)=Wprof(:)
       call corrector_set_profile(rlat,rlon,Force_Qprof,Wprof,pver)
       Force_Qtau(icol,:,lchnk)=Wprof(:)

       Force_PStau(icol,lchnk)=corrector_set_PSprofile(rlat,rlon,Force_PSprof)
     end do
     Force_Utau(:ncol,:pver,lchnk) =                             &
     Force_Utau(:ncol,:pver,lchnk) * Force_Ucoef/float(Force_Step)
     Force_Vtau(:ncol,:pver,lchnk) =                             &
     Force_Vtau(:ncol,:pver,lchnk) * Force_Vcoef/float(Force_Step)
     Force_Stau(:ncol,:pver,lchnk) =                             &
     Force_Stau(:ncol,:pver,lchnk) * Force_Tcoef/float(Force_Step)
     Force_Qtau(:ncol,:pver,lchnk) =                             &
     Force_Qtau(:ncol,:pver,lchnk) * Force_Qcoef/float(Force_Step)
     Force_PStau(:ncol,lchnk)=                             &
     Force_PStau(:ncol,lchnk)* Force_PScoef/float(Force_Step)

     Force_Ustep(:pcols,:pver,lchnk)=0._r8
     Force_Vstep(:pcols,:pver,lchnk)=0._r8
     Force_Sstep(:pcols,:pver,lchnk)=0._r8
     Force_Qstep(:pcols,:pver,lchnk)=0._r8
     Force_PSstep(:pcols,lchnk)=0._r8
     Target_U(:pcols,:pver,lchnk)=0._r8
     Target_V(:pcols,:pver,lchnk)=0._r8
     Target_S(:pcols,:pver,lchnk)=0._r8
     Target_Q(:pcols,:pver,lchnk)=0._r8
     Target_PS(:pcols,lchnk)=0._r8

      Model_state_U(:pcols,:pver,lchnk)=0._r8
      Model_state_V(:pcols,:pver,lchnk)=0._r8
      Model_state_T(:pcols,:pver,lchnk)=0._r8
      Model_state_Q(:pcols,:pver,lchnk)=0._r8
      Model_state_QLIQ(:pcols,:pver,lchnk)=0._r8
      Model_state_QICE(:pcols,:pver,lchnk)=0._r8
      Model_state_OMEGA(:pcols,:pver,lchnk)=0._r8

      Model_state_PS(:pcols,lchnk)=0._r8
      Model_state_SOLIN(:pcols,lchnk)=0._r8
      Model_state_LHFLX(:pcols,lchnk)=0._r8
      Model_state_SHFLX(:pcols,lchnk)=0._r8
      Model_state_SNOWHLND(:pcols,lchnk)=0._r8
      Model_state_PHIS(:pcols,lchnk)=0._r8
      Model_state_TAUX(:pcols,lchnk)=0._r8
      Model_state_TAUY(:pcols,lchnk)=0._r8
      Model_state_TS(:pcols,lchnk)=0._r8
      Model_state_ICEFRAC(:pcols,lchnk)=0._r8
      Model_state_LANDFRAC(:pcols,lchnk)=0._r8
      Model_state_lat(:pcols,lchnk)=0._r8
      Model_state_lon(:pcols,lchnk)=0._r8
      Model_state_tod(:pcols,lchnk)=0._r8
      Model_state_toy(:pcols,lchnk)=0._r8

   end do

   ! End Routine
   !------------
   return
  end subroutine ! corrector_init
  !================================================================



  !================================================================
  subroutine corrector_timestep_init(phys_state)
   ! 
   ! corrector_TIMESTEP_INIT: 
   !                 Check the current time and update corrector 
   !                 arrays when necessary. Toggle the corrector flag
   !                 when the time is withing the corrector window.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind
   use dycore       ,only: dycore_is
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use filenames    ,only: interpret_filename_spec
   use ESMF

   ! Arguments
   !-----------
   type(physics_state),intent(in):: phys_state(begchunk:endchunk)

   ! Local values
   !----------------
   integer Year,Month,Day,Sec
   integer YMD1,YMD2,YMD
   logical Update_Force,Sync_Error
   logical After_Beg   ,Before_End
   integer lchnk,ncol,indw

   type(ESMF_Time)         Date1,Date2
   type(ESMF_TimeInterval) DateDiff
   integer                 DeltaT
   real(r8)                Tscale
   real(r8)                Tfrac
   integer                 rc
   integer                 nn
   integer                 kk
   real(r8)                Sbar,Qbar,Wsum
   integer                 dtime

   ! Check if corrector is initialized
   !---------------------------------
   if(.not.Force_Initialized) then
     call endrun('corrector_timestep_init:: corrector NOT Initialized')
   endif

   ! Get time step size
   !--------------------
   dtime = get_step_size()

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(Force_Beg_Year*10000) + (Force_Beg_Month*100) + Force_Beg_Day
   call timemgr_time_ge(YMD1,Force_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(Force_End_Year*10000) + (Force_End_Month*100) + Force_End_Day
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Force_End_Sec,Before_End)

   !----------------------------------------------------------------
   ! When past the NEXT time, Update corrector Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Force_Next_Year*10000) + (Force_Next_Month*100) + Force_Next_Day
   call timemgr_time_ge(YMD1,Force_Next_Sec,            &
                        YMD ,Sec           ,Update_Force)

   if((Before_End).and.(Update_Force)) then

     ! Increment the Force times by the current interval
      !---------------------------------------------------
      Force_Curr_Year =Force_Next_Year
      Force_Curr_Month=Force_Next_Month
      Force_Curr_Day  =Force_Next_Day
      Force_Curr_Sec  =Force_Next_Sec
      YMD1=(Force_Curr_Year*10000) + (Force_Curr_Month*100) + Force_Curr_Day
      call timemgr_time_inc(YMD1,Force_Curr_Sec,              &
                            YMD2,Force_Next_Sec,Force_Step,0,0)
      Force_Next_Year =(YMD2/10000)
      YMD2            = YMD2-(Force_Next_Year*10000)
      Force_Next_Month=(YMD2/100)
      Force_Next_Day  = YMD2-(Force_Next_Month*100)

     ! Read the analysis file from the current timestep
     !---------------------------------------------------------------
     Force_File=interpret_filename_spec(Force_File_Template      , &
                                         yr_spec=Force_Curr_Year , &
                                        mon_spec=Force_Curr_Month, &
                                        day_spec=Force_Curr_Day  , &
                                        sec_spec=Force_Curr_Sec    )
     if(masterproc) then
      write(iulog,*) 'corrector: Reading forcing:',trim(Force_Path)//trim(Force_File)
     endif

       call corrector_update_analyses_fv (trim(Force_Path)//trim(Force_File))

   endif ! ((Before_End).and.(Update_Force)) then

   !----------------------------------------------------------------
   ! Toggle corrector flag when the time interval is between 
   ! beginning and ending times, and all of the analyses files exist.
   !----------------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
       Force_ON=Force_File_Present
   else
     Force_ON=.false.
   endif

   !---------------------------------------------------
   ! If Data arrays have changed update stepping arrays
   !---------------------------------------------------
   if((Before_End).and.(Update_Force)) then

     ! Update the corrector tendencies
     !--------------------------------
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Force_Ustep(:ncol,:pver,lchnk)=Target_U(:ncol,:pver,lchnk)*Force_Utau(:ncol,:pver,lchnk)
       Force_Vstep(:ncol,:pver,lchnk)=Target_V(:ncol,:pver,lchnk)*Force_Vtau(:ncol,:pver,lchnk)
       Force_Sstep(:ncol,:pver,lchnk)=Target_S(:ncol,:pver,lchnk)*Force_Stau(:ncol,:pver,lchnk)
       Force_Qstep(:ncol,:pver,lchnk)=Target_Q(:ncol,:pver,lchnk)*Force_Qtau(:ncol,:pver,lchnk)
       Force_PSstep(:ncol,     lchnk)=Target_PS(:ncol,lchnk)*Force_PStau(:ncol,lchnk)
     end do

     if (masterproc) then
        write(iulog,*)  'day, sec', Force_Curr_Day, Force_Curr_Sec
        write(iulog,*) 'Force_Utau(1,20,1) = ', Force_Utau(1,20,begchunk)
        write(iulog,*) 'Target_U(1,20,1) = ', Target_U(1,20,begchunk)
        write(iulog,*) 'Force_Ustep(1,20,1) = ', Force_Ustep(1,20,begchunk)
     end if

   endif ! ((Before_End).and.(Update_Force)) then

   ! End Routine
   !------------
   return
  end subroutine ! corrector_timestep_init
  !================================================================

!================================================================
  subroutine nncorrector_timestep_init(phys_state, cam_in)
    ! 
    ! nncorrector_TIMESTEP_INIT: 
    ! Zeyuan Hu 12/23/2024: 
    !                 This subroutine is inherited from corrector_TIMESTEP_INIT subroutine to use the neural network to update bias correctors
    !                 Check the current time and update corrector 
    !                 arrays when necessary. Toggle the corrector flag
    !                 when the time is withing the corrector window.
    !===============================================================
    use physconst    ,only: cpair
    use physics_types,only: physics_state
    use constituents ,only: cnst_get_ind
    use dycore       ,only: dycore_is
    use ppgrid       ,only: pver,pcols,begchunk,endchunk
    use filenames    ,only: interpret_filename_spec
    use ESMF
    use camsrfexch     ,only: cam_in_t,cam_out_t
 
    ! Arguments
    !-----------
    type(physics_state),intent(in):: phys_state(begchunk:endchunk)
    type(cam_in_t),intent(in):: cam_in(begchunk:endchunk)
 
    ! Local values
    !----------------
    integer Year,Month,Day,Sec
    integer YMD1,YMD2,YMD
    logical Update_Force,Sync_Error
    logical After_Beg   ,Before_End
    integer lchnk,ncol,indw
 
    type(ESMF_Time)         Date1,Date2
    type(ESMF_TimeInterval) DateDiff
    integer                 DeltaT
    real(r8)                Tscale
    real(r8)                Tfrac
    integer                 rc
    integer                 nn
    integer                 kk
    real(r8)                Sbar,Qbar,Wsum
    integer                 dtime
 
    ! Check if corrector is initialized
    !---------------------------------
    if(.not.Force_Initialized) then
      call endrun('nncorrector_timestep_init:: corrector NOT Initialized')
    endif
 
    ! Get time step size
    !--------------------
    dtime = get_step_size()
 
    ! Get Current time
    !--------------------
    call get_curr_date(Year,Month,Day,Sec)
    YMD=(Year*10000) + (Month*100) + Day
 
    !-------------------------------------------------------
    ! Determine if the current time is AFTER the begining time
    ! and if it is BEFORE the ending time.
    !-------------------------------------------------------
    YMD1=(Force_Beg_Year*10000) + (Force_Beg_Month*100) + Force_Beg_Day
    call timemgr_time_ge(YMD1,Force_Beg_Sec,         &
                         YMD ,Sec          ,After_Beg)
 
    YMD1=(Force_End_Year*10000) + (Force_End_Month*100) + Force_End_Day
    call timemgr_time_ge(YMD ,Sec,                    &
                         YMD1,Force_End_Sec,Before_End)
 
    !----------------------------------------------------------------
    ! When past the NEXT time, Update corrector Arrays and time indices
    !----------------------------------------------------------------
    YMD1=(Force_Next_Year*10000) + (Force_Next_Month*100) + Force_Next_Day
    call timemgr_time_ge(YMD1,Force_Next_Sec,            &
                         YMD ,Sec           ,Update_Force)
 
    if((Before_End).and.(Update_Force)) then
 
      ! Increment the Force times by the current interval
       !---------------------------------------------------
       Force_Curr_Year =Force_Next_Year
       Force_Curr_Month=Force_Next_Month
       Force_Curr_Day  =Force_Next_Day
       Force_Curr_Sec  =Force_Next_Sec
       YMD1=(Force_Curr_Year*10000) + (Force_Curr_Month*100) + Force_Curr_Day
       call timemgr_time_inc(YMD1,Force_Curr_Sec,              &
                             YMD2,Force_Next_Sec,Force_Step,0,0)
       Force_Next_Year =(YMD2/10000)
       YMD2            = YMD2-(Force_Next_Year*10000)
       Force_Next_Month=(YMD2/100)
       Force_Next_Day  = YMD2-(Force_Next_Month*100)
 
      ! Read the analysis file from the current timestep
      !---------------------------------------------------------------
      Force_File=interpret_filename_spec(Force_File_Template      , &
                                          yr_spec=Force_Curr_Year , &
                                         mon_spec=Force_Curr_Month, &
                                         day_spec=Force_Curr_Day  , &
                                         sec_spec=Force_Curr_Sec    )
      if(masterproc) then
       write(iulog,*) 'corrector: Reading forcing:',trim(Force_Path)//trim(Force_File)
      endif
 
        ! call corrector_update_analyses_fv (trim(Force_Path)//trim(Force_File))
      call nncorrector_update(trim(Force_Path)//trim(Force_File), phys_state, cam_in)
 
    endif ! ((Before_End).and.(Update_Force)) then
 
    !----------------------------------------------------------------
    ! Toggle corrector flag when the time interval is between 
    ! beginning and ending times, and all of the analyses files exist.
    !----------------------------------------------------------------
    if((After_Beg).and.(Before_End)) then
        Force_ON=Force_File_Present
      !Force_ON = .true. ! always turn on corrector when it is within the time window ! Zeyuan Hu 12/23/2024
    else
      Force_ON=.false.
    endif
 
    !---------------------------------------------------
    ! If Data arrays have changed update stepping arrays
    !---------------------------------------------------
    if((Before_End).and.(Update_Force)) then
 
      ! Update the corrector tendencies
      !--------------------------------
      do lchnk=begchunk,endchunk
        ncol=phys_state(lchnk)%ncol
        Force_Ustep(:ncol,:pver,lchnk)=Target_U(:ncol,:pver,lchnk)*Force_Utau(:ncol,:pver,lchnk)
        Force_Vstep(:ncol,:pver,lchnk)=Target_V(:ncol,:pver,lchnk)*Force_Vtau(:ncol,:pver,lchnk)
        Force_Sstep(:ncol,:pver,lchnk)=Target_S(:ncol,:pver,lchnk)*Force_Stau(:ncol,:pver,lchnk)
        Force_Qstep(:ncol,:pver,lchnk)=Target_Q(:ncol,:pver,lchnk)*Force_Qtau(:ncol,:pver,lchnk)
        Force_PSstep(:ncol,     lchnk)=Target_PS(:ncol,lchnk)*Force_PStau(:ncol,lchnk)
      end do
 
      if (masterproc) then
         write(iulog,*)  'day, sec', Force_Curr_Day, Force_Curr_Sec
         write(iulog,*) 'Force_Utau(1,20,1) = ', Force_Utau(1,20,begchunk)
         write(iulog,*) 'Target_U(1,20,1) = ', Target_U(1,20,begchunk)
         write(iulog,*) 'Force_Ustep(1,20,1) = ', Force_Ustep(1,20,begchunk)
      end if
 
    endif ! ((Before_End).and.(Update_Force)) then
 
    ! End Routine
    !------------
    return
   end subroutine ! nncorrector_timestep_init

  !================================================================
  subroutine corrector_timestep_tend(phys_state,phys_tend)
   ! 
   ! corrector_TIMESTEP_TEND: 
   !                If corrector is ON, return the corrector contributions 
   !                to forcing using the current contents of the Nudge 
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk
   logical lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'corrector',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if(Force_ON) then
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
     phys_tend%u(:ncol,:pver)     =Force_Ustep(:ncol,:pver,lchnk)
     phys_tend%v(:ncol,:pver)     =Force_Vstep(:ncol,:pver,lchnk)
     phys_tend%s(:ncol,:pver)     =Force_Sstep(:ncol,:pver,lchnk)
     phys_tend%q(:ncol,:pver,indw)=Force_Qstep(:ncol,:pver,lchnk)

     if (masterproc) then
      write(iulog,*) 'phys_tend%u = ', phys_tend%u(1,20)
     end if

     call outfld( 'Force_U',phys_tend%u                ,pcols,lchnk)
     call outfld( 'Force_V',phys_tend%v                ,pcols,lchnk)
     call outfld( 'Force_T',phys_tend%s/cpair          ,pcols,lchnk)
     call outfld( 'Force_Q',phys_tend%q(1,1,indw)      ,pcols,lchnk)

   endif

   ! End Routine
   !------------
   return
  end subroutine ! corrector_timestep_tend
  !================================================================


  !================================================================
  subroutine corrector_update_analyses_fv(anal_file)
   ! 
   ! corrector_UPDATE_ANALYSES_FV: 
   !                 Open the given analyses data file, read in 
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Force_nlon,Force_nlat,Force_nlev)
   real(r8) PSanal(Force_nlon,Force_nlat)
   real(r8) Lat_anal(Force_nlat)
   real(r8) Lon_anal(Force_nlon)
   real(r8) Xtrans(Force_nlon,Force_nlev,Force_nlat)
   integer  nn,Nindex

   ! Rotate Force_ObsInd() indices, then check the existence of the analyses 
   ! file; broadcast the updated indices and file status to all the other MPI nodes. 
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_file),EXIST=Force_File_Present)
     write(iulog,*)'corrector: Force_File_Present=',Force_File_Present
   endif
#ifdef SPMD
   call mpibcast(Force_File_Present, 1, mpilog, 0, mpicom)
#endif
   if(.not.Force_File_Present) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then
   
     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
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

     if((Force_nlon.ne.nlon).or.(Force_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: corrector_update_analyses_fv: nlon=',nlon,' Force_nlon=',Force_nlon
      write(iulog,*) 'ERROR: corrector_update_analyses_fv: nlat=',nlat,' Force_nlat=',Force_nlat
      write(iulog,*) 'ERROR: corrector_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('corrector_update_analyses_fv: analyses dimension mismatch')
     endif

     ! Read in, transpose lat/lev indices, 
     ! and scatter data arrays
     !----------------------------------
     istat=nf90_inq_varid(ncid,'UDIFF',varid)
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
   call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
                               Target_U(1,1,begchunk))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'VDIFF',varid)
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
   call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
                               Target_V(1,1,begchunk))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'SDIFF',varid)
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
   call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
                               Target_S(1,1,begchunk))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'QDIFF',varid)
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
   call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
                               Target_Q(1,1,begchunk))

   ! End Routine
   !------------
   return
  end subroutine ! corrector_update_analyses_fv
  !================================================================


  subroutine nncorrector_update(anal_file, phys_state, cam_in)
    ! 
    ! nncorrector_UPDATE: 
    !                 generate NN predicted bias correctors of 
    !                 U,V,T,Q, and PS values and then distribute
    !                 the values to all of the chunks.
    !===============================================================
    use ppgrid ,only: pver,pcols,begchunk,endchunk
    use netcdf
    use constituents ,only: cnst_get_ind
    use physics_types,only: physics_state
    use camsrfexch     ,only: cam_in_t,cam_out_t
    use radconstants,     only: nswbands, get_ref_solar_band_irrad
    use rad_solar_var,    only: get_variability
    use time_manager,     only: get_curr_calday
    use phys_grid,        only: get_rlat_all_p, get_rlon_all_p
    use cam_control_mod,  only: lambm0, obliqr, eccen, mvelpp
    use shr_orb_mod,      only: shr_orb_decl, shr_orb_cosz

    ! Arguments
    !-------------
    character(len=*),intent(in):: anal_file
    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
    type(cam_in_t),intent(in):: cam_in(begchunk:endchunk)

    ! real(r8),allocatable::Model_state_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_QLIQ     (:,:,:)  !(pcols,pver,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_QICE     (:,:,:)  !(pcols,pver,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_OMEGA     (:,:,:)  !(pcols,pver,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_PS    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_SOLIN    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_LHFLX    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_SHFLX    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_SNOWHLND    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_PHIS    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_TAUX    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_TAUY    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_TS    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_ICEFRAC    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_LANDFRAC    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_lat    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_lon    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_tod    (:,:)    !(pcols,begchunk:endchunk)
    ! real(r8),allocatable::Model_state_toy    (:,:)    !(pcols,begchunk:endchunk)

    ! Local values
    !-------------
    integer lev
    integer nlon,nlat,plev,istat
    integer ncid,varid
    integer ilat,ilon,ilev
    real(r8) Xanal(Force_nlon,Force_nlat,Force_nlev)
    real(r8) Uanal(Force_nlon,Force_nlat,Force_nlev)
    real(r8) Vanal(Force_nlon,Force_nlat,Force_nlev)
    real(r8) Tanal(Force_nlon,Force_nlat,Force_nlev)
    real(r8) Qanal(Force_nlon,Force_nlat,Force_nlev)
    real(r8) QLIQanal(Force_nlon,Force_nlat,Force_nlev)
    real(r8) QICEanal(Force_nlon,Force_nlat,Force_nlev)
    real(r8) OMEGAanal(Force_nlon,Force_nlat,Force_nlev)
    real(r8) PSanal(Force_nlon,Force_nlat)
    real(r8) SOLINanal(Force_nlon,Force_nlat)
    real(r8) LHFLXanal(Force_nlon,Force_nlat)
    real(r8) SHFLXanal(Force_nlon,Force_nlat)
    real(r8) SNOWHLNDanal(Force_nlon, Force_nlat)
    real(r8) PHISanal(Force_nlon, Force_nlat)
    real(r8) TAUXanal(Force_nlon, Force_nlat)
    real(r8) TAUYanal(Force_nlon, Force_nlat)
    real(r8) TSanal(Force_nlon, Force_nlat)
    real(r8) ICEFRACanal(Force_nlon, Force_nlat)
    real(r8) LANDFRACanal(Force_nlon, Force_nlat)
    real(r8) tod_anal(Force_nlon, Force_nlat)
    real(r8) toy_anal(Force_nlon, Force_nlat)
    real(r8) clat_anal(Force_nlon, Force_nlat)
    real(r8) clon_anal(Force_nlon, Force_nlat)

    real(r8) Lat_anal(Force_nlat)
    real(r8) Lon_anal(Force_nlon)
    real(r8) Xtrans(Force_nlon,Force_nlev,Force_nlat)
    real(r8) Xtransf(1,Force_nlon,Force_nlat)
    integer  nn,Nindex
    integer  lchnk,ncol,indw, ixcldice,ixcldliq
    real(r8) pi

    real(r8) :: sfac(1:nswbands)  ! time varying scaling factors due to Solar Spectral Irrad at 1 A.U. per band
    real(r8) :: solar_band_irrad(1:nswbands) ! rrtmg-assumed solar irradiance in each sw band
    real(r8) :: delta    ! Solar declination angle  in radians
    real(r8) :: dt_avg = 0.0_r8   ! time step to use for the shr_orb_cosz calculation, if use_rad_dt_cosz set to true
    real(r8) :: eccf     ! Earth orbit eccentricity factor
    real(r8) :: calday       ! current calendar day
    real(r8) :: clat(pcols)  ! current latitudes(radians)
    real(r8) :: clon(pcols)  ! current longitudes(radians)
    real(r8), dimension(pcols,begchunk:endchunk) :: coszrs  ! Cosine solar zenith angle
    ! real(r8), dimension(pcols,begchunk:endchunk) :: solin   ! Insolation

    integer :: n,i,j,k
    type(torch_tensor_wrap) :: input_tensors
    type(torch_tensor) :: out_tensor
    real(real32) :: input_torch(144, 96, nn_inputlength, 1)
    real(real32), pointer :: output_torch(:, :, :, :)
    
    ! Data arrays for output netcdf file
    integer :: varid_input, varid_output
    integer :: dimids_input(4), dimids_output(4)
    integer :: retval
    real, dimension(1, nn_inputlength, Force_nlat, Force_nlon) :: data_input
    real, dimension(1, nn_outputlength, Force_nlat, Force_nlon) :: data_output
    character(len=100) :: nc_filename
    ! character(:), allocatable :: filename
    ! character(len=50) :: outputfile
    ! integer :: arglen, stat
    integer :: unit


    pi = 3.14159265358979323846_r8
    call cnst_get_ind('Q',indw)
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    ! Rotate Force_ObsInd() indices, then check the existence of the analyses 
    ! file; broadcast the updated indices and file status to all the other MPI nodes. 
    ! If the file is not there, then just return.
    !------------------------------------------------------------------------
    if(masterproc) then
      inquire(FILE=trim(anal_file),EXIST=Force_File_Present)
      write(iulog,*)'corrector: Force_File_Present=',Force_File_Present
    endif
 #ifdef SPMD
    call mpibcast(Force_File_Present, 1, mpilog, 0, mpicom)
 #endif
    if(.not.Force_File_Present) return
 
    ! masterporc does all of the work here
    !-----------------------------------------
    if(masterproc) then
    
      ! Open the given file
      !-----------------------
      istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
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
 
      if((Force_nlon.ne.nlon).or.(Force_nlat.ne.nlat).or.(plev.ne.pver)) then
       write(iulog,*) 'ERROR: corrector_update_analyses_fv: nlon=',nlon,' Force_nlon=',Force_nlon
       write(iulog,*) 'ERROR: corrector_update_analyses_fv: nlat=',nlat,' Force_nlat=',Force_nlat
       write(iulog,*) 'ERROR: corrector_update_analyses_fv: plev=',plev,' pver=',pver
       call endrun('corrector_update_analyses_fv: analyses dimension mismatch')
      endif
    endif ! (masterproc) 
    
    ! Zeyuan Hu 12/23/2024: gather global state variables
    !---------------------------------------------------
    do lchnk=begchunk,endchunk
        ncol=phys_state(lchnk)%ncol
        Model_state_U(:ncol,:pver,lchnk)=phys_state(lchnk)%u(:ncol,:pver)
        Model_state_V(:ncol,:pver,lchnk)=phys_state(lchnk)%v(:ncol,:pver)
        Model_state_T(:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
        Model_state_Q(:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,indw)
        Model_state_QLIQ(:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,ixcldliq)
        Model_state_QICE(:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,ixcldice)
        Model_state_OMEGA(:ncol,:pver,lchnk)=phys_state(lchnk)%omega(:ncol,:pver)
        Model_state_PS(:ncol,lchnk)=phys_state(lchnk)%ps(:ncol)
        ! Model_state_SOLIN(:ncol,lchnk)=0.0 ! Zeyuan Hu 12/23/2024: set to 0 for now
        Model_state_LHFLX(:ncol,lchnk)=cam_in(lchnk)%lhf(:ncol)
        Model_state_SHFLX(:ncol,lchnk)=cam_in(lchnk)%shf(:ncol)
        Model_state_SNOWHLND(:ncol,lchnk)=cam_in(lchnk)%snowhland(:ncol)
        Model_state_PHIS(:ncol,lchnk)=phys_state(lchnk)%phis(:ncol)
        Model_state_TAUX(:ncol,lchnk)=cam_in(lchnk)%wsx(:ncol)
        Model_state_TAUY(:ncol,lchnk)=cam_in(lchnk)%wsy(:ncol)
        Model_state_TS(:ncol,lchnk)=cam_in(lchnk)%ts(:ncol)
        Model_state_ICEFRAC(:ncol,lchnk)=cam_in(lchnk)%icefrac(:ncol)
        Model_state_LANDFRAC(:ncol,lchnk)=cam_in(lchnk)%landfrac(:ncol)
        Model_state_lat(:ncol,lchnk)=phys_state(lchnk)%lat(:ncol)*(180./pi)
        Model_state_lon(:ncol,lchnk)=phys_state(lchnk)%lon(:ncol)*(180./pi)
        Model_state_tod(:ncol,lchnk)=Force_Curr_Sec/3600. ! in hours
        Model_state_toy(:ncol,lchnk)=Force_Curr_Day ! in day
    end do

    call get_ref_solar_band_irrad( solar_band_irrad ) ! this can move to init subroutine
    call get_variability(sfac)                        ! "
    do lchnk=begchunk,endchunk
      ncol = phys_state(lchnk)%ncol
      calday = get_curr_calday() ! get current calendar day; no time offset as was in E3SM, need to double check!
      ! coszrs
      call get_rlat_all_p(lchnk, ncol, clat)
      call get_rlon_all_p(lchnk, ncol, clon)
      call shr_orb_decl(calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                        delta   ,eccf      )
      ! call zenith(calday, clat, clon, coszrs(:,lchnk), ncol, dt_avg)
      do i = 1, ncol
        coszrs(i,lchnk) = shr_orb_cosz(calday, clat(i), clon(i), delta, dt_avg)
      end do
      ! solin(:,lchnk) = sum(sfac(:)*solar_band_irrad(:)) * eccf * coszrs(:,lchnk)
      Model_state_SOLIN(:ncol,lchnk) = sum(sfac(:)*solar_band_irrad(:)) * eccf * coszrs(:,lchnk)
    end do

    call gather_chunk_to_field(1,Force_nlev,1,Force_nlon,Model_state_U,Xtrans)
    if (masterproc) then
      do ilat=1,nlat
      do ilev=1,plev
      do ilon=1,nlon
        Uanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
      end do
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,Force_nlev,1,Force_nlon,Model_state_V,Xtrans)
    if (masterproc) then
      do ilat=1,nlat
      do ilev=1,plev
      do ilon=1,nlon
        Vanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
      end do
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,Force_nlev,1,Force_nlon,Model_state_T,Xtrans)
    if (masterproc) then
      do ilat=1,nlat
      do ilev=1,plev
      do ilon=1,nlon
        Tanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
      end do
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,Force_nlev,1,Force_nlon,Model_state_Q,Xtrans)
    if (masterproc) then
      do ilat=1,nlat
      do ilev=1,plev
      do ilon=1,nlon
        Qanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
      end do
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,Force_nlev,1,Force_nlon,Model_state_QLIQ,Xtrans)
    if (masterproc) then
      do ilat=1,nlat
      do ilev=1,plev
      do ilon=1,nlon
        QLIQanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
      end do
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,Force_nlev,1,Force_nlon,Model_state_QICE,Xtrans)
    if (masterproc) then
      do ilat=1,nlat
      do ilev=1,plev
      do ilon=1,nlon
        QICEanal(ilon,ilat,ilev)=Xtrans(ilon,ilat,ilev)
      end do
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,Force_nlev,1,Force_nlon,Model_state_OMEGA,Xtrans)
    if (masterproc) then
      do ilat=1,nlat
      do ilev=1,plev
      do ilon=1,nlon
        OMEGAanal(ilon,ilat,ilev)=Xtrans(ilon,ilat,ilev)
      end do
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_PS,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        PSanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then
    
    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_SOLIN,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        SOLINanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_LHFLX,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        LHFLXanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_SHFLX,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        SHFLXanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_SNOWHLND,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        SNOWHLNDanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_PHIS,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        PHISanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_TAUX,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        TAUXanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_TAUY,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        TAUYanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_TS,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        TSanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then
    
    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_ICEFRAC,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        ICEFRACanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_LANDFRAC,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        LANDFRACanal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then
    
    ! call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_tod,Xtransf)
    ! if (masterproc) then
    !   do ilat=1,nlat
    !   do ilon=1,nlon
    !     tod_anal(ilon,ilat)=Xtransf(1,ilon,ilat)
    !   end do
    !   end do
    ! endif ! (masterproc) then

    ! call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_toy,Xtransf)
    ! if (masterproc) then
    !   do ilat=1,nlat
    !   do ilon=1,nlon
    !     toy_anal(ilon,ilat)=Xtransf(1,ilon,ilat)
    !   end do
    !   end do
    ! endif ! (masterproc) then
    tod_anal(:,:) = Force_Curr_Sec/3600. ! in hours
    toy_anal(:,:) = Force_Curr_Day ! in day
    
    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_lat,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        clat_anal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then

    call gather_chunk_to_field(1,1,1,Force_nlon,Model_state_lon,Xtransf)
    if (masterproc) then
      do ilat=1,nlat
      do ilon=1,nlon
        clon_anal(ilon,ilat)=Xtransf(1,ilon,ilat)
      end do
      end do
    endif ! (masterproc) then
    
    ! collect and prepare the input data for the neural network
    if (masterproc) then   
      do ilon=1,nlon
        do ilat=1,nlat

          input_torch(ilon,ilat,0*plev+1:1*plev,1) = Tanal(ilon,ilat,1:plev)
          input_torch(ilon,ilat,1*plev+1:2*plev,1) = Qanal(ilon,ilat,1:plev)
          input_torch(ilon,ilat,2*plev+1:3*plev,1) = Uanal(ilon,ilat,1:plev)
          input_torch(ilon,ilat,3*plev+1:4*plev,1) = Vanal(ilon,ilat,1:plev)
          input_torch(ilon,ilat,4*plev+1:5*plev,1) = QLIQanal(ilon,ilat,1:plev)
          input_torch(ilon,ilat,5*plev+1:6*plev,1) = QICEanal(ilon,ilat,1:plev)
          input_torch(ilon,ilat,6*plev+1:7*plev,1) = OMEGAanal(ilon,ilat,1:plev)
          input_torch(ilon,ilat,7*plev+1,1) = PSanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+2,1) = SOLINanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+3,1) = LHFLXanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+4,1) = SHFLXanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+5,1) = SNOWHLNDanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+6,1) = PHISanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+7,1) = TAUXanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+8,1) = TAUYanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+9,1) = TSanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+10,1) = ICEFRACanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+11,1) = LANDFRACanal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+12,1) = clat_anal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+13,1) = clon_anal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+14,1) = tod_anal(ilon,ilat)
          input_torch(ilon,ilat,7*plev+15,1) = toy_anal(ilon,ilat)

        end do
      end do
    else
      input_torch(:,:,:,:) = 0.0_r8
    endif ! (masterproc) then
    ! run the NN inference
  
    call input_tensors%create
    call input_tensors%add_array(input_torch)
    call torch_mod(1)%forward(input_tensors, out_tensor, flags=module_use_inference_mode)
    call out_tensor%to_array(output_torch)
      ! integer :: n

      ! integer :: i
      ! integer :: use_gpu
      ! type(torch_module) :: torch_mod
      ! type(torch_tensor_wrap) :: input_tensors
      ! type(torch_tensor) :: out_tensor
   
      ! real(real32) :: input(124, 1)
      ! real(real32), pointer :: output(:, :)
   
      ! ! character(:), allocatable :: filename
      ! ! character(len=50) :: outputfile
      ! ! integer :: arglen, stat
      ! integer :: unit
   
      ! unit = 20
   
      !    if (masterproc) then
      !       write(iulog, *)  "Reading input from test_input.txt"
      !    end if
      !    open(unit=unit, file="/n/home00/zeyuanhu/spcam_ml_sourcemode/test_files/test_input.txt", status="old", action="read")
      !    do i = 1, size(input, 1)
      !       read(unit,*) input(i,1)
      !    end do
      !    close(unit)
   
      !    use_gpu = 0 !module_use_device
   
      !    ! write(iulog, *) "Creating input tensor"
      !    call input_tensors%create
      !    ! write(iulog, *) "Adding input data"
      !    call input_tensors%add_array(input)
      !    ! write(iulog, *) "Loading model"
      !    call torch_mod%load("/n/home00/zeyuanhu/spcam_ml_sourcemode/test_files/final_hsr_wrapped.pt", use_gpu)
      !    ! write(iulog, *) "Running forward pass"
      !    call torch_mod%forward(input_tensors, out_tensor, flags=module_use_inference_mode)
      !    ! write(iulog, *) "Getting output data"
      !    call out_tensor%to_array(output)
         
   
      !    if (masterproc) then
      !       write(iulog, *) "torch Output data:"
      !       do i = 1, size(output, 1)
      !          write(iulog, *) output(i,1)
      !       end do
      !    end if


    ! a placeholder for now to set 0 for correctors and scatter to all chunks
    if (masterproc) then 
      do ilat=1,nlat
        do ilev=1,plev
        do ilon=1,nlon
          Xtrans(ilon,ilev,ilat)=output_torch(ilon,ilat,ilev,1)
        end do
        end do
      end do
    endif ! (masterproc) then
    call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
    Target_S(1,1,begchunk))

    if (masterproc) then 
      do ilat=1,nlat
        do ilev=1,plev
        do ilon=1,nlon
          Xtrans(ilon,ilev,ilat)=output_torch(ilon,ilat,1*plev+ilev,1)
        end do
        end do
      end do
    endif ! (masterproc) then
    call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
    Target_Q(1,1,begchunk))

    if (masterproc) then 
      do ilat=1,nlat
        do ilev=1,plev
        do ilon=1,nlon
          Xtrans(ilon,ilev,ilat)=output_torch(ilon,ilat,2*plev+ilev,1)
        end do
        end do
      end do
    endif ! (masterproc) then
    call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
    Target_U(1,1,begchunk))

    if (masterproc) then 
      do ilat=1,nlat
        do ilev=1,plev
        do ilon=1,nlon
          Xtrans(ilon,ilev,ilat)=output_torch(ilon,ilat,4*plev+ilev,1)
        end do
        end do
      end do
    endif ! (masterproc) then
    call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
    Target_V(1,1,begchunk))
    
  if (masterproc) then  
    ! output the data of input/output to a netcdf file
    do i=1,nn_inputlength
      do ilat=1,nlat
        do ilon=1,nlon
          data_input(1,i,ilat,ilon) = input_torch(ilon,ilat,i,1)
        end do
      end do
    end do

    do i=1,nn_outputlength
      do ilat=1,nlat
        do ilon=1,nlon
          data_output(1,i,ilat,ilon) = output_torch(ilon,ilat,i,1)
        end do
      end do
    end do

      ! Create filename with time information
    write(nc_filename, '(A,I4.4,A,I2.2,A,I2.2,A,I5.5,A)') &
    "nn_verification_", Force_Curr_Year, "-", &
    Force_Curr_Month, "-", Force_Curr_Day, "-", &
    Force_Curr_Sec, ".nc"

    print *, "Filename: ", trim(nc_filename)

    ! Create NetCDF file
  istat = nf90_create(trim(nc_filename), nf90_clobber, ncid)
  if (istat /= nf90_noerr) then
    write(*, *) 'NF90_CREATE: failed for file ', trim(nc_filename)
    write(*, *) nf90_strerror(istat)
    call endrun('CREATE_NETCDF')
  endif


  ! Define dimensions
  istat = nf90_def_dim(ncid, "time", 1, dimids_input(1))
  istat = nf90_def_dim(ncid, "input_length", nn_inputlength, dimids_input(2))
  istat = nf90_def_dim(ncid, "lat", nlat, dimids_input(3))
  istat = nf90_def_dim(ncid, "lon", nlon, dimids_input(4))
  if (istat /= nf90_noerr) then
    write(*, *) 'NF90_DEF_DIM: failed'
    write(*, *) nf90_strerror(istat)
    call endrun('DEFINE_DIMENSIONS')
  endif

    dimids_output = dimids_input
    istat = nf90_def_dim(ncid, "output_length", nn_outputlength, dimids_output(2))
    if (istat /= nf90_noerr) then
      write(*, *) 'NF90_DEF_DIM (output_length): failed'
      write(*, *) nf90_strerror(istat)
      call endrun('DEFINE_OUTPUT_DIMENSIONS')
    endif
      
      ! Define variables for input and output data
    istat = nf90_def_var(ncid, "data_input", nf90_real, dimids_input, varid_input)
    if (istat /= nf90_noerr) then
      write(*, *) 'NF90_DEF_VAR (data_input): failed'
      write(*, *) nf90_strerror(istat)
      call endrun('DEFINE_VAR_INPUT')
    endif

    istat = nf90_def_var(ncid, "data_output", nf90_real, dimids_output, varid_output)
    if (istat /= nf90_noerr) then
      write(*, *) 'NF90_DEF_VAR (data_output): failed'
      write(*, *) nf90_strerror(istat)
      call endrun('DEFINE_VAR_OUTPUT')
    endif

    ! End define mode
    istat = nf90_enddef(ncid)
    if (istat /= nf90_noerr) then
      write(*, *) 'NF90_ENDDEF: failed'
      write(*, *) nf90_strerror(istat)
      call endrun('ENDDEF')
    endif

    ! Write data
    istat = nf90_put_var(ncid, varid_input, data_input)
    if (istat /= nf90_noerr) then
      write(*, *) 'NF90_PUT_VAR (data_input): failed'
      write(*, *) nf90_strerror(istat)
      call endrun('PUT_VAR_INPUT')
    endif

    istat = nf90_put_var(ncid, varid_output, data_output)
    if (istat /= nf90_noerr) then
      write(*, *) 'NF90_PUT_VAR (data_output): failed'
      write(*, *) nf90_strerror(istat)
      call endrun('PUT_VAR_OUTPUT')
    endif

    ! Close file
    istat = nf90_close(ncid)
    if (istat /= nf90_noerr) then
      write(*, *) 'NF90_CLOSE: failed'
      write(*, *) nf90_strerror(istat)
      call endrun('CLOSE_NETCDF')
    endif

    print *, "Successfully written input and output arrays to ", trim(nc_filename)

  endif ! (masterproc) then

    ! if (masterproc) then
    !   ! Read in, transpose lat/lev indices, 
    !   ! and scatter data arrays
    !   !----------------------------------
    !   istat=nf90_inq_varid(ncid,'UDIFF',varid)
    !   if(istat.ne.NF90_NOERR) then
    !     write(iulog,*) nf90_strerror(istat)
    !     call endrun ('UPDATE_ANALYSES_FV')
    !   endif
    !   istat=nf90_get_var(ncid,varid,Xanal)
    !   if(istat.ne.NF90_NOERR) then
    !     write(iulog,*) nf90_strerror(istat)
    !     call endrun ('UPDATE_ANALYSES_FV')
    !   endif
    !   do ilat=1,nlat
    !   do ilev=1,plev
    !   do ilon=1,nlon
    !     Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    !   end do
    !   end do
    !   end do
    ! endif ! (masterproc) then
    ! call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
    !                             Target_U(1,1,begchunk))
 
    ! if(masterproc) then
    !   istat=nf90_inq_varid(ncid,'VDIFF',varid)
    !   if(istat.ne.NF90_NOERR) then
    !     write(iulog,*) nf90_strerror(istat)
    !     call endrun ('UPDATE_ANALYSES_FV')
    !   endif
    !   istat=nf90_get_var(ncid,varid,Xanal)
    !   if(istat.ne.NF90_NOERR) then
    !     write(iulog,*) nf90_strerror(istat)
    !     call endrun ('UPDATE_ANALYSES_FV')
    !   endif
    !   do ilat=1,nlat
    !   do ilev=1,plev
    !   do ilon=1,nlon
    !     Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    !   end do
    !   end do
    !   end do
    ! endif ! (masterproc) then
    ! call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
    !                             Target_V(1,1,begchunk))
 
    ! if(masterproc) then
    !   istat=nf90_inq_varid(ncid,'SDIFF',varid)
    !   if(istat.ne.NF90_NOERR) then
    !     write(iulog,*) nf90_strerror(istat)
    !     call endrun ('UPDATE_ANALYSES_FV')
    !   endif
    !   istat=nf90_get_var(ncid,varid,Xanal)
    !   if(istat.ne.NF90_NOERR) then
    !     write(iulog,*) nf90_strerror(istat)
    !     call endrun ('UPDATE_ANALYSES_FV')
    !   endif
    !   do ilat=1,nlat
    !   do ilev=1,plev
    !   do ilon=1,nlon
    !     Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    !   end do
    !   end do
    !   end do
    ! endif ! (masterproc) then
    ! call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
    !                             Target_S(1,1,begchunk))
 
    ! if(masterproc) then
    !   istat=nf90_inq_varid(ncid,'QDIFF',varid)
    !   if(istat.ne.NF90_NOERR) then
    !     write(iulog,*) nf90_strerror(istat)
    !     call endrun ('UPDATE_ANALYSES_FV')
    !   endif
    !   istat=nf90_get_var(ncid,varid,Xanal)
    !   if(istat.ne.NF90_NOERR) then
    !     write(iulog,*) nf90_strerror(istat)
    !     call endrun ('UPDATE_ANALYSES_FV')
    !   endif
    !   do ilat=1,nlat
    !   do ilev=1,plev
    !   do ilon=1,nlon
    !     Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    !   end do
    !   end do
    !   end do
    ! endif ! (masterproc) then
    ! call scatter_field_to_chunk(1,Force_nlev,1,Force_nlon,Xtrans,   &
    !                             Target_Q(1,1,begchunk))
 
    ! End Routine
    !------------
    return
   end subroutine ! nncorrector_update

  subroutine init_neural_net()

    implicit none

    integer :: i, k

    allocate(torch_mod (1))
    call torch_mod(1)%load(trim(torch_model), 0) !0 is not using gpu, for now just use cpu for NN inference
    !call torch_mod(1)%load(trim(cb_torch_model), module_use_device) will use gpu if available
    
  end subroutine init_neural_net

  !================================================================
  subroutine corrector_set_profile(rlat,rlon,Force_prof,Wprof,nlev)
   ! 
   ! corrector_SET_PROFILE: for the given lat,lon, and corrector_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on corrector strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Force_prof
   real(r8) rlat,rlon
   real(r8) Wprof(nlev)

   ! Local values
   !----------------
   integer  ilev
   real(r8) Hcoef,latx,lonx,Vmax,Vmin
   real(r8) lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(Force_prof.eq.0) then
     ! No corrector
     !-------------
     Wprof(:)=0.0_r8
   elseif(Force_prof.eq.1) then
     ! Uniform corrector
     !-----------------
     Wprof(:)=1.0_r8
   elseif(Force_prof.eq.2) then
     ! Localized corrector with specified Heaviside window function
     !------------------------------------------------------------
     if(Force_Hwin_max.le.Force_Hwin_min) then
       ! For a constant Horizontal window function, 
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Force_Hwin_lo,Force_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Force_Hwin_lat0
       lonx=rlon-Force_Hwin_lon0
       if(lonx.gt. 180._r8) lonx=lonx-360._r8
       if(lonx.le.-180._r8) lonx=lonx+360._r8

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Force_Hwin_lonWidthH+lonx)/Force_Hwin_lonDelta
       lon_hi=(Force_Hwin_lonWidthH-lonx)/Force_Hwin_lonDelta
       lat_lo=(Force_Hwin_latWidthH+latx)/Force_Hwin_latDelta
       lat_hi=(Force_Hwin_latWidthH-latx)/Force_Hwin_latDelta
       Hcoef=((1._r8+tanh(lon_lo))/2._r8)*((1._r8+tanh(lon_hi))/2._r8) &
            *((1._r8+tanh(lat_lo))/2._r8)*((1._r8+tanh(lat_hi))/2._r8)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Force_Hwin_min)/(Force_Hwin_max-Force_Hwin_min)
       Hcoef=(1._r8-Hcoef)*Force_Hwin_lo + Hcoef*Force_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Force_Vwin_Lindex)/Force_Vwin_Ldelta
       lev_hi=(Force_Vwin_Hindex-float(ilev))/Force_Vwin_Hdelta
       Wprof(ilev)=((1._r8+tanh(lev_lo))/2._r8)*((1._r8+tanh(lev_hi))/2._r8)
     end do 

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if((Vmax.le.Vmin).or.((Force_Vwin_Hindex.ge.(nlev+1)).and. &
                           (Force_Vwin_Lindex.le. 0      )     )) then
       ! For a constant Vertical window function, 
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Force_Vwin_lo,Force_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Force_Vwin_lo + Wprof(:)*(Force_Vwin_hi-Force_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile 
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('corrector_set_profile:: Unknown Force_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! corrector_set_profile
  !================================================================


  !================================================================
  real(r8) function corrector_set_PSprofile(rlat,rlon,Force_PSprof)
   ! 
   ! corrector_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on corrector strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Force_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Force_PSprof.eq.0) then
     ! No corrector
     !-------------
     corrector_set_PSprofile=0.0_r8
   elseif(Force_PSprof.eq.1) then
     ! Uniform corrector
     !-----------------
     corrector_set_PSprofile=1.0_r8
   else
     call endrun('corrector_set_PSprofile:: Unknown Force_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! corrector_set_PSprofile
  !================================================================


end module corrector
