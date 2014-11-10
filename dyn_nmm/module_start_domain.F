module module_start_domain

  use module_domain, only: domain, get_ijk_from_grid
  use module_model_constants, only: r_d, cp, g
  use module_configure, only: model_to_grid_config_rec, &
       grid_config_rec_type, in_use_for_config
  use module_mpp, only: nest, restrt, anal

  implicit none

  private
  public :: start_domain_main

  logical :: allowed_to_read,is_real_nmm

  real, parameter :: capa=r_d/cp, g_inv=1./g, rog=r_d*g_inv
  logical :: S_BDY, N_BDY, W_BDY, E_BDY
  integer :: LM, &
       IDS,IDE,JDS,JDE,KDS,KDE,&
       IMS,IME,JMS,JME,KMS,KME,&
       IPS,IPE,JPS,JPE,KPS,KPE,&
       ITS,ITE,JTS,JTE,KTS,KTE,&
       IFS,IFE,JFS,JFE,KFS,KFE,&
       JTPAD2,JBPAD2,IRPAD2,ILPAD2
  real :: tstart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ITEMP,LOWLYR
    REAL,ALLOCATABLE,DIMENSION(:) :: SFULL,SMID
    REAL,ALLOCATABLE,DIMENSION(:) :: DZS,ZS
    REAL,ALLOCATABLE,DIMENSION(:,:,:) :: RQCBLTEN,RQIBLTEN            &
         &                                    ,RQVBLTEN,RTHBLTEN            &
         &                                    ,RUBLTEN,RVBLTEN              &
         &                                    ,RQCCUTEN,RQICUTEN,RQRCUTEN   &
         &                                    ,RQSCUTEN,RQVCUTEN,RTHCUTEN   &
         &                                    ,RUSHTEN,RVSHTEN              &
         &                                    ,RQCSHTEN,RQISHTEN,RQRSHTEN   &
         &                                    ,RQSSHTEN,RQVSHTEN,RTHSHTEN   &
         &                                    ,RQGSHTEN                     &
         &                                    ,RTHRATEN                     &
         &                                    ,RTHRATENLW,RTHRATENSW
    REAL,ALLOCATABLE,DIMENSION(:,:) :: EMISS,EMTEMP,GLW,HFX           &
         &                                  ,NCA                            &
         &                                  ,QFX,RAINBL,RAINC,RAINNC        &
         &                                  ,RAINNCV                        &
         &                                  ,SNOWNC,SNOWNCV                 &
         &                                  ,GRAUPELNC,GRAUPELNCV           &
         &                                  ,SNOWC,THC,TMN,TSFC

    REAL,ALLOCATABLE,DIMENSION(:,:) :: Z0_DUM, ALBEDO_DUM
    !
    REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ZINT,RRI,CONVFAC,ZMID
    REAL,ALLOCATABLE,DIMENSION(:,:,:) :: T_TRANS,PINT_TRANS
    REAL,ALLOCATABLE,DIMENSION(:,:,:) :: CLDFRA_TRANS
    REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: SCALAR_TRANS
    REAL,ALLOCATABLE,DIMENSION(:,:,:) :: CLDFRA_OLD
    !..Need to fill special height var for setting up initial condition.  G. Thompson
    REAL,ALLOCATABLE,DIMENSION(:,:,:) :: z_at_q

contains

  subroutine sd_dealloc()
#define dealloc_if(X) if(allocated(x)) deallocate(x)
    dealloc_if(ITEMP)
    dealloc_if(LOWLYR)
    dealloc_if(SFULL)
    dealloc_if(SMID)
    dealloc_if(DZS)
    dealloc_if(ZS)
    dealloc_if(RQCBLTEN)
    dealloc_if(RQIBLTEN)
    dealloc_if(RQVBLTEN)
    dealloc_if(RTHBLTEN)
    dealloc_if(RUBLTEN)
    dealloc_if(RVBLTEN)
    dealloc_if(RQCCUTEN)
    dealloc_if(RQICUTEN)
    dealloc_if(RQRCUTEN)
    dealloc_if(RQSCUTEN)
    dealloc_if(RQVCUTEN)
    dealloc_if(RTHCUTEN)
    dealloc_if(RUSHTEN)
    dealloc_if(RVSHTEN)
    dealloc_if(RQCSHTEN)
    dealloc_if(RQISHTEN)
    dealloc_if(RQRSHTEN)
    dealloc_if(RQSSHTEN)
    dealloc_if(RQVSHTEN)
    dealloc_if(RTHSHTEN)
    dealloc_if(RQGSHTEN)
    dealloc_if(RTHRATEN)
    dealloc_if(RTHRATENLW)
    dealloc_if(RTHRATENSW)
    dealloc_if(EMISS)
    dealloc_if(EMTEMP)
    dealloc_if(GLW)
    dealloc_if(HFX)
    dealloc_if(NCA)
    dealloc_if(QFX)
    dealloc_if(RAINBL)
    dealloc_if(RAINC)
    dealloc_if(RAINNC)
    dealloc_if(RAINNCV)
    dealloc_if(SNOWNC)
    dealloc_if(SNOWNCV)
    dealloc_if(GRAUPELNC)
    dealloc_if(GRAUPELNCV)
    dealloc_if(SNOWC)
    dealloc_if(THC)
    dealloc_if(TMN)
    dealloc_if(TSFC)
    dealloc_if(Z0_DUM)
    dealloc_if(ALBEDO_DUM)
    dealloc_if(ZINT)
    dealloc_if(RRI)
    dealloc_if(CONVFAC)
    dealloc_if(ZMID)
    dealloc_if(T_TRANS)
    dealloc_if(PINT_TRANS)
    dealloc_if(CLDFRA_TRANS)
    dealloc_if(SCALAR_TRANS)
    dealloc_if(CLDFRA_OLD)
    dealloc_if(z_at_q)
#undef dealloc_if
  end subroutine sd_dealloc

#if defined( HWRF )
  ! HWRF-specific initialization
  subroutine sd_coupling(gridid,DTC)
    implicit none
    real, intent(out) :: DTC
    integer, intent(in) :: gridid
    type(grid_config_rec_type) :: config_flags
    character*255 :: message

    real :: dt
    integer :: nphs,movemin

    if(grid%id==3) then
       grid%force_sst=1
    else
       grid%force_sst=2
    endif
    grid%dtc=0

    if(model_config_rec%max_dom>=2) then
       ! Normally, in an HWRF simulation, domain 2 has
       ! dt*nphs*movemin = coupling timestep
       ! So if there is a domain 2, we'll use that to
       ! guess the coupling timestep.
       if(gridid==1) then
          ! The first time this routine is called, it is for domain 1,
          ! and domain 2 has not been created yet.  That means domain
          ! 2 still has its dt set to the internal WRF default of 2.
          ! To get domain 2's correct dt, we get domain 1's dt and
          ! divide by 3
          call wrf_message('Guessing coupling timestep from d01''s dt/3 and d02''s nphs and movemin...')
          call model_to_grid_config_rec(1,model_config_rec,config_flags)
          dt=config_flags%dt/3
          call model_to_grid_config_rec(2,model_config_rec,config_flags)
       else
          call wrf_message('Guessing coupling timestep from d02 settings...')
          call model_to_grid_config_rec(2,model_config_rec,config_flags)
          dt=config_flags%dt
       endif
    else
       ! Someone is doing a single domain run, so guess the
       ! coupling timestep using the
       call wrf_message('Guessing coupling timestep from d01 settings because there is no d02...')
       call model_to_grid_config_rec(1,model_config_rec,config_flags)
       dt=config_flags%dt
    endif
    nphs=config_flags%nphs
    movemin=config_flags%movemin

    dtc = dt*nphs*movemin
388 format("dtc=dt*nphs*movemin = ",F0.3,"=",F0.3,"*",I0,"*",I0)
    write(message,388) dtc,dt,nphs,movemin
    call wrf_message(message)
  end subroutine sd_coupling

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_clear_halos(grid,config_flags)
    use module_domain_type, only: domain
    use module_clear_halos, only: clear_ij_halos
    type(domain) :: grid
    type(grid_config_rec_type) :: config_flags
    call clear_ij_halos(grid,config_flags)
  end subroutine sd_clear_halos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_grid_bounds(grid)
    use module_domain_type, only: domain
    type(domain), intent(inout) :: grid
    character*255 :: message
    real :: tempdx, tempdy

    call get_ijk_from_grid(grid,                    &
         IDS,IDE,JDS,JDE,KDS,KDE, &
         IMS,IME,JMS,JME,KMS,KME, &
         IPS,IPE,JPS,JPE,KPS,KPE)

    ! no openmp -> tile=patch
    ITS=IPS
    ITE=IPE
    JTS=JPS
    JTE=JPE
    KTS=KPS
    KTE=KPE

    ! *_BDY: boundary flags.  Are we at this boundary?
    S_BDY=(JPS==JDS)
    N_BDY=(JPE==JDE)
    W_BDY=(IPS==IDS)
    E_BDY=(IPE==IDE)

    ! F = region of tile that is used:
    IFS=ITS
    JFS=JTS
    KFS=KTS
    IFE=min(ITE,IDE-1)
    JFE=min(JTE,JDE-1)
    KFE=min(KTE,KDE-1)

!
    JTPAD2=2
    JBPAD2=2
    IRPAD2=2
    ILPAD2=2
    if(S_BDY)         JBPAD2=0
    if(N_BDY)         JTPAD2=0
    if(W_BDY)         ILPAD2=0
    if(E_BDY)         IRPAD2=0

    LM=KPE-KPS+1

    ! Recalculate grid bounds, regardless of reason for
    ! start_domain_nmm call, or version of NMM:
    has_parent: if(grid%id /= 1) then ! NOTE REQUIREMENT: MOAD IS GRID 1!!
       ! Nest gets dx & dy from parent and grid ratio.
3302   format('Grid ',I0,' calculating west/south bounds relative to parent grid ',I0)
       write(message,3302) grid%id, grid%parents(1)%ptr%id
       call wrf_debug(2,message)
       tempdx=grid%parents(1)%ptr%dx/grid%parent_grid_ratio
       tempdy=grid%parents(1)%ptr%dy/grid%parent_grid_ratio
       grid%wbd0var = grid%parents(1)%ptr%wbd0var + (grid%i_parent_start-1)*2.*grid%parents(1)%ptr%dx + mod(grid%j_parent_start+1,2)*grid%parents(1)%ptr%dx
       grid%sbd0var = grid%parents(1)%ptr%sbd0var + (grid%j_parent_start-1)*grid%parents(1)%ptr%dy
3303   format('Parent wbd0=',F0.3,' sbd0=',F0.3,' i_parent_start=',I0,' j_parent_start=',I0)
       write(message,3303) grid%parents(1)%ptr%wbd0var, &
            grid%parents(1)%ptr%sbd0var, grid%i_parent_start,grid%j_parent_start
       call wrf_debug(2,message)
    else
       ! MOAD gets dx & dy from namelist.
3305   format('Grid ',I0,' calculating west/south bounds as MOAD')
       write(message,3305) grid%id
       call wrf_debug(2,message)
       call nl_get_dx(grid%id,tempdx)
       call nl_get_dy(grid%id,tempdy)
       grid%wbd0var = -(IDE-2)*tempdx
       grid%sbd0var = -((JDE-1)/2)*tempdy
#if ! ( NMM_NEST == 1 )
       ! When NMM is compiled without nesting support, we need to
       ! update the bound metadata here:
       grid%wbd0 = grid%wbd0var
       grid%sbd0 = grid%sbd0var
#endif
    endif has_parent
    if(tempdx<1e-5 .or. tempdy<1e-5) then
       ! Should never reach here unless someone adds a bug to NMM.
       ! The meaning of grid%dx varies during nest initialization
       ! (sometimes it is MOAD dx, sometimes nest dx), so this check
       ! is here just in case someone screws up the code later on.
       write(message,1045) tempdx,tempdy
       call wrf_message('WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING ')
       call wrf_message('Warning: dx or dy are invalid after parent calculation or namelist check.')
       call wrf_message(message)
       call wrf_message('WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING ')
1045   format('Must be >1e-5.  dx=',F0.7,' dy=',F0.7,'.  Grid bounds will be wrong.')
    endif
    write(message,3011) grid%id,tempdx,tempdy,grid%wbd0var,grid%sbd0var
    call wrf_debug(2,message)
3011 format('Grid ',I0,': dx=',F0.3,' dy=',F0.3,' wbd0=',F0.3,' sbd0=',F0.3)
  end subroutine sd_grid_bounds

  subroutine sd_grid_sanity(grid)
    use module_CTLBLK
    use module_domain_type, only: domain
    type(domain), intent(inout) :: grid
    character*255 :: message
    integer :: interval_seconds, restart_interval

    IF(IME>NMM_MAX_DIM )THEN
       WRITE(wrf_err_message,*)                                       &
            'start_domain_nmm ime (',ime,') > ',NMM_MAX_DIM,    &
            '. Increase NMM_MAX_DIM in configure.wrf, clean, and recompile.'
       CALL WRF_ERROR_FATAL(wrf_err_message)
    ENDIF
    !
    IF(JME>NMM_MAX_DIM )THEN
       WRITE(wrf_err_message,*)                                       &
            'start_domain_nmm jme (',jme,') > ',NMM_MAX_DIM,    &
            '. Increase NMM_MAX_DIM in configure.wrf, clean, and recompile.'
       CALL WRF_ERROR_FATAL(wrf_err_message)
    ENDIF

    WRITE(message,196)IHRST,IDAT
    CALL wrf_message(trim(message))
196 FORMAT(' FORECAST BEGINS ',I2,' GMT ',2(I2,'/'),I4)


    !!    Restarts must be made from times for which boundary data is available

    CALL nl_get_interval_seconds(GRID%ID, interval_seconds)
    CALL nl_get_restart_interval(GRID%ID, restart_interval)
    IF (MOD(restart_interval*60,interval_seconds) /= 0) THEN
       WRITE(wrf_err_message,*)' restart_interval is not integer multiplier of interval_seconds'
       CALL WRF_ERROR_FATAL(wrf_err_message)
    END IF

  end subroutine sd_grid_sanity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_set_mpp_bounds(grid)
    use module_mpp
    use module_domain_type, only: domain
    type(domain), intent(inout) :: grid
    character*255 :: message

!!!! For now, set NPES to 1
    NPES=1
!!!!!!tlb
    MY_IS_GLB=IPS
    MY_IE_GLB=IPE-1
    MY_JS_GLB=JPS
    MY_JE_GLB=JPE-1
    !
    IM=IPE-1
    JM=JPE-1
!!!!!!!!!
    !! All "my" variables defined below have had the IDE or JDE specification
    !! reduced by 1
!!!!!!!!!!!

    MYIS=MAX(IDS,IPS)
    MYIE=MIN(IDE-1,IPE)
    MYJS=MAX(JDS,JPS)
    MYJE=MIN(JDE-1,JPE)

    MYIS1  =MAX(IDS+1,IPS)
    MYIE1  =MIN(IDE-2,IPE)
    MYJS2  =MAX(JDS+2,JPS)
    MYJE2  =MIN(JDE-3,JPE)
    !
    MYIS_P1=MAX(IDS,IPS-1)
    MYIE_P1=MIN(IDE-1,IPE+1)
    MYIS_P2=MAX(IDS,IPS-2)
    MYIE_P2=MIN(IDE-1,IPE+2)
    MYIS_P3=MAX(IDS,IPS-3)
    MYIE_P3=MIN(IDE-1,IPE+3)
    MYJS_P3=MAX(JDS,JPS-3)
    MYJE_P3=MIN(JDE-1,JPE+3)
    MYIS_P4=MAX(IDS,IPS-4)
    MYIE_P4=MIN(IDE-1,IPE+4)
    MYJS_P4=MAX(JDS,JPS-4)
    MYJE_P4=MIN(JDE-1,JPE+4)
    MYIS_P5=MAX(IDS,IPS-5)
    MYIE_P5=MIN(IDE-1,IPE+5)
    MYJS_P5=MAX(JDS,JPS-5)
    MYJE_P5=MIN(JDE-1,JPE+5)
    !
    MYIS1_P1=MAX(IDS+1,IPS-1)
    MYIE1_P1=MIN(IDE-2,IPE+1)
    MYIS1_P2=MAX(IDS+1,IPS-2)
    MYIE1_P2=MIN(IDE-2,IPE+2)
    !
    MYJS1_P1=MAX(JDS+1,JPS-1)
    MYJS2_P1=MAX(JDS+2,JPS-1)
    MYJE1_P1=MIN(JDE-2,JPE+1)
    MYJE2_P1=MIN(JDE-3,JPE+1)
    MYJS1_P2=MAX(JDS+1,JPS-2)
    MYJE1_P2=MIN(JDE-2,JPE+2)
    MYJS2_P2=MAX(JDS+2,JPS-2)
    MYJE2_P2=MIN(JDE-3,JPE+2)
    MYJS1_P3=MAX(JDS+1,JPS-3)
    MYJE1_P3=MIN(JDE-2,JPE+3)
    MYJS2_P3=MAX(JDS+2,JPS-3)
    MYJE2_P3=MIN(JDE-3,JPE+3)
  end subroutine sd_set_mpp_bounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_halos(grid,all)
    use module_mpp, only: mype
    use module_dm, only: ntasks_x, ntasks_y, mytask, ntasks, &
         local_communicator, wrf_get_myproc
    use module_domain_type, only: domain
    type(domain), intent(inout) :: grid
    integer :: myproc
    character*255 :: message

#ifdef DM_PARALLEL

    CALL WRF_GET_MYPROC(MYPROC)
    MYPE=MYPROC
#  include <HALO_NMM_INIT_1.inc>
#  include <HALO_NMM_INIT_2.inc>
#  include <HALO_NMM_INIT_3.inc>
#  include <HALO_NMM_INIT_4.inc>
#  include <HALO_NMM_INIT_5.inc>
#  include <HALO_NMM_INIT_6.inc>
#  include <HALO_NMM_INIT_7.inc>
#  include <HALO_NMM_INIT_8.inc>
#  include <HALO_NMM_INIT_9.inc>
#  include <HALO_NMM_INIT_10.inc>
#  include <HALO_NMM_INIT_11.inc>
#  include <HALO_NMM_INIT_12.inc>
#  include <HALO_NMM_INIT_13.inc>
#  include <HALO_NMM_INIT_14.inc>
#  include <HALO_NMM_INIT_15.inc>
    if(all) then
#  include <HALO_NMM_INIT_15B.inc>
    endif
#  include <HALO_NMM_INIT_16.inc>
#  include <HALO_NMM_INIT_17.inc>
#  include <HALO_NMM_INIT_18.inc>
#  include <HALO_NMM_INIT_19.inc>
#  include <HALO_NMM_INIT_20.inc>
#  include <HALO_NMM_INIT_21.inc>
#  include <HALO_NMM_INIT_22.inc>
#  include <HALO_NMM_INIT_23.inc>
#  include <HALO_NMM_INIT_24.inc>
#  include <HALO_NMM_INIT_25.inc>
#  include <HALO_NMM_INIT_26.inc>
#  include <HALO_NMM_INIT_27.inc>
#  include <HALO_NMM_INIT_28.inc>
#  include <HALO_NMM_INIT_29.inc>
#  include <HALO_NMM_INIT_30.inc>
#  include <HALO_NMM_INIT_31.inc>
#  include <HALO_NMM_INIT_32.inc>
#  include <HALO_NMM_INIT_33.inc>
#  include <HALO_NMM_INIT_34.inc>
#  include <HALO_NMM_INIT_35.inc>
#  include <HALO_NMM_INIT_36.inc>
#  include <HALO_NMM_INIT_37.inc>
#  include <HALO_NMM_INIT_38.inc>
#  include <HALO_NMM_INIT_39.inc>
#endif
  end subroutine sd_halos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_refl(grid)
    use module_domain_type, only: domain
    type(domain), intent(inout) :: grid
    integer :: i, k, j

    call wrf_message("Fill REFL_10CM with -35 dBZ")
    if(size(grid%refl_10cm,1)*size(grid%refl_10cm,3)>1) then
       do J=JPS,JPE
          do K=KPS,KPE
             do I=IPS,IPE
                grid%refl_10cm(i,k,j)=-35.0
             enddo
          enddo
          do I=IPS,IPE
             grid%refd_max(i,j)=-35.0
          enddo
       enddo
    endif
  end subroutine sd_refl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_j_const(grid)
    use module_mpp
    use module_domain_type, only: domain
    type(domain), intent(inout) :: grid
    integer :: j

    DO J=MYJS_P4,MYJE_P4
       grid%iheg(J)=MOD(J+1,2)
       grid%ihwg(J)=grid%iheg(J)-1
       grid%iveg(J)=MOD(J,2)
       grid%ivwg(J)=grid%iveg(J)-1
    ENDDO
    !
    DO J=MYJS_P4,MYJE_P4
       grid%ivw(J)=grid%ivwg(J)
       grid%ive(J)=grid%iveg(J)
       grid%ihe(J)=grid%iheg(J)
       grid%ihw(J)=grid%ihwg(J)
    ENDDO
  end subroutine sd_j_const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_random(grid)
    use module_domain_type, only: domain
    type(domain), intent(inout) :: grid
    character*255 :: message
    integer :: one

    randif: IF(in_use_for_config(grid%id,'random'))THEN
       ! Reset random number generator at first initialization,
       ! or after a nest move.  (Need to reset after a nest move
       ! or the leading edge will be filled with 3x3 areas with
       ! the same random number generator state.)
       random_seed=config_flags%random_seed + grid%ntsd

       write(message,'(A,I0,A,I0)') 'Resetting random number for domain ',grid%id,' with seed ',random_seed
       call wrf_message(message)
       one=1
       call srand_grid(grid%randstate1,grid%randstate2, &
            grid%randstate3,grid%randstate4, &
            IDS,IDE,JDS,JDE,one,one, &
            IMS,IME,JMS,JME,one,one, &
            IPS,IPE,JPS,JPE,one,one,random_seed)
       call rand_grid_r4(grid%randstate1,grid%randstate2, &
            grid%randstate3,grid%randstate4, &
            grid%random, &
            IDS,IDE,JDS,JDE,one,one, &
            IMS,IME,JMS,JME,one,one, &
            IPS,IPE,JPS,JPE,one,one)
    else
       grid%random = 0.0
    endif randif
  end subroutine sd_random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef HWRF
  subroutine sd_hifreq(grid,config_flags)
    use module_hifreq, only : hifreq_open

    type(domain) :: grid
    type(grid_config_rec_type) :: config_flags
    if(grid%id==config_flags%high_dom) then
       ! Open HTCF LUN:
       call HIFREQ_OPEN(grid,config_flags)
       ! Open per-nest-move ATCF LUN:
       call HIFREQ_OPEN(grid,config_flags,atcf=.true.)
    elseif(config_flags%high_dom==-99) then
       nestless=.true.
       CALL nl_get_max_dom( 1, max_dom )
       nestdo: do nestid=2,max_dom
          call nl_get_grid_allowed(nestid,grid_allowed)
          if(grid_allowed) then
             call nl_get_parent_id(nestid,parent_id)
             if(parent_id==grid%id) then
                nestless=.false.
                exit nestdo
             endif
          endif
       enddo nestdo
       if(nestless) then
          ! Open HTCF LUN:
          call HIFREQ_OPEN(grid,config_flags)
          ! Open per-nest-move ATCF LUN:
          call HIFREQ_OPEN(grid,config_flags,atcf=.true.)
       endif
    endif
  end subroutine sd_hifreq
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_hires_diag(grid,config_flags)
    type(domain) :: grid
    type(grid_config_rec_type) :: config_flags

#ifdef HWRF
    grid%nomove_freq_hr=config_flags%nomove_freq
    if(allowed_to_read .and. config_flags%high_freq) then
       call sd_hifreq(grid,config_flags)
    end if
#endif

    if(anal .or. (allowed_to_read .and. .not. restrt)) then
       call init_tornado_genesis(grid,config_flags)
    endif

#ifdef HWRF
    if ( .not.is_real_nmm ) then
       call vorttrak_init(grid,config_flags,       &
            (allowed_to_read .and. .not. restrt), &
            IDS,IDE,JDS,JDE,KDS,KDE, &
            IMS,IME,JMS,JME,KMS,KME, &
            ITS,ITE,JTS,JTE,KTS,KTE)
       call init_swath(grid, config_flags, &
            (allowed_to_read .and. .not. restrt) )
    endif
    if(anal .or. allowed_to_read) then
       grid%update_interest=.true.
    endif
#endif
  end subroutine sd_hires_diag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_clean_start(grid)
    integer :: i, j, k
    real :: rg, ht, pm1, apem1, tsfck, psfck, term1, eps, cloges, ese
    DO J=JFS,JFE
       DO I=IFS,IFE
          grid%pdsl(I,J)  =grid%pd(I,J)*grid%res(I,J)
          grid%prec(I,J)  =0.
          IF(allowed_to_read)grid%acprec(I,J)=0.  ! This is gopal's inclusion for moving nest
          grid%cuprec(I,J)=0.
          rg=1./g
          ht=grid%fis(i,j)*rg
!!!       fisx=ht*g
          !          fisx=max(grid%fis(i,j),0.)
          !          prodx=grid%z0(I,J)*Z0MAX
          !          grid%z0(I,J)    =grid%sm(I,J)*Z0SEA+(1.-grid%sm(I,J))*                      &
          !     &                (grid%z0(I,J)*Z0MAX+FISx    *FCM+Z0LAND)
!!!  &                (prodx        +FISx    *FCM+Z0LAND)
          grid%qsh(I,J)   =0.
          grid%akms(I,J)  =0.
          grid%akhs(I,J)  =0.
          grid%twbs(I,J)  =0.
          grid%qwbs(I,J)  =0.
          IF(allowed_to_read)THEN       ! This is gopal's inclusion for moving nest
             grid%cldefi(I,J)=1.
             grid%htop(I,J)  =REAL(KTS)
             grid%htopd(I,J) =REAL(KTS)
             grid%htops(I,J) =REAL(KTS)
             grid%hbot(I,J)  =REAL(KTE)
             grid%hbotd(I,J) =REAL(KTE)
             grid%hbots(I,J) =REAL(KTE)
          ENDIF
          !***
          !***  AT THIS POINT, WE MUST CALCULATE THE INITIAL POTENTIAL TEMPERATURE
          !***  OF THE SURFACE AND OF THE SUBGROUND.
          !***  EXTRAPOLATE DOWN FOR INITIAL SURFACE POTENTIAL TEMPERATURE.
          !***  ALSO DO THE SHELTER PRESSURE.
          !***
          !
          !***  BECAUSE WE REINITIALIZE TOPOGRAPHY, LAND SEA MASK AND FIND THE TEMPERATURE
          !***  FIELD OVER THE NEW TOPOGRAPHY, AFTER THE MOVE, I THINK IT MORE APPROPRIATE
          !***  TO USE grid%nmm_tsk OR grid%sst TO RE-DERIVE grid%ths AND QS (AND CONSEQUENTLY grid%thz0 AND grid%qz0).
          !***  THIS MAY BE MORE CONSISTENT WITH THE PSEUDO-HYDROSTATIC BALANCING THAT IS
          !***  DONE OVER THE NEW TERRAIN (AND WITH NEW grid%sm). gopal!
          !***
          !***

          IF(allowed_to_read)THEN       ! This is gopal's inclusion for moving nest

             PM1=grid%aeta1(KTS)*grid%pdtop+grid%aeta2(KTS)*grid%pdsl(I,J)+grid%pt
             APEM1=(1.E5/PM1)**CAPA

             IF(grid%nmm_tsk(I,J)>=200.)THEN         ! have a specific skin temp, use it
#ifdef HWRF
                grid%ths(I,J)=grid%nmm_tsk(I,J)*(1.+P608*grid%q(I,J,KTS+1))*APEM1
                TSFCK=grid%nmm_tsk(I,J)*(1.+P608*grid%q(I,J,KTS+1))
#else
                grid%ths(I,J)=grid%nmm_tsk(I,J)*APEM1
                TSFCK=grid%nmm_tsk(I,J)
#endif

             ELSE                               ! use lowest layer as a proxy
#ifdef HWRF
                grid%ths(I,J)=grid%t(I,J,KTS)*(1.+P608*grid%q(I,J,KTS+1))*APEM1
                TSFCK=grid%t(I,J,KTS)*(1.+P608*grid%q(I,J,KTS+1))
#else
                grid%ths(I,J)=grid%t(I,J,KTS)*APEM1
                TSFCK=grid%t(I,J,KTS)
#endif
             ENDIF

             PSFCK=grid%pd(I,J)+grid%pdtop+grid%pt
             !
             IF(grid%sm(I,J)<0.5) THEN
                grid%qsh(I,J)=PQ0/PSFCK*EXP(A2*(TSFCK-A3)/(TSFCK-A4))
             ELSEIF(grid%sm(I,J)>0.5) THEN
                grid%ths(I,J)=grid%sst(I,J)*(1.E5/(grid%pd(I,J)+grid%pdtop+grid%pt))**CAPA
             ENDIF
             !
             TERM1=-0.068283/grid%t(I,J,KTS)
             grid%pshltr(I,J)=(grid%pd(I,J)+grid%pdtop+grid%pt)*EXP(TERM1)
             !
             grid%ustar(I,J)=0.1
             grid%thz0(I,J)=grid%ths(I,J)
             grid%qz0(I,J)=grid%qsh(I,J)
             grid%uz0(I,J)=0.
             grid%vz0(I,J)=0.

          ENDIF  ! endif for allowed to read
          ! 
       ENDDO
    ENDDO

    !***
    !***  INITIALIZE CLOUD FIELDS
    !***
    IF (MAXVAL(grid%cwm(ips:ipe,jps:jpe,:)) .gt. 0. .and. MAXVAL(grid%cwm(ips:ipe,jps:jpe,:)) .lt. 1.) then
       CALL wrf_message('appear to have grid%cwm values...do not zero')
    ELSE
       IF(allowed_to_read)THEN       ! This is gopal's inclusion for moving nest
          CALL wrf_message('zeroing grid%cwm')
          DO K=KPS,KPE
             DO J=JFS,JFE
                DO I=IFS,IFE
                   grid%cwm(I,J,K)=0.
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    !***
    !***  INITIALIZE ACCUMULATOR ARRAYS TO ZERO.
    !***
    grid%ardsw=0.0
    grid%ardlw=0.0
    grid%asrfc=0.0
    grid%avrain=0.0
    grid%avcnvc=0.0
    !
    DO J=JFS,JFE
       DO I=IFS,IFE
          grid%acfrcv(I,J)=0.
          grid%ncfrcv(I,J)=0
          grid%acfrst(I,J)=0.
          grid%ncfrst(I,J)=0
          grid%acsnow(I,J)=0.
          grid%acsnom(I,J)=0.
          grid%ssroff(I,J)=0.
          grid%bgroff(I,J)=0.
          grid%alwin(I,J) =0.
          grid%alwout(I,J)=0.
          grid%alwtoa(I,J)=0.
          grid%aswin(I,J) =0.
          grid%aswout(I,J)=0.
          grid%aswtoa(I,J)=0.
          grid%sfcshx(I,J)=0.
          grid%sfclhx(I,J)=0.
          grid%subshx(I,J)=0.
          grid%snopcx(I,J)=0.
          grid%sfcuvx(I,J)=0.
          grid%sfcevp(I,J)=0.
          grid%potevp(I,J)=0.
          grid%potflx(I,J)=0.
       ENDDO
    ENDDO
    !***
    !***  INITIALIZE SATURATION SPECIFIC HUMIDITY OVER THE WATER.
    !***
    EPS=R_D/R_V
    !
    IF(allowed_to_read)THEN       ! This is gopal's inclusion for moving nest
       DO J=JFS,JFE
          DO I=IFS,IFE
             IF(grid%sm(I,J)>0.5)THEN
                CLOGES =-CM1/grid%sst(I,J)-CM2*ALOG10(grid%sst(I,J))+CM3
                ESE    = 10.**(CLOGES+2.)
                grid%qsh(I,J)= grid%sm(I,J)*EPS*ESE/(grid%pd(I,J)+grid%pdtop+grid%pt-ESE*(1.-EPS))
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !***  
    !***  INITIALIZE TURBULENT KINETIC ENERGY (TKE) TO A SMALL
    !***  VALUE (EPSQ2) ABOVE GROUND.  SET TKE TO ZERO IN THE
    !***  THE LOWEST MODEL LAYER.  IN THE LOWEST TWO ATMOSPHERIC
    !***  ETA LAYERS SET TKE TO A SMALL VALUE (Q2INI).
    !***
    !***EROGERS: add check for realistic values of grid%q2
    !
    IF (MAXVAL(grid%q2(ips:ipe,jps:jpe,:)) .gt. epsq2 .and. MAXVAL(grid%q2(ips:ipe,jps:jpe,:)) .lt. 200.) then
       CALL wrf_message('appear to have grid%q2 values...do not zero')
    ELSE
       IF(allowed_to_read)THEN       ! This is gopal's inclusion for moving nest
          CALL wrf_message('zeroing grid%q2')
          DO K=KPS,KPE-1
             DO J=JFS,JFE
                DO I=IFS,IFE
#ifdef HWRF
                   grid%q2(I,J,K)=0.
#else
                   grid%q2(I,J,K)=grid%hbm2(I,J)*EPSQ2
#endif
                ENDDO
             ENDDO
          ENDDO
          !
          DO J=JFS,JFE
             DO I=IFS,IFE
                grid%q2(I,J,LM)    = 0.
#ifdef HWRF
                grid%q2(I,J,KTE-2)= 0.
                grid%q2(I,J,KTE-1)= 0.
#else
                grid%q2(I,J,KTE-2)= grid%hbm2(I,J)*Q2INI
                grid%q2(I,J,KTE-1)= grid%hbm2(I,J)*Q2INI
#endif
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    !***  
    !***  PAD ABOVE GROUND SPECIFIC HUMIDITY IF IT IS TOO SMALL.
    !***  INITIALIZE LATENT HEATING ACCUMULATION ARRAYS.
    !***
    DO K=KPS,KPE
       DO J=JFS,JFE
          DO I=IFS,IFE
             IF(grid%q(I,J,K)<EPSQ)grid%q(I,J,K)=EPSQ
             grid%train(I,J,K)=0.
             grid%tcucn(I,J,K)=0.
          ENDDO
       ENDDO
    ENDDO
    !
    !***
    !***  INITIALIZE MAX/MIN TEMPERATURES.
    !***
    DO J=JFS,JFE
       DO I=IFS,IFE
          grid%tlmax(I,J)=grid%t(I,J,KPS)
          grid%tlmin(I,J)=grid%t(I,J,KPS)
       ENDDO
    ENDDO
    !
    !----------------------------------------------------------------------
    !***  END OF SCRATCH START INITIALIZATION BLOCK.
    !----------------------------------------------------------------------
    !
    CALL wrf_message('INIT:  INITIALIZED ARRAYS FOR CLEAN START')

  end subroutine sd_clean_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_fix_0_T(grid,config_flags)
    type(domain) :: grid
    type(grid_config_rec_type) :: config_flags
    integer :: i, j
       DO J=JFS,JFE
          DO I=IFS,IFE
             !
             IF(grid%t(I,J,KTS)==0.)THEN
                grid%t(I,J,KTS)=grid%t(I,J,KTS+1)
             ENDIF
             !
             TERM1=-0.068283/grid%t(I,J,KTS)
             grid%pshltr(I,J)=(grid%pd(I,J)+grid%pdtop+grid%pt)*EXP(TERM1)
          ENDDO
       ENDDO
    ENDIF

  end subroutine sd_clear_halos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_restart_initializing(grid)
    use module_MPP
    use module_CTLBLK
    implicit none
    integer :: i, j
    TSPH=3600./GRID%DT ! needed?
    grid%nphs0=GRID%NPHS
#ifdef HWRF
    !zhang's doing
    tstart = grid%TSTART  ! NOTE: module-level variable
    !zhang's doing ends
#endif
    
      IF(MYPE==0)THEN
        WRITE( wrf_err_message, * )' start_nmm TSTART=',grid%tstart
        CALL wrf_debug( 1, TRIM(wrf_err_message) )
        WRITE( wrf_err_message, * )' start_nmm TPREC=',grid%tprec
        CALL wrf_debug( 1, TRIM(wrf_err_message) )
        WRITE( wrf_err_message, * )' start_nmm THEAT=',grid%theat
        CALL wrf_debug( 1, TRIM(wrf_err_message) )
        WRITE( wrf_err_message, * )' start_nmm TCLOD=',grid%tclod
        CALL wrf_debug( 1, TRIM(wrf_err_message) )
        WRITE( wrf_err_message, * )' start_nmm TRDSW=',grid%trdsw
        CALL wrf_debug( 1, TRIM(wrf_err_message) )
        WRITE( wrf_err_message, * )' start_nmm TRDLW=',grid%trdlw
        CALL wrf_debug( 1, TRIM(wrf_err_message) )
        WRITE( wrf_err_message, * )' start_nmm TSRFC=',grid%tsrfc
        CALL wrf_debug( 1, TRIM(wrf_err_message) )
        WRITE( wrf_err_message, * )' start_nmm PCPFLG=',grid%pcpflg
        CALL wrf_debug( 1, TRIM(wrf_err_message) )
      ENDIF

      NSTART = INT(grid%TSTART*TSPH+0.5)
!
      grid%ntsd = NSTART

!! want non-zero values for grid%nprec, grid%nheat type vars to avoid problems
!! with mod statements below.

      grid%nprec  = INT(grid%TPREC *TSPH+0.5)
      grid%nheat  = INT(grid%THEAT *TSPH+0.5)
      grid%nclod  = INT(grid%TCLOD *TSPH+0.5)
      grid%nrdsw  = INT(grid%TRDSW *TSPH+0.5)
      grid%nrdlw  = INT(grid%TRDLW *TSPH+0.5)
      grid%nsrfc  = INT(grid%TSRFC *TSPH+0.5)
#ifdef HWRF
!zhang's dong for analysis option:
      grid%NCNVC0  = grid%NCNVC
      grid%NPHS0   = grid%NPHS
#endif
!
!----------------------------------------------------------------------
!
!***  FLAG FOR INITIALIZING ARRAYS, LOOKUP TABLES, & CONSTANTS USED IN
!***  MICROPHYSICS AND RADIATION
!
!----------------------------------------------------------------------
!
      grid%micro_start=.TRUE.
!
!----------------------------------------------------------------------
!***
!***  INITIALIZE ADVECTION TENDENCIES TO ZERO SO THAT
!***  BOUNDARY POINTS WILL ALWAYS BE ZERO
!***
      DO J=JFS,JFE
      DO I=IFS,IFE
        grid%adt(I,J)=0.
        grid%adu(I,J)=0.
        grid%adv(I,J)=0.
      ENDDO
      ENDDO
!----------------------------------------------------------------------
!***
!***  SET INDEX ARRAYS FOR UPSTREAM ADVECTION
!***
!----------------------------------------------------------------------
      DO J=JFS,JFE
        grid%n_iup_h(J)=0
        grid%n_iup_v(J)=0
        grid%n_iup_adh(J)=0
        grid%n_iup_adv(J)=0
!
        DO I=IFS,IFE
          grid%iup_h(I,J)=-999
          grid%iup_v(I,J)=-999
          grid%iup_adh(I,J)=-999
          grid%iup_adv(I,J)=-999
        ENDDO
!
      ENDDO

  end subroutine sd_restart_initializing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef NO_UPSTREAM_ADVECTION
  subroutine sd_upstream_advection(grid)
    use module_mpp
    implicit none
    INTEGER :: ILPAD2,IRPAD2,JBPAD2,JTPAD2, I,J,K, II,JJ,KK
    INTEGER :: KNTI, ISTART,IEND, JAM, KNT
!
!***  n_iup_h HOLDS THE NUMBER OF MASS POINTS NEEDED IN EACH ROW
!***  FOR UPSTREAM ADVECTION (FULL ROWS IN THE 3RD THROUGH 7TH
!***  ROWS FROM THE SOUTH AND NORTH GLOBAL BOUNDARIES AND 
!***  FOUR POINTS ADJACENT TO THE WEST AND EAST GLOBAL BOUNDARIES
!***  ON ALL OTHER INTERNAL ROWS).  SIMILARLY FOR n_iup_v.
!***  BECAUSE OF HORIZONTAL OPERATIONS, THESE POINTS EXTEND OUTSIDE 
!***  OF THE UPSTREAM REGION SOMEWHAT.
!***  n_iup_adh HOLDS THE NUMBER OF MASS POINTS NEEDED IN EACH ROW
!***  FOR THE COMPUTATION OF THE TENDENCIES THEMSELVES (adt, ADQ2M
!***  AND ADQ2L); SPECIFICALLY THESE TENDENCIES ARE ONLY DONE IN
!***  THE UPSTREAM REGION.
!***  n_iup_adv HOLDS THE NUMBER OF MASS POINTS NEEDED IN EACH ROW
!***  FOR THE VELOCITY POINT TENDENCIES.
!***  iup_h AND iup_v HOLD THE ACTUAL I VALUES USED IN EACH ROW.
!***  LIKEWISE FOR iup_adh AND iup_adv. 
!***  ALSO, SET upstrm FOR THOSE TASKS AROUND THE GLOBAL EDGE.
!
      grid%upstrm=.FALSE.
!
      IF(S_BDY)THEN
        grid%upstrm=.TRUE.

!
        DO JJ=1,7
          J=JJ      ! -MY_JS_GLB+1
          KNTI=0
          DO I=MYIS_P2,MYIE_P2
            grid%iup_h(IMS+KNTI,J)=I
            grid%iup_v(IMS+KNTI,J)=I
            KNTI=KNTI+1
          ENDDO
          grid%n_iup_h(J)=KNTI
          grid%n_iup_v(J)=KNTI
        ENDDO
!
        DO JJ=3,5
          J=JJ      ! -MY_JS_GLB+1
          KNTI=0
          ISTART=MYIS1_P2
          IEND=MYIE1_P2
          IF(E_BDY)IEND=IEND-MOD(JJ+1,2)
          DO I=ISTART,IEND
            grid%iup_adh(IMS+KNTI,J)=I
            KNTI=KNTI+1
          ENDDO
          grid%n_iup_adh(J)=KNTI
!
          KNTI=0
          ISTART=MYIS1_P2
          IEND=MYIE1_P2
          IF(E_BDY)IEND=IEND-MOD(JJ,2)
          DO I=ISTART,IEND
            grid%iup_adv(IMS+KNTI,J)=I
            KNTI=KNTI+1
          ENDDO
          grid%n_iup_adv(J)=KNTI
        ENDDO
      ENDIF
!
      IF(N_BDY)THEN
        grid%upstrm=.TRUE.
!
        DO JJ=JDE-7, JDE-1 ! JM-6,JM
          J=JJ      ! -MY_JS_GLB+1
          KNTI=0
          DO I=MYIS_P2,MYIE_P2
            grid%iup_h(IMS+KNTI,J)=I
            grid%iup_v(IMS+KNTI,J)=I
            KNTI=KNTI+1
          ENDDO
          grid%n_iup_h(J)=KNTI
          grid%n_iup_v(J)=KNTI
        ENDDO
!
        DO JJ=JDE-5, JDE-3 ! JM-4,JM-2
          J=JJ      ! -MY_JS_GLB+1
          KNTI=0
          ISTART=MYIS1_P2
          IEND=MYIE1_P2
          IF(E_BDY)IEND=IEND-MOD(JJ+1,2)
          DO I=ISTART,IEND
            grid%iup_adh(IMS+KNTI,J)=I
            KNTI=KNTI+1
          ENDDO
          grid%n_iup_adh(J)=KNTI
!
          KNTI=0
          ISTART=MYIS1_P2
          IEND=MYIE1_P2
          IF(E_BDY)IEND=IEND-MOD(JJ,2)
          DO I=ISTART,IEND
            grid%iup_adv(IMS+KNTI,J)=I
            KNTI=KNTI+1
          ENDDO
          grid%n_iup_adv(J)=KNTI
        ENDDO
      ENDIF
!
      IF(W_BDY)THEN
        grid%upstrm=.TRUE.
        DO JJ=8,JDE-8   ! JM-7
          IF(JJ>=MY_JS_GLB-2.AND.JJ<=MY_JE_GLB+2)THEN
            J=JJ      ! -MY_JS_GLB+1
!
            DO I=1,4
              grid%iup_h(IMS+I-1,J)=I
              grid%iup_v(IMS+I-1,J)=I
            ENDDO
            grid%n_iup_h(J)=4
            grid%n_iup_v(J)=4
          ENDIF
        ENDDO
!
        DO JJ=6,JDE-6   ! JM-5
          IF(JJ>=MY_JS_GLB-2.AND.JJ<=MY_JE_GLB+2)THEN
            J=JJ      ! -MY_JS_GLB+1
            KNTI=0
            IEND=2+MOD(JJ,2)
            DO I=2,IEND
              grid%iup_adh(IMS+KNTI,J)=I
              KNTI=KNTI+1
            ENDDO
            grid%n_iup_adh(J)=KNTI
!
            KNTI=0
            IEND=2+MOD(JJ+1,2)
            DO I=2,IEND
              grid%iup_adv(IMS+KNTI,J)=I
              KNTI=KNTI+1
            ENDDO
            grid%n_iup_adv(J)=KNTI
!
          ENDIF
        ENDDO
      ENDIF
!
      CALL WRF_GET_NPROCX(INPES)
!
      IF(E_BDY)THEN
        grid%upstrm=.TRUE.
        DO JJ=8,JDE-8   ! JM-7
          IF(JJ>=MY_JS_GLB-2.AND.JJ<=MY_JE_GLB+2)THEN
            J=JJ      ! -MY_JS_GLB+1
            IEND=IM-MOD(JJ+1,2)
            ISTART=IEND-3
!
!***  IN CASE THERE IS ONLY A SINGLE GLOBAL TASK IN THE
!***  I DIRECTION THEN WE MUST ADD THE WESTSIDE UPSTREAM
!***  POINTS TO THE EASTSIDE POINTS IN EACH ROW.
!
            KNTI=0
            IF(INPES.EQ.1)KNTI=grid%n_iup_h(J)
!
            DO II=ISTART,IEND
              I=II      ! -MY_IS_GLB+1
              grid%iup_h(IMS+KNTI,J)=I
              KNTI=KNTI+1
            ENDDO
            grid%n_iup_h(J)=KNTI
          ENDIF
        ENDDO
!
        DO JJ=6,JDE-6   ! JM-5
          IF(JJ>=MY_JS_GLB-2.AND.JJ<=MY_JE_GLB+2)THEN
            J=JJ      ! -MY_JS_GLB+1
            IEND=IM-1-MOD(JJ+1,2)
            ISTART=IEND-MOD(JJ,2)
            KNTI=0
            IF(INPES==1)KNTI=grid%n_iup_adh(J)
            DO II=ISTART,IEND
              I=II      ! -MY_IS_GLB+1
              grid%iup_adh(IMS+KNTI,J)=I
              KNTI=KNTI+1
            ENDDO
            grid%n_iup_adh(J)=KNTI
          ENDIF
        ENDDO
!***
        DO JJ=8,JDE-8  ! JM-7
          IF(JJ>=MY_JS_GLB-2.AND.JJ<=MY_JE_GLB+2)THEN
            J=JJ      ! -MY_JS_GLB+1
            IEND=IM-MOD(JJ,2)
            ISTART=IEND-3
            KNTI=0
            IF(INPES==1)KNTI=grid%n_iup_v(J)
!
            DO II=ISTART,IEND
              I=II      ! -MY_IS_GLB+1
              grid%iup_v(IMS+KNTI,J)=I
              KNTI=KNTI+1
            ENDDO
            grid%n_iup_v(J)=KNTI
          ENDIF
        ENDDO
!
        DO JJ=6,JDE-6  !  JM-5
          IF(JJ>=MY_JS_GLB-2.AND.JJ<=MY_JE_GLB+2)THEN
            J=JJ      ! -MY_JS_GLB+1
            IEND=IM-1-MOD(JJ,2)
            ISTART=IEND-MOD(JJ+1,2)
            KNTI=0
            IF(INPES==1)KNTI=grid%n_iup_adv(J)
            DO II=ISTART,IEND
              I=II      ! -MY_IS_GLB+1
              grid%iup_adv(IMS+KNTI,J)=I
              KNTI=KNTI+1
            ENDDO
            grid%n_iup_adv(J)=KNTI
          ENDIF
        ENDDO
      ENDIF
!----------------------------------------------------------------------
      jam=6+2*(JDE-JDS-1-9)
!
!***  EXTRACT em AND emt FOR THE LOCAL SUBDOMAINS
!
      DO J=MYJS_P5,MYJE_P5
        grid%em_loc(J)=-9.E9
        grid%emt_loc(J)=-9.E9
      ENDDO
!!!   IF(IBROW==1)THEN
      IF(S_BDY)THEN
        DO J=3,5
          grid%em_loc(J)=grid%em(J-2)
          grid%emt_loc(J)=grid%emt(J-2)
        ENDDO
      ENDIF
!!!   IF(ITROW==1)THEN
      IF(N_BDY)THEN
        KNT=3
        DO JJ=JDE-5,JDE-3 ! JM-4,JM-2
          KNT=KNT+1
          J=JJ      ! -MY_JS_GLB+1
          grid%em_loc(J)=grid%em(KNT)
          grid%emt_loc(J)=grid%emt(KNT)
        ENDDO
      ENDIF
!!!   IF(ILCOL==1)THEN
      IF(W_BDY)THEN
        KNT=6
        DO JJ=6,JDE-6 ! JM-5
          KNT=KNT+1
          IF(JJ>=MY_JS_GLB-2.AND.JJ<=MY_JE_GLB+2)THEN
            J=JJ      ! -MY_JS_GLB+1
            grid%em_loc(J)=grid%em(KNT)
            grid%emt_loc(J)=grid%emt(KNT)
          ENDIF
        ENDDO
      ENDIF
!!!   IF(IRCOL==1)THEN
      IF(E_BDY)THEN
        KNT=6+JDE-11 ! JM-10
        DO JJ=6,JDE-6 ! JM-5
          KNT=KNT+1
          IF(JJ>=MY_JS_GLB-2.AND.JJ<=MY_JE_GLB+2)THEN
            J=JJ      ! -MY_JS_GLB+1
            grid%em_loc(J)=grid%em(KNT)
            grid%emt_loc(J)=grid%emt(KNT)
          ENDIF
        ENDDO
      ENDIF
  end subroutine sd_upstream_advection
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_zero_diag_arrays(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    integer :: i,j,ns
         GRID%NSOIL= GRID%NUM_SOIL_LAYERS
        DO J=JFS,JFE
        DO I=IFS,IFE
          grid%pctsno(I,J)=-999.0
          IF(grid%sm(I,J)<0.5)THEN
              grid%cmc(I,J)=0.0
!              grid%cmc(I,J)=grid%canwat(i,j)   ! tgs
            IF(grid%sice(I,J)>0.5)THEN
!***
!***  SEA-ICE CASE
!***
              grid%smstav(I,J)=1.0
              grid%smstot(I,J)=1.0
              grid%ssroff(I,J)=0.0
              grid%bgroff(I,J)=0.0
              grid%cmc(I,J)=0.0
              DO NS=1,GRID%NSOIL
                grid%smc(I,NS,J)=1.0
!               grid%sh2o(I,NS,J)=0.05
                grid%sh2o(I,NS,J)=1.0
              ENDDO
            ENDIF
          ELSE
!***
!***  WATER CASE
!***
            grid%smstav(I,J)=1.0
            grid%smstot(I,J)=1.0
            grid%ssroff(I,J)=0.0
            grid%bgroff(I,J)=0.0
            grid%soiltb(I,J)=273.16
            grid%grnflx(I,J)=0.
            grid%subshx(I,J)=0.0
            grid%acsnow(I,J)=0.0
            grid%acsnom(I,J)=0.0
            grid%snopcx(I,J)=0.0
            grid%cmc(I,J)=0.0
            grid%sno(I,J)=0.0
            DO NS=1,GRID%NSOIL
              grid%smc(I,NS,J)=1.0
              grid%stc(I,NS,J)=273.16
!             grid%sh2o(I,NS,J)=0.05
              grid%sh2o(I,NS,J)=1.0
            ENDDO
          ENDIF
!
        ENDDO
        ENDDO
!
        grid%aphtim=0.0
        grid%aratim=0.0
        grid%acutim=0.0
!
  end subroutine sd_zero_diag_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_init_radtn(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    integer,dimension(3) :: LPTOP
    integer :: i,j,k
    integer :: K400,KCCO2
    real,allocatable,dimension(:) :: PHALF
    real :: PSUM, SLPM, PDIF, PSS, TIME, ADDL
    !
    !----------------------------------------------------------------------
    !***  INITIALIZE RADTN VARIABLES
    !***  CALCULATE THE NUMBER OF STEPS AT EACH POINT.
    !***  THE ARRAY 'lvl' WILL COORDINATE VERTICAL LOCATIONS BETWEEN
    !***  THE LIFTED WORKING ARRAYS AND THE FUNDAMENTAL MODEL ARRAYS.
    !***  lvl HOLDS THE HEIGHT (IN MODEL LAYERS) OF THE TOPOGRAPHY AT
    !***  EACH GRID POINT.
    !----------------------------------------------------------------------
    !   
    DO J=JFS,JFE
       DO I=IFS,IFE
          grid%lvl(I,J)=LM-KTE
       ENDDO
    ENDDO
    !***
    !***  DETERMINE MODEL LAYER LIMITS FOR HIGH(3), MIDDLE(2),
    !***  AND LOW(1) CLOUDS.  ALSO FIND MODEL LAYER THAT IS JUST BELOW
    !***  (HEIGHT-WISE) 400 MB. (K400)
    !*** 
    K400=0
    PSUM=grid%pt
    SLPM=101325.
    PDIF=SLPM-grid%pt
    DO K=1,LM
       PSUM=PSUM+grid%deta(K)*PDIF
       IF(LPTOP(3)==0)THEN
          IF(PSUM>PHITP)LPTOP(3)=K
       ELSEIF(LPTOP(2)==0)THEN
          IF(PSUM>PMDHI)LPTOP(2)=K
       ELSEIF(K400==0)THEN
          IF(PSUM>P400)K400=K
       ELSEIF(LPTOP(1)==0)THEN
          IF(PSUM>PLOMD)LPTOP(1)=K
       ENDIF
    ENDDO
    !***
    !*** CALL GRADFS ONCE TO CALC. CONSTANTS AND GET O3 DATA
    !***
    KCCO2=0
    !***
    !*** CALCULATE THE MIDLAYER PRESSURES IN THE STANDARD ATMOSPHERE
    !***
    PSS=101325.
    PDIF=PSS-grid%pt
    !
    ALLOCATE(PHALF(LM+1),STAT=I)
    !
    DO K=KPS,KPE-1
       PHALF(K+1)=grid%aeta(K)*PDIF+grid%pt
    ENDDO

    !
    PHALF(1)=0.
    PHALF(LM+1)=PSS
    !***
!!!   CALL GRADFS(PHALF,KCCO2,NUNIT_CO2)
    !***
    !***  CALL SOLARD TO CALCULATE NON-DIMENSIONAL SUN-EARTH DISTANCE
    !***
!!!   IF(MYPE==0)CALL SOLARD(SUN_DIST)
!!!   CALL MPI_BCAST(SUN_DIST,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)

    !***
    !***  CALL ZENITH SIMPLY TO GET THE DAY OF THE YEAR FOR
    !***  THE SETUP OF THE OZONE DATA
    !***
    TIME=(grid%ntsd-1)*GRID%DT
    !
!!!   CALL ZENITH(TIME,DAYI,HOUR)
    !
    ADDL=0.
    IF(MOD(IDAT(3),4)==0)ADDL=1.
    !
!!!   CALL O3CLIM
    !
    !
    DEALLOCATE(PHALF)

  end subroutine sd_init_radtn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_init_turbl(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    integer :: i, j
    real :: ulm,vlm,tlm,plm,qlm,apelm,term1,apelmnw,thlm,dplm,dzlm,fac1,fac2
    logical :: do_me
    real :: zoqing

       DO J=JFS,JFE
          DO I=IFS,IFE
             !***
             !***  TRY A SIMPLE LINEAR INTERP TO GET 2/10 M VALUES
             !***
#ifdef HWRF
             !zhang's doing
             IF(.NOT.RESTRT .OR. .NOT.allowed_to_read) then
                grid%PDSL(I,J)=grid%PD(I,J)*grid%RES(I,J)
             endif
             !end of zhang's doing
#else
             grid%PDSL(I,J)=grid%PD(I,J)*grid%RES(I,J)
#endif
             !
             ULM=grid%u(I,J,KTS)
             VLM=grid%v(I,J,KTS)
             TLM=grid%t(I,J,KTS)
             QLM=grid%q(I,J,KTS)
             PLM=grid%aeta1(KTS)*grid%pdtop+grid%aeta2(KTS)*grid%pdsl(I,J)+grid%pt
             APELM=(1.0E5/PLM)**CAPA
             TERM1=-0.068283/grid%t(I,J,KTS)
             grid%pshltr(I,J)=(grid%pd(I,J)+grid%pdtop+grid%pt)*EXP(TERM1)
             APELMNW=(1.0E5/grid%pshltr(I,J))**CAPA
             THLM=TLM*APELM
             DPLM=(grid%deta1(KTS)*grid%pdtop+grid%deta2(KTS)*grid%pdsl(I,J))*0.5
             DZLM=R_D*DPLM*TLM*(1.+P608*QLM)/(G*PLM)
             FAC1=10./DZLM
             FAC2=(DZLM-10.)/DZLM
             IF(DZLM<=10.)THEN
                FAC1=1.
                FAC2=0.
             ENDIF
             !
#ifdef HWRF
             !zhang's doing
             do_me=(.NOT.RESTRT .OR. .NOT.allowed_to_read)
             !end of zhang's doing
#else
             do_me=(.NOT.RESTRT)
#endif
             if(do_me) then
                grid%th10(I,J)=FAC2*grid%ths(I,J)+FAC1*THLM
                grid%q10(I,J)=FAC2*grid%qsh(I,J)+FAC1*QLM
#ifdef HWRF
                IF(grid%sm(I,J).LT.0.5)THEN
                   grid%u10(I,J)=ULM*(log(10./grid%z0(I,J))/log(DZLM/grid%z0(I,J)))      ! this is all Qingfu's doing
                   grid%v10(I,J)=VLM*(log(10./grid%z0(I,J))/log(DZLM/grid%z0(I,J)))
                   ZOQING=1.944*SQRT(grid%u10(I,J)*grid%u10(I,J)+grid%v10(I,J)*grid%v10(I,J))
                   IF(ZOQING.GT.60.)THEN
                      grid%u10(I,J)=grid%u10(I,J)*(1.12-7.2/ZOQING)
                      grid%v10(I,J)=grid%v10(I,J)*(1.12-7.2/ZOQING)
                   ENDIF
                ELSE
                   ZOQING=(0.074*SQRT(ULM*ULM+VLM*VLM)-0.58)*1.0e-3
                   ZOQING=MAX(ZOQING,grid%z0(I,J))          ! for winds greater than 12.5 m/s
                   grid%u10(I,J)=ULM*(log(10./ZOQING))/log(DZLM/ZOQING)      ! this is all Qingfu's doing
                   grid%v10(I,J)=VLM*(log(10./ZOQING))/log(DZLM/ZOQING)
                   ZOQING=1.944*SQRT(grid%u10(I,J)*grid%u10(I,J)+grid%v10(I,J)*grid%v10(I,J))
                   IF(ZOQING.GT.60.)THEN
                      grid%u10(I,J)=grid%u10(I,J)*(1.12-7.2/ZOQING)
                      grid%v10(I,J)=grid%v10(I,J)*(1.12-7.2/ZOQING)
                   END IF
                ENDIF
#else
                grid%u10(I,J)=ULM
                grid%v10(I,J)=VLM
#endif
             ENDIF
             !
             !        FAC1=2./DZLM
             !        FAC2=(DZLM-2.)/DZLM
             !        IF(DZLM.LE.2.)THEN
             !          FAC1=1.
             !          FAC2=0.
             !        ENDIF
             !
             IF(.NOT.RESTRT.OR.NEST)THEN
                !
                IF ( (THLM-grid%ths(I,J))>2.0) THEN  ! weight differently in different scenarios
                   FAC1=0.3
                   FAC2=0.7
                ELSE
                   FAC1=0.8
                   FAC2=0.2
                ENDIF

#ifdef HWRF
                grid%tshltr(I,J)=0.2*grid%ths(I,J)+0.8*THLM
                grid%qshltr(I,J)=0.2*grid%qsh(I,J)+0.8*QLM
#else
                grid%tshltr(I,J)=FAC2*grid%ths(I,J)+FAC1*THLM
                grid%qshltr(I,J)=FAC2*grid%qsh(I,J)+FAC1*QLM
#endif
             ENDIF
             !***
             !***  NEED TO CONVERT TO THETA IF IS THE RESTART CASE
             !***  AS CHKOUT.f WILL CONVERT TO TEMPERATURE
             !***
             !EROGERS: COMMENT OUT IN WRF-NMM
             !***
             !       IF(RESTRT)THEN
             !         grid%tshltr(I,J)=grid%tshltr(I,J)*APELMNW
             !       ENDIF
          ENDDO
       ENDDO
  end subroutine sd_init_turbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_init_old(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    integer :: i,j,k

    DO K=KPS,KPE
       DO J=JFS,JFE
          DO I=ifs,ife
             grid%told(I,J,K)=grid%t(I,J,K)   ! grid%t AT TAU-1
             grid%uold(I,J,K)=grid%u(I,J,K)   ! grid%u AT TAU-1
             grid%vold(I,J,K)=grid%v(I,J,K)   ! grid%v AT TAU-1
          ENDDO
       ENDDO
    ENDDO

  end subroutine sd_init_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_reset_nonhydro(grid,nest)
    use module_domain_type, only: domain
    implicit none
    logical, intent(in) :: nest
    type(domain), intent(inout) :: grid
    logical :: do_me
    real :: apelp, rtopp, dz, sfct, a
    integer :: i,j,k
    !----------------------------------------------------------------------
    !***  INITIALIZE NONHYDROSTATIC QUANTITIES
    !----------------------------------------------------------------------
    !
!!!!	SHOULD grid%dwdt BE REDEFINED IF RESTRT?

    IF((.NOT.RESTRT.OR.NEST).AND. allowed_to_read)THEN ! This is gopal's inclusion for moving nest
       DO K=KPS,KPE
          DO J=JFS,JFE
             DO I=IFS,IFE
                grid%dwdt(I,J,K)=1.
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !***
#ifdef HWRF
    IF(.NOT.RESTRT .OR. .NOT.allowed_to_read) THEN !zhang's doing
#endif
       IF(GRID%SIGMA==1)THEN
          DO J=JFS,JFE
             DO I=IFS,IFE
                grid%pdsl(I,J)=grid%pd(I,J)
             ENDDO
          ENDDO
       ELSE
          DO J=JFS,JFE
             DO I=IFS,IFE
                grid%pdsl(I,J)=grid%res(I,J)*grid%pd(I,J)
             ENDDO
          ENDDO
       ENDIF
#ifdef HWRF
    ENDIF !zhang's doing
#endif
    !
    !***
    !
    !
!!!!	SHOULD pint,z,w BE REDEFINED IF RESTRT?

    WRITE( wrf_err_message, * )' restrt=',restrt,' nest=',nest
    CALL wrf_debug( 0, TRIM(wrf_err_message) )
    WRITE( wrf_err_message, * )' grid%pdtop=',grid%pdtop,' grid%pt=',grid%pt
    CALL wrf_debug( 0, TRIM(wrf_err_message) )
#ifdef HWRF
    !zhang's doing
    do_me=(.NOT.RESTRT.OR.NEST .OR. .NOT.allowed_to_read)
    !end of zhang's doing
#else
    do_me=(.NOT.RESTRT.OR.NEST)
#endif
    if(do_me) then
       DO K=KPS,KPE
          DO J=JFS,JFE
             DO I=IFS,IFE
                grid%pint(I,J,K)=grid%eta1(K)*grid%pdtop+grid%eta2(K)*grid%pdsl(I,J)+grid%pt
                grid%z(I,J,K)=grid%pint(I,J,K)
                grid%w(I,J,K)=0.
#ifdef IDEAL_NMM_TC
                grid%f(I,J)=0.5*GRID%DT*3.15656e-5 ! IDEAL CASE 0.5*DT*f (and not f!)
#endif
             ENDDO
          ENDDO
       ENDDO
    ENDIF
#ifdef HWRF
    !zhang's doing
    IF(.NOT.RESTRT.OR.NEST .OR. .NOT.allowed_to_read)THEN
#endif

       DO K=KTS,KTE-1
          DO J=JFS,JFE
             DO I=IFS,IFE
                grid%rtop(I,J,K)=(grid%q(I,J,K)*P608-grid%cwm(I,J,K)+1.)*grid%t(I,J,K)*R_D/ &
                     ((grid%pint(I,J,K+1)+grid%pint(I,J,K))*0.5)
             ENDDO
          ENDDO
       ENDDO
#ifdef HWRF
    ENDIF    !zhang 
#endif

#ifdef HWRF
    hwrfx_mslp: if(grid%vortex_tracker /= 1) then
       DO J=JFS,JFE
          DO I=IFS,IFE
             grid%Z(I,J,KFS)=grid%FIS(I,J)*GI
          ENDDO
       ENDDO

       ! Z now correctly calculated after nest move.  Needed by membrane MSLP.
       DO K=KFS,KFE
          DO J=JFS,JFE
             DO I=IFS,IFE
                APELP      = (grid%PINT(I,J,K+1)+grid%PINT(I,J,K))
                RTOPP      = TRG*grid%T(I,J,K)*(1.0+grid%Q(I,J,K)*P608)/APELP
                DZ         = RTOPP*(grid%DETA1(K)*grid%PDTOP+grid%DETA2(K)*grid%PD(I,J))
                grid%Z(I,J,K+1) = grid%Z(I,J,K) + DZ
             ENDDO
          ENDDO
       ENDDO

       DO K=KFS,KFE
          DO J=JFS,JFE
             DO I=IFS,IFE
                grid%Z(i,j,k)=(grid%Z(i,j,k)+grid%Z(i,j,k+1))*0.5
             ENDDO
          ENDDO
       ENDDO

       grid%MSLP=-9999.99
       DO J=JFS,JFE
          DO I=IFS,IFE
             SFCT      = grid%T(I,J,1)*(1.+D608*grid%Q(I,J,1)) + LAPSR*grid%Z(I,J,1)
             A         = LAPSR*grid%FIS(i,j)*gi/SFCT
             grid%MSLP(I,J) = grid%PINT(I,J,1)*(1-A)**COEF2
          ENDDO
       ENDDO

    endif hwrfx_mslp
#endif

  end subroutine sd_reset_nonhydro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_restrict_accel(grid)
    use module_domain_type, only: domain
    use module_MPP
    implicit none
    type(domain), intent(inout) :: grid
    integer :: i,j,jhh,jx,ihl,ix,ihh
    REAL, parameter :: EPSB=0.1,EPSIN=9.8
    !----------------------------------------------------------------------
    !***  RESTRICTING THE ACCELERATION ALONG THE BOUNDARIES
    !----------------------------------------------------------------------
    !
    DO J=JFS,JFE
       DO I=IFS,IFE
          grid%dwdtmn(I,J)=-EPSIN
          grid%dwdtmx(I,J)= EPSIN
       ENDDO
    ENDDO
    !
    !***
    IF(JHL>1)THEN
       JHH=JDE-1-JHL+1 ! JM-JHL+1
       IHL=JHL/2+1
       !
       DO J=1,JHL
          IF(J>=MY_JS_GLB-JBPAD2.AND.J<=MY_JE_GLB+JTPAD2)THEN
             JX=J      ! -MY_JS_GLB+1
             DO I=1,IDE-1 ! IM
                IF(I>=MY_IS_GLB-ILPAD2.AND.I<=MY_IE_GLB+IRPAD2)THEN
                   IX=I      ! -MY_IS_GLB+1
                   grid%dwdtmn(IX,JX)=-EPSB
                   grid%dwdtmx(IX,JX)= EPSB
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       !
       DO J=JHH,JDE-1   ! JM
          IF(J>=MY_JS_GLB-JBPAD2.AND.J<=MY_JE_GLB+JTPAD2)THEN
             JX=J      ! -MY_JS_GLB+1
             DO I=1,IDE-1 ! IM
                IF(I>=MY_IS_GLB-ILPAD2.AND.I<=MY_IE_GLB+IRPAD2)THEN
                   IX=I      ! -MY_IS_GLB+1
                   grid%dwdtmn(IX,JX)=-EPSB
                   grid%dwdtmx(IX,JX)= EPSB
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       !
       DO J=1,JDE-1 ! JM
          IF(J>=MY_JS_GLB-JBPAD2.AND.J<=MY_JE_GLB+JTPAD2)THEN
             JX=J      ! -MY_JS_GLB+1
             DO I=1,IHL
                IF(I>=MY_IS_GLB-ILPAD2.AND.I<=MY_IE_GLB+IRPAD2)THEN
                   IX=I      ! -MY_IS_GLB+1
                   grid%dwdtmn(IX,JX)=-EPSB
                   grid%dwdtmx(IX,JX)= EPSB
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       !
       DO J=1,JDE-1 ! JM
          IF(J>=MY_JS_GLB-JBPAD2.AND.J<=MY_JE_GLB+JTPAD2)THEN
             JX=J      ! -MY_JS_GLB+1
             ! moved this line to inside the J-loop, 20030429, jm
             IHH=IDE-1-IHL+MOD(J,2) ! IM-IHL+MOD(J,2)
             DO I=IHH,IDE-1 ! IM
                IF(I>=MY_IS_GLB-ILPAD2.AND.I<=MY_IE_GLB+IRPAD2)THEN
                   IX=I      ! -MY_IS_GLB+1
                   grid%dwdtmn(IX,JX)=-EPSB
                   grid%dwdtmx(IX,JX)= EPSB
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       !
    ENDIF


  end subroutine sd_restrict_accel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_physics_init(grid)
    use module_domain_type, only: domain
    use module_ext_internal
    implicit none
    type(grid_config_rec_type) :: config_flags
    type(domain), intent(inout) :: grid

    TYPE(WRFU_Time) :: currentTime
    integer :: i, j, k, n
    real :: sfczpsurf, exnsfc, thsij, plyr, tl, cwml, pdtot, rpdtot, pdbot
    real :: dsigsum, dsig
    character*255 :: message
    integer :: num_ozmixm, num_aerosolc
    logical :: start_of_simulation
    INTEGER :: hr, mn, sec, ms, rc

    ALLOCATE(SFULL(KMS:KME))           ; SFULL    = 0.
    ALLOCATE(SMID(KMS:KME))            ; SMID     = 0.
    ALLOCATE(EMISS(IMS:IME,JMS:JME))   ; EMISS    = 0.
    ALLOCATE(EMTEMP(IMS:IME,JMS:JME))  ; EMTEMP   = 0.
    ALLOCATE(GLW(IMS:IME,JMS:JME))     ; GLW      = 0.
    ALLOCATE(HFX(IMS:IME,JMS:JME))     ; HFX      = 0.
    ALLOCATE(LOWLYR(IMS:IME,JMS:JME))  ; LOWLYR   = 0.
    !     ALLOCATE(grid%mavail(IMS:IME,JMS:JME))  ; grid%mavail   = 0.
    ALLOCATE(NCA(IMS:IME,JMS:JME))     ; NCA      = 0.
    ALLOCATE(QFX(IMS:IME,JMS:JME))     ; QFX      = 0.
    ALLOCATE(RAINBL(IMS:IME,JMS:JME))  ; RAINBL   = 0.
    ALLOCATE(RAINC(IMS:IME,JMS:JME))   ; RAINC    = 0.
    ALLOCATE(RAINNC(IMS:IME,JMS:JME))  ; RAINNC   = 0.
    ALLOCATE(RAINNCV(IMS:IME,JMS:JME)) ; RAINNCV  = 0.
    ALLOCATE(SNOWNC(IMS:IME,JMS:JME))  ; SNOWNC   = 0.
    ALLOCATE(SNOWNCV(IMS:IME,JMS:JME)) ; SNOWNCV  = 0.
    ALLOCATE(GRAUPELNC(IMS:IME,JMS:JME))  ; GRAUPELNC   = 0.
    ALLOCATE(GRAUPELNCV(IMS:IME,JMS:JME)) ; GRAUPELNCV  = 0.

    ALLOCATE(ZS(KMS:KME))              ; ZS       = 0.
    ALLOCATE(SNOWC(IMS:IME,JMS:JME))   ; SNOWC    = 0.
    ALLOCATE(THC(IMS:IME,JMS:JME))     ; THC      = 0.
    ALLOCATE(TMN(IMS:IME,JMS:JME))     ; TMN      = 0.
    ALLOCATE(TSFC(IMS:IME,JMS:JME))    ; TSFC     = 0.
    ALLOCATE(Z0_DUM(IMS:IME,JMS:JME))  ; Z0_DUM   = 0.
    ALLOCATE(ALBEDO_DUM(IMS:IME,JMS:JME))  ; ALBEDO_DUM   = 0.

    ALLOCATE(DZS(KMS:KME))                         ; DZS = 0.
    ALLOCATE(RQCBLTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQCBLTEN = 0.
    ALLOCATE(RQIBLTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQIBLTEN = 0.
    ALLOCATE(RQVBLTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQVBLTEN =  0.
    ALLOCATE(RTHBLTEN(IMS:IME,KMS:KME,JMS:JME))    ; RTHBLTEN =  0.
    ALLOCATE(RUBLTEN(IMS:IME,KMS:KME,JMS:JME))     ; RUBLTEN = 0.
    ALLOCATE(RVBLTEN(IMS:IME,KMS:KME,JMS:JME))     ; RVBLTEN = 0.
    ALLOCATE(RQCCUTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQCCUTEN = 0.
    ALLOCATE(RQICUTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQICUTEN  = 0.
    ALLOCATE(RQRCUTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQRCUTEN = 0.
    ALLOCATE(RQSCUTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQSCUTEN = 0.
    ALLOCATE(RQVCUTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQVCUTEN = 0.
    ALLOCATE(RTHCUTEN(IMS:IME,KMS:KME,JMS:JME))    ; RTHCUTEN = 0.
    ALLOCATE(RUSHTEN(IMS:IME,KMS:KME,JMS:JME))     ; RUSHTEN = 0.
    ALLOCATE(RVSHTEN(IMS:IME,KMS:KME,JMS:JME))     ; RVSHTEN = 0.
    ALLOCATE(RQCSHTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQCSHTEN = 0.
    ALLOCATE(RQISHTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQISHTEN  = 0.
    ALLOCATE(RQRSHTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQRSHTEN = 0.
    ALLOCATE(RQSSHTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQSSHTEN = 0.
    ALLOCATE(RQGSHTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQGSHTEN = 0.
    ALLOCATE(RQVSHTEN(IMS:IME,KMS:KME,JMS:JME))    ; RQVSHTEN = 0.
    ALLOCATE(RTHSHTEN(IMS:IME,KMS:KME,JMS:JME))    ; RTHSHTEN = 0.
    ALLOCATE(RTHRATEN(IMS:IME,KMS:KME,JMS:JME))    ; RTHRATEN  = 0.
    ALLOCATE(RTHRATENLW(IMS:IME,KMS:KME,JMS:JME))  ; RTHRATENLW = 0.
    ALLOCATE(RTHRATENSW(IMS:IME,KMS:KME,JMS:JME))  ; RTHRATENSW = 0.
    ALLOCATE(ZINT(IMS:IME,KMS:KME,JMS:JME))        ; ZINT = 0.
    ALLOCATE(CONVFAC(IMS:IME,KMS:KME,JMS:JME))     ; CONVFAC = 0.
    ALLOCATE(PINT_TRANS(IMS:IME,KMS:KME,JMS:JME))  ; PINT_TRANS = 0.
    ALLOCATE(T_TRANS(IMS:IME,KMS:KME,JMS:JME))     ;  T_TRANS = 0.
    ALLOCATE(RRI(IMS:IME,KMS:KME,JMS:JME))         ;  RRI = 0.
    ALLOCATE(CLDFRA_TRANS(IMS:IME,KMS:KME,JMS:JME)); CLDFRA_TRANS = 0.
    ALLOCATE(SCALAR_TRANS(IMS:IME,KMS:KME,JMS:JME,NUM_SCALAR))
    ALLOCATE(CLDFRA_OLD(IMS:IME,KMS:KME,JMS:JME))  ; CLDFRA_OLD = 0.
    ALLOCATE(Z_AT_Q(IMS:IME,KMS:KME,JMS:JME))      ; z_at_q = 0.

    GRID%RADT=GRID%NRADS*GRID%DT/60.
    GRID%BLDT=GRID%NPHS*GRID%DT/60.
    GRID%CUDT=GRID%NCNVC*GRID%DT/60.
    GRID%GSMDT=GRID%NPHS*GRID%DT/60.
    !
    ! translate scalar(i,j,k,n) to scalar_trans(i,k,j,n)
    DO N=1,NUM_SCALAR
       !$omp parallel do                                                       &
       !$omp& private(i,j,k)
       DO K=KMS,KME
          DO J=JMS,JME
             DO I=IMS,IME
                SCALAR_TRANS(I,K,J,N)=SCALAR(I,J,K,N)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    DO J=MYJS,MYJE
       DO I=MYIS,MYIE
          SFCZ=grid%fis(I,J)*G_INV
          ZINT(I,KTS,J)=SFCZ
#ifdef HWRF
          !zhang's doing
          IF(.NOT.RESTRT .OR. .NOT.allowed_to_read) then
             grid%PDSL(I,J)=grid%PD(I,J)*grid%RES(I,J)
          endif
          !end of zhang's doing
#else
          grid%pdsl(I,J)=grid%pd(I,J)*grid%res(I,J)
#endif
          PSURF=grid%pint(I,J,KTS)
          EXNSFC=(1.E5/PSURF)**CAPA
          grid%xland(I,J)=grid%sm(I,J)+1.
          THSIJ=(grid%sst(I,J)*EXNSFC)*(grid%xland(I,J)-1.)                         &
               &        +grid%ths(I,J)*(2.-grid%sm(I,J))
          TSFC(I,J)=THSIJ/EXNSFC
          !
          DO K=KTS,KTE-1
             PLYR=(grid%pint(I,J,K)+grid%pint(I,J,K+1))*0.5
             TL=grid%t(I,J,K)
             CWML=grid%cwm(I,J,K)
             RRI(I,K,J)=R_D*TL*(1.+P608*grid%q(I,J,K))/PLYR
             ZINT(I,K+1,J)=ZINT(I,K,J)+TL/PLYR                             & 
                  *(grid%deta1(K)*grid%pdtop+grid%deta2(K)*grid%pdsl(I,J))*ROG        & 
                  *(grid%q(I,J,K)*P608-CWML+1.)
          ENDDO
          !
          !        DO K=KTS,KTE
!!!       ZMID(I,K,J)=0.5*(ZINT(I,K,J)+ZINT(I,K+1,J))
          !        ENDDO
       ENDDO
    ENDDO
    !
    !-----------------------------------------------------------------------
    !***  RECREATE SIGMA VALUES AT LAYER INTERFACES FOR THE FULL VERTICAL
    !***  DOMAIN FROM THICKNESS VALUES FOR THE TWO SUBDOMAINS.
    !***  NOTE: KTE=NUMBER OF LAYERS PLUS ONE
    !-----------------------------------------------------------------------
    !
    PDTOT=101325.-grid%pt
    RPDTOT=1./PDTOT
    PDBOT=PDTOT-grid%pdtop
    SFULL(KTS)=1.
    SFULL(KTE)=0.
    DSIGSUM = 0.
    DO K=KTS+1,KTE
       DSIG=(grid%deta1(K-1)*grid%pdtop+grid%deta2(K-1)*PDBOT)*RPDTOT
       DSIGSUM=DSIGSUM+DSIG
       SFULL(K)=SFULL(K-1)-DSIG
       SMID(K-1)=0.5*(SFULL(K-1)+SFULL(K))
    ENDDO
    DSIG=(grid%deta1(KTE-1)*grid%pdtop+grid%deta2(KTE-1)*PDBOT)*RPDTOT
    DSIGSUM=DSIGSUM+DSIG
    SMID(KTE-1)=0.5*(SFULL(KTE-1)+SFULL(KTE))
#ifdef HWRF
    !zhang's doing
    if(.NOT.RESTRT .OR. .NOT.allowed_to_read)grid%LU_INDEX=grid%IVGTYP
    !end of zhang's doing
#else
    grid%lu_index=grid%ivgtyp
#endif
    IF(.NOT.RESTRT)THEN
       DO J=MYJS,MYJE
          DO I=MYIS,MYIE
             Z0_DUM(I,J)=grid%z0(I,J) ! hold
             ALBEDO_DUM(I,J)=grid%albedo(I,J) ! Save albedos
          ENDDO
       ENDDO
    ENDIF
    !
    !***  Always define the quantity grid%z0base

4041 format('Bounds: ip=',I0,',',I0,' jp=',I0,',',I0,' myi=',I0,',',I0,&
         ' myj=',I0,',',I0)
    write(message,4041) ips,ipe,jps,jpe,myis,myie,myjs,myje
    call wrf_message(message)
    IF(.NOT.RESTRT)THEN
       DO J=MYJS,MYJE
          DO I=MYIS,MYIE
             !
             IF(grid%sm(I,J)==0)then
                grid%z0base(I,J)=VZ0TBL_24(grid%ivgtyp(I,J))+Z0LAND
             ELSE
                grid%z0base(I,J)=VZ0TBL_24(grid%ivgtyp(I,J))+Z0SEA
             ENDIF
             !
          ENDDO
       ENDDO
    ENDIF
    !
    ! when allocating CAM radiation 4d arrays (ozmixm, aerosolc) these are not needed
    num_ozmixm=1
    num_aerosolc=1

    ! Several arguments are RCONFIG entries in Registry.NMM. Registry no longer
    ! includes these as dummy arguments or declares them.  Access them from 
    ! GRID.  JM 20050819
#ifndef WRF_NMM_NEST
    ! NOTE: we always get here for all NMM configurations because the
    !       #if condition is wrong.  Leaving this as is, just in case
    !       there was a good reason for it.
    grid%moved = .FALSE.
#endif

    IF (GRID%RESTART) THEN
       LRESTART = GRID%RESTART
    ELSE
       IF (grid%moved) THEN
          LRESTART = .TRUE.
       ELSE
          LRESTART = .FALSE.
       ENDIF
    END IF

    CALL PHY_INIT(GRID%ID,CONFIG_FLAGS,GRID%DT,LRESTART,SFULL,SMID    &
         &             ,grid%pt,TSFC,GRID%RADT,GRID%BLDT,GRID%CUDT,GRID%GSMDT    &
         &             ,grid%DUCUDT, grid%DVCUDT                            &
         &             ,RTHCUTEN, RQVCUTEN, RQRCUTEN                        &
         &             ,RQCCUTEN, RQSCUTEN, RQICUTEN                        &
         &             ,RUSHTEN,  RVSHTEN,  RTHSHTEN                        &
         &             ,RQVSHTEN, RQRSHTEN, RQCSHTEN                        &
         &             ,RQSSHTEN, RQISHTEN, RQGSHTEN                        &
         &             ,RUBLTEN,RVBLTEN,RTHBLTEN                            &
         &             ,RQVBLTEN,RQCBLTEN,RQIBLTEN                          &
         &             ,RTHRATEN,RTHRATENLW,RTHRATENSW                      &
         &             ,STEPBL,STEPRA,STEPCU                                &
         &             ,grid%w0avg, RAINNC, RAINC, grid%raincv, RAINNCV               &
         &             ,SNOWNC, SNOWNCV, GRAUPELNC, GRAUPELNCV              &
         &             ,z_at_q, grid%qnwfa2d                                &
         &             ,scalar_trans(ims,kms,jms,1),num_scalar              &
         &             ,grid%re_cloud, grid%re_ice, grid%re_snow            & ! G. Thompson
         &             ,grid%has_reqc,grid%has_reqi,grid%has_reqs           & ! G. Thompson
         &             ,NCA,GRID%SWRAD_SCAT                                 &
         &             ,grid%cldefi,LOWLYR                                       &
         &             ,grid%mass_flux                                           &
         &             ,grid%rthften, grid%rqvften                                    &
         &             ,CLDFRA_TRANS,CLDFRA_OLD,GLW,grid%gsw,EMISS,EMTEMP,grid%lu_index&
         &             ,GRID%LANDUSE_ISICE, GRID%LANDUSE_LUCATS             &
         &             ,GRID%LANDUSE_LUSEAS, GRID%LANDUSE_ISN               &
         &             ,GRID%LU_STATE                                       &
         &             ,grid%xlat,grid%xlong,grid%glat,grid%glon&
         &             ,grid%albedo,grid%albbck                             &
         &             ,GRID%GMT,GRID%JULYR,GRID%JULDAY                     &
         &             ,GRID%LEVSIZ, NUM_OZMIXM, NUM_AEROSOLC, GRID%PAERLEV &
         &             ,grid%alevsiz, grid%no_src_types                     &
         &             ,TMN,grid%xland,grid%znt,grid%z0,grid%ustar,grid%mol,grid%pblh,grid%tke_pbl             &
         &             ,grid%exch_h,THC,SNOWC,grid%mavail,HFX,QFX,RAINBL              &
         &             ,grid%stc,grid%sldpth,grid%DZSoil,GRID%NUM_SOIL_LAYERS,WARM_RAIN           &
         &             ,ADV_MOIST_COND,IS_CAMMGMP_USED                      &
         &             ,grid%apr_gr,grid%apr_w,grid%apr_mc,grid%apr_st,grid%apr_as                   &
         &             ,grid%apr_capma,grid%apr_capme,grid%apr_capmi                       &
         &             ,grid%xice,grid%xice,grid%vegfra,grid%snow,grid%canwat,grid%smstav                 &
         &             ,grid%smstot, grid%sfcrunoff,grid%udrunoff,grid%grdflx,grid%acsnow            &
         &             ,grid%acsnom,grid%ivgtyp,grid%isltyp,grid%sfcevp,grid%smc                     &
         &             ,grid%sh2o, grid%snowh, grid%smfr3d                                 &  ! temporary
         &             ,grid%SNOALB                                         &
         &             ,GRID%DX,GRID%DY,grid%f_ice_phy,grid%f_rain_phy,grid%f_rimef_phy    &
         &             ,grid%mp_restart_state,grid%tbpvs_state,grid%tbpvs0_state           &
         &             ,ALLOWED_TO_READ,grid%moved,START_OF_SIMULATION                    &
         &             ,1                                                   & ! lagday
         &             ,IDS, IDE, JDS, JDE, KDS, KDE                        &
         &             ,IMS, IME, JMS, JME, KMS, KME                        &
         &             ,ITS, ITE, JTS, JTE, KTS, KTE                        &
         &             ,NUM_URBAN_LAYERS,NUM_URBAN_HI                       &
         &             ,GRID%RAINCV_A,GRID%RAINCV_B                         &
         &               ,ISNOWXY=grid%ISNOWXY, ZSNSOXY=grid%ZSNSOXY, TSNOXY=grid%TSNOXY,    & ! Optional Noah-MP
         &                SNICEXY=grid%SNICEXY, SNLIQXY=grid%SNLIQXY, TVXY=grid%TVXY,        & ! Optional Noah-MP
         &                TGXY=grid%TGXY, CANICEXY=grid%CANICEXY,                            & ! Optional Noah-MP
         &                CANLIQXY=grid%CANLIQXY, EAHXY=grid%EAHXY,                          & ! Optional Noah-MP
         &                TAHXY=grid%TAHXY, CMXY=grid%CMXY,                                  & ! Optional Noah-MP
         &                CHXY=grid%CHXY, FWETXY=grid%FWETXY, SNEQVOXY=grid%SNEQVOXY,        & ! Optional Noah-MP
         &                ALBOLDXY=grid%ALBOLDXY, QSNOWXY=grid%QSNOWXY,                      & ! Optional Noah-MP
         &                WSLAKEXY=grid%WSLAKEXY, ZWTXY=grid%ZWTXY, WAXY=grid%WAXY,          & ! Optional Noah-MP
         &                WTXY=grid%WTXY, LFMASSXY=grid%LFMASSXY, RTMASSXY=grid%RTMASSXY,    & ! Optional Noah-MP
         &                STMASSXY=grid%STMASSXY, WOODXY=grid%WOODXY,                        & ! Optional Noah-MP
         &                STBLCPXY=grid%STBLCPXY, FASTCPXY=grid%FASTCPXY,                    & ! Optional Noah-MP
         &                XSAIXY=grid%XSAIXY,                                                & ! Optional Noah-MP
         &                T2MVXY=grid%T2MVXY, T2MBXY=grid%T2MBXY, CHSTARXY=grid%CHSTARXY,    & ! Optional Noah-MP
         &                smoiseq=grid%smoiseq, smcwtdxy=grid%smcwtdxy, rechxy=grid%rechxy,   &
         &                deeprechxy=grid%deeprechxy,                                         &
         &                lakedepth2d=grid%lakedepth2d,  savedtke12d=grid%savedtke12d,  snowdp2d=grid%snowdp2d,   h2osno2d=grid%h2osno2d,               & !lake
         &                snl2d= grid%snl2d, t_grnd2d=grid%t_grnd2d, t_lake3d=grid%t_lake3d, lake_icefrac3d=grid%lake_icefrac3d,                        & !lake
         &                z_lake3d=grid%z_lake3d, dz_lake3d=grid%dz_lake3d, t_soisno3d=grid%t_soisno3d, h2osoi_ice3d=grid%h2osoi_ice3d,                 & !lake
         &                h2osoi_liq3d=grid%h2osoi_liq3d, h2osoi_vol3d=grid%h2osoi_vol3d, z3d=grid%z3d, dz3d=grid%dz3d,                                 & !lake
         &                zi3d=grid%zi3d, watsat3d=grid%watsat3d, csol3d=grid%csol3d, tkmg3d=grid%tkmg3d,                                               & !lake
         &                tkdry3d=grid%tkdry3d, tksatu3d=grid%tksatu3d, lake2d=grid%lake2d,                                                             & !lake
         &                lakedepth_default=config_flags%lakedepth_default, lake_min_elev=config_flags%lake_min_elev, lake_depth=grid%lake_depth,               & !lake
         &                lake_depth_flag=grid%LAKE_DEPTH_FLAG, use_lakedepth=grid%use_lakedepth,                                                      & !lake
         &                sf_surface_mosaic=config_flags%sf_surface_mosaic, mosaic_cat=config_flags%mosaic_cat, nlcat=1,     & ! Noah tiling
         &                MAXPATCH=1,ccn_conc=config_flags%ccn_conc & ! CLM
         &                )


    !
    !
    !mp replace F*_PHY with values defined in module_initialize_real.F?
#ifdef HWRF
    IF (.NOT. RESTRT .and. ALLOWED_TO_READ) THEN   !zhang
       moist = 0.0
       if(size(grid%f_ice)>1) grid%f_ice = grid%f_ice_phy
       if(size(grid%f_rimef)>1) grid%f_rimef = grid%f_rimef_phy
       if(size(grid%f_rain)>1) grid%f_rain = grid%f_rain_phy
    ENDIF                  !zhang
#endif

  end subroutine sd_physics_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_copy_or_reset_phys(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    
    !-----------------------------------------------------------------------
    !
#ifdef HWRF
    do_me=(NSTART.EQ.0 .or. .not.allowed_to_read )
#else
    do_me=(NSTART==0)
#endif
    if(do_me) then

       DO J=JMS,JME
          DO I=IMS,IME
             grid%z0(I,J)=grid%z0base(I,J)
          ENDDO
       ENDDO

       DO K=KMS,KME
          DO J=JMS,JME
             DO I=IMS,IME
                grid%cldfra(I,J,K)=CLDFRA_TRANS(I,K,J)
             ENDDO
          ENDDO
       ENDDO

       !
       !
       !mp replace F*_PHY with values defined in module_initialize_real.F?
#ifdef HWRF
       IF (.NOT. RESTRT .and. ALLOWED_TO_READ) THEN   !zhang
          moist = 0.0
          if(size(grid%f_ice)>1) grid%f_ice = grid%f_ice_phy
          if(size(grid%f_rimef)>1) grid%f_rimef = grid%f_rimef_phy
          if(size(grid%f_rain)>1) grid%f_rain = grid%f_rain_phy
       ENDIF                  !zhang
#endif

    ENDIF

  end subroutine sd_copy_or_reset_phys
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_calc_micro(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    integer :: i,j,k,i_m
    LOGICAL :: ETAMP_Regional, ICE1_indx, ICE2_indx
    ! Also sum cloud water, ice, rain, snow, graupel into Ferrier
    ! cwm array (if any hydrometeors found and non-zero from
    ! initialization package).  Then, determine fractions ice and
    ! rain from species.
    have_valid_cwm: IF ( &
         .not. (MAXVAL(grid%cwm(ips:ipe,jps:jpe,:)).gt.0. &
         .and. MAXVAL(grid%cwm(ips:ipe,jps:jpe,:)).lt.1.) ) then    
       do i_m = 2, num_moist
          if (i_m.ne.p_qv) &
               &       CALL wrf_message(' summing moist(:,:,:,i_m) into cwm array')
          DO K=KPS,KPE
             DO J=JFS,JFE
                DO I=IFS,IFE
                   IF ( (moist(I,J,K,i_m).gt.EPSQ) .and. (i_m.ne.p_qv) ) THEN  
                      grid%cwm(I,J,K) = grid%cwm(I,J,K) + moist(I,J,K,i_m)               
                   ENDIF
                enddo
             enddo
          enddo
       enddo

       IF (size(grid%f_ice)>1 .and. size(grid%f_rain)>1) then
          IF( .not. ( (maxval(grid%f_ice(ips:ipe,:,jps:jpe)) &
               +maxval(grid%f_rain(ips:ipe,:,jps:jpe))) .gt. EPSQ) ) THEN
             ETAMP_Regional=.FALSE.    !-- Regional NAM or HRW (Ferrier) microphysics
             if (model_config_rec%mp_physics(grid%id).EQ.ETAMPOLD .OR.          &
                  &          model_config_rec%mp_physics(grid%id).EQ.ETAMPNEW )             &
                  &          ETAMP_Regional=.TRUE.
             CALL wrf_message(' computing grid%f_ice')
             do i_m = 2, num_moist
                ICE1_indx=.FALSE.
                IF (i_m==P_qi .or. i_m==P_qg ) ICE1_indx=.TRUE.
                ICE2_indx=ICE1_indx
                IF (i_m==P_qs) ICE2_indx=.TRUE.
                IF (ETAMP_Regional .AND. ICE1_indx) THEN
                   DO J=JFS,JFE
                      DO K=KPS,KPE
                         DO I=IFS,IFE
                            moist(I,J,K,p_qs)=moist(I,J,K,p_qs)+moist(I,J,K,i_m)
                            moist(I,J,K,i_m) =0.
                         enddo
                      enddo
                   enddo
                   if(size(grid%f_ice)>1) then
                      DO J=JFS,JFE
                         DO K=KPS,KPE
                            DO I=IFS,IFE
                               grid%f_ice(I,K,J) = grid%f_ice(I,K,J) + moist(I,J,K,i_m)
                            enddo
                         enddo
                      enddo
                   endif
                ENDIF
             enddo
             CALL wrf_message(' computing f_rain')
             !
             if(size(grid%f_ice)>1 .and. size(grid%f_rain)>1) then
                DO J=JFS,JFE
                   DO K=KPS,KPE
                      DO I=IFS,IFE
                         IF(grid%f_ice(i,k,j)<=EPSQ)THEN
                            grid%f_ice(I,K,J)=0.
                         ELSE
                            grid%f_ice(I,K,J) = grid%f_ice(I,K,J)/grid%cwm(I,J,K)
                         ENDIF
                         IF ( (moist(I,J,K,p_qr)+moist(I,J,K,p_qc)).gt.EPSQ) THEN
                            IF(moist(i,j,k,p_qr)<=EPSQ)THEN
                               grid%f_rain(I,K,J)=0.
                            ELSE
                               grid%f_rain(I,K,J) = moist(i,j,k,p_qr) &
                                    &                    / (moist(i,j,k,p_qr)+moist(i,j,k,p_qc))
                            ENDIF
                         ENDIF
                      enddo
                   enddo
                enddo
             endif
          ENDIF
       ENDIF
    ENDIF have_valid_cwm
    ! End addition by Greg Thompson

    copy_f_ice: if(size(grid%f_ice)>1) then
       IF (maxval(grid%f_ice(ips:ipe,:,jps:jpe)) .gt. 0.) THEN
          do J=JMS,JME
             do K=KMS,KME
                do I=IMS,IME
                   grid%f_ice_phy(I,K,J)=grid%f_ice(I,K,J)
                enddo
             enddo
          enddo
       ENDIF
    endif copy_f_ice

    copy_f_rain: if(size(grid%f_rain)>1) then
       IF (maxval(grid%f_rain(ips:ipe,:,jps:jpe)) .gt. 0.) THEN
          do J=JMS,JME
             do K=KMS,KME
                do I=IMS,IME
                   grid%f_rain_phy(I,K,J)=grid%f_rain(I,K,J)
                enddo
             enddo
          enddo
       ENDIF
    endif copy_f_rain

    copy_f_rimef: if(size(grid%f_rimef)>1) then
       IF (maxval(grid%f_rimef(ips:ipe,:,jps:jpe)) .gt. 0.) THEN
          do J=JMS,JME
             do K=KMS,KME
                do I=IMS,IME
                   grid%f_rimef_phy(I,K,J)=grid%f_rimef(I,K,J)
                enddo
             enddo
          enddo
       ENDIF
    endif copy_f_rimef

  end subroutine sd_calc_micro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_replace_data(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    IF (.NOT. RESTRT) THEN
       !-- Replace albedos if original albedos are nonzero
       IF(MAXVAL(ALBEDO_DUM(ips:ipe,jps:jpe))>0.)THEN
          DO J=JMS,JME
             DO I=IMS,IME
                grid%albedo(I,J)=ALBEDO_DUM(I,J)
             ENDDO
          ENDDO
       ENDIF
    ENDIF

  end subroutine sd_replace_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_reset_precip(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
#ifdef HWRF
    do_me=(.NOT. RESTRT .OR. .NOT.allowed_to_read) 
    !zhang's doing
#else
    do_me=(.NOT.RESTRT)
#endif
    if(do_me) then
       DO J=jps,min(jpe,jde-1)
          DO I=ips,min(ipe,ide-1)
             grid%aprec(I,J)=RAINNC(I,J)*1.E-3
             grid%cuprec(I,J)=grid%raincv(I,J)*1.E-3
          ENDDO
       ENDDO
    endif
  end subroutine sd_zero_diag_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_compute_ct_ct2(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    !
    ! computes ct and ct2 for topo_wind=1 or topo_wind=2
    ! NOTE: not yet implemented for NMM: =1.0 will have no effect
    !
    grid%ctopo=1.
    grid%ctopo2=1.
    !
  end subroutine sd_compute_ct_ct2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd_init_gwd(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    integer :: i,j
    real :: dtphs,cen_lat,cen_lon
    !-----------------------------------------------------------------------
    !---- Initialization for gravity wave drag (GWD) & mountain blocking (MB)
    !
#if (HWRF == 1)
    do_me=(grid%gwd_opt .eq. 2 .AND. grid%id .eq. 1 .AND. allowed_to_read) 
#else
    do_me=(grid%gwd_opt .eq. 2 )
#endif
    if(do_me) then
       CALL nl_get_cen_lat(GRID%ID, CEN_LAT)    !-- CEN_LAT in deg
       CALL nl_get_cen_lon(GRID%ID, CEN_LON)    !-- CEN_LON in deg
       DTPHS=grid%dt*grid%nphs
       CALL GWD_init(DTPHS,GRID%DX,GRID%DY,CEN_LAT,CEN_LON,RESTRT        &
            &              ,grid%glat,grid%glon,grid%crot,grid%srot,grid%hangl                          &
            &              ,IDS,IDE,JDS,JDE,KDS,KDE                            &
            &              ,IMS,IME,JMS,JME,KMS,KME                            &
            &              ,ITS,ITE,JTS,JTE,KTS,KTE )
    ENDIF
    IF(.NOT.RESTRT)THEN
       DO J=MYJS,MYJE
          DO I=MYIS,MYIE
             grid%ugwdsfc(I,J)=0.
             grid%vgwdsfc(I,J)=0.
          ENDDO
       ENDDO
    ENDIF
  end subroutine sd_init_gwd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_set_time(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid

    ! Set GMT, JULDAY, and JULYR outside of phy_init because it is no
    ! longer called inside phy_init due to moving nest changes.  (When
    ! nests move phy_init may not be called on a process if, for
    ! example, it is a moving nest and if this part of the domain is
    ! not being initialized (not the leading edge).)  Calling
    ! domain_setgmtetc() here will avoid this problem when NMM moves
    ! to moving nests.
    CALL domain_setgmtetc( GRID, START_OF_SIMULATION )

    if(restrt) then
#ifdef HWRF
       ! FIXME: REMOVE THIS - THE CODE IN THE "ELSE" BLOCK IS CORRECT 
       !zhang 
       CALL nl_get_julyr (grid%id, grid%julyr)
       CALL nl_get_julday (grid%id, grid%julday)
       CALL nl_get_gmt (grid%id, grid%gmt)
       !zhang end
#else
       ! FIXME: THIS SHOULD BE USED FOR HWRF TOO
       CALL domain_clock_get( grid, current_time=currentTime )
       CALL WRFU_TimeGet( currentTime, YY=grid%julyr, dayOfYear=grid%julday, &
            H=hr, M=mn, S=sec, MS=ms, rc=rc)
       grid%gmt=hr+real(mn)/60.+real(sec)/3600.+real(ms)/(1000*3600)
       WRITE( wrf_err_message , * ) 'DEBUG start_domain_nmm():  gmt = ',grid%gmt
       CALL wrf_debug( 150, TRIM(wrf_err_message) )
#endif
    endif

  end subroutine sd_set_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sd_set_flags_after_phys_init(grid)
    use module_domain_type, only: domain
    implicit none
    type(domain), intent(inout) :: grid
    ! Several arguments are RCONFIG entries in Registry.NMM. Registry no longer
    ! includes these as dummy arguments or declares them.  Access them from 
    ! GRID.  JM 20050819
#ifndef WRF_NMM_NEST
    ! NOTE: we always get here for all NMM configurations because the
    !       #if condition is wrong.  Leaving this as is, just in case
    !       there was a good reason for it.
    grid%moved = .FALSE.
#endif

    IF (GRID%RESTART) THEN
       LRESTART = GRID%RESTART
    ELSE
       IF (grid%moved) THEN
          LRESTART = .TRUE.
       ELSE
          LRESTART = .FALSE.
       ENDIF
    END IF
    
  end subroutine sd_set_flags_after_phys_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine start_domain_main( &
       grid, allowed_to_read_ &
                                !
#include <dummy_new_args.inc>
       !
    )
    use module_mpp
    use module_configure, only: model_to_grid_config_rec, grid_config_rec_type
#ifdef HWRF
    USE module_tornado_genesis, only: init_tornado_genesis
    use module_stats_for_move, only: vorttrak_init
    use module_swath, only: init_swath
#endif
    type(domain), intent(inout) :: grid
    logical, intent(in) :: allowed_to_read

    ! -- locals --
    type(grid_config_rec_type) :: config_flags
    logical :: doit

    call start_timing
    CALL model_to_grid_config_rec(grid%id,model_config_rec,config_flags)

    restrt=config_flags%restart
    anal=config_flags%analysis                
    allowed_to_read=allowed_to_read_
    is_real_nmm=program_name(1:8) == 'REAL_NMM'

    call sd_clear_halos(grid,config_flags)
    call sd_coupling(grid)
    call sd_grid_bounds(grid)
    call sd_grid_sanity(grid)
    call sd_set_mpp_bounds(grid)
    call sd_halos(grid,all=.false.)
    if((allowed_to_read .and. .not. restrt) .or. anal) then
       call sd_refl(grid)
    end if
    call sd_j_const(grid)

    if((allowed_to_read.and..not.(restrt)) .or. .not.allowed_to_read) then
       call sd_random(grid)
    endif

    call sd_hires_diag(grid,config_flags)

#ifdef HWRF
    doit=(.NOT.RESTRT .AND. .NOT.ANAL) .OR. .NOT.allowed_to_read
#else
    doit=.not.RESTRT
#endif
    dome: IF(doit)THEN
       call sd_clean_start(grid)
    endif dome

    if(nest) call sd_fix_0_T(grid)

    grid%nphs0=GRID%NPHS

    call sd_restart_initializing(grid)
#ifndef NO_UPSTREAM_ADVECTION
    call sd_upstream_advection(grid)
#else
    CALL wrf_message( 'start_domain_nmm: upstream advection commented out')
#endif    

#ifdef HWRF
    doit= (nstart==0 .or. .not. allowed_to_read)
#else
    doit= (nstart==0)
#endif
    if(doit) call sd_zero_diag_arrays(grid)

    call sd_init_radtn(grid)
    if(allowed_to_read .and. (.not.RESTRT)) then
       call sd_init_turbl(grid)
    endif

#ifdef HWRF
    doit = .not.restrt .or. .not.allowed_to_read
#else
    doit= .not. restrt
#endif
    if(doit) call sd_init_old(grid)

    call sd_reset_nonhydro(grid,nest)

#ifndef NO_RESTRICT_ACCEL
    call sd_restrict_accel(grid)
#endif

    !------------------------------------------------
    ! allocate massive amounts of data in this region
    call sd_physics_init(grid)
    call sd_copy_or_reset_phys(grid)
    IF (.NOT. RESTRT .and. ALLOWED_TO_READ) THEN
       call sd_calc_micro(grid)
    ENDIF

#ifdef HWRF
    ! FIXME: REMOVE THIS.  IT DOES NOTHING
    !zhang's doing
    grid%julyr_rst=grid%julyr_rst
    grid%julday_rst=grid%julday_rst
    grid%gmt_rst=grid%gmt_rst
    !end of zhang's doing
#endif

    call sd_replace_data(grid)
    call sd_reset_precip(grid)
    call sd_compte_ct_ct2(grid)
    call sd_init_gwd(grid)
    call sd_set_time(grid)
    call sd_set_flags_after_phys_init(grid)

    call sd_dealloc(grid)
    ! have deallocated huge data arrays
    !------------------------------------------------

    ! do halos
    call sd_halos(grid,all=.true.) ! all = extra halo for MOIST & SCALAR

    write(message,*) "Timing for start_domain on d",grid%id
    call end_timing(message)

  END SUBROUTINE start_domain_main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module module_start_domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
