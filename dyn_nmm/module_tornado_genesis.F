module module_tornado_genesis
  implicit none
  private

  public :: init_tornado_genesis, calc_tornado_genesis, reset_tornado_genesis

  real, parameter :: wwind_cutoff = 40000.0 ! pascals

contains

  subroutine update_tg_time(grid,init)
    ! Helper function that updates the three time interval variables
    ! based on the grid's clock.  If init=.true. then both times
    ! (interval start and end) are set to the current time, otherwise
    ! only the interval end is updated.  In either case, tg_duration
    ! is set to the length in seconds of the interval.
    use module_domain, only: domain, domain_get_time_since_sim_start
    use module_symbols_util, only: WRFU_TimeIntervalGet, WRFU_TimeInterval
    type(domain), intent(inout) :: grid
    type(WRFU_TimeInterval) :: since_start
    logical, intent(in) :: init
    integer :: s_i, s_n, s_d

    since_start=domain_get_time_since_sim_start(grid)
    s_i=0
    s_n=0
    s_d=1
    call WRFU_TimeIntervalGet(since_start,S=s_i,Sn=s_n,Sd=s_d)
    if(s_d==0) s_d=1
    grid%tg_interval_end=real(s_i) + real(s_n)/real(s_d)
    if(init)  grid%tg_interval_start=grid%tg_interval_end
    grid%tg_duration=grid%tg_interval_end-grid%tg_interval_start
  end subroutine update_tg_time

  subroutine init_tg_vars(grid,config_flags, &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         ips, ipe, jps, jpe, kps, kpe    )
    ! Helper function that resets all min/max accumulation arrays to 0
    use module_domain, only: domain, get_ijk_from_grid
    use module_configure, only : grid_config_rec_type
    use module_state_description, only: tg_emc2014spc
    type(domain), intent(inout) :: grid
    type(grid_config_rec_type), intent(in) :: config_flags
    integer, intent(in) :: IDS,IDE,JDS,JDE,KDS,KDE
    integer, intent(in) :: IMS,IME,JMS,JME,KMS,KME
    integer, intent(in) :: IPS,IPE,JPS,JPE,KPS,KPE
    integer :: i,j, istart,iend, jstart,jend

    if(config_flags%tg_option/=tg_emc2014spc) return
    jstart=max(jds,jps)
    jend=min(jpe,jde-1)
    istart=max(ids,ips)
    iend=min(ipe,ide-1)

    do j=jstart,jend
       do i=istart,iend
          grid%tg_max_m10wind(i,j)=0
          grid%tg_max_wwind(i,j)=0
          grid%tg_min_wwind(i,j)=0
          grid%tg_max_zhel_25(i,j)=0
          grid%tg_min_zhel_25(i,j)=0
          grid%tg_max_zhel_03(i,j)=0
          grid%tg_min_zhel_03(i,j)=0
          grid%tg_total_precip(i,j)=0
       enddo
    enddo
    
    grid%tg_interval_end=grid%tg_interval_start
    grid%tg_duration=0.0
  end subroutine init_tg_vars

  subroutine init_tornado_genesis(grid,config_flags)
    ! Called to initialize tornado genesis data arrays.  Should only
    ! be called at initial time.
    use module_domain, only: domain, get_ijk_from_grid
    use module_state_description, only: tg_emc2014spc
    use module_configure, only : grid_config_rec_type
    type(domain), intent(inout) :: grid
    type(grid_config_rec_type), intent(in) :: config_flags
    integer :: IDS,IDE,JDS,JDE,KDS,KDE
    integer :: IMS,IME,JMS,JME,KMS,KME
    integer :: IPS,IPE,JPS,JPE,KPS,KPE

    if(config_flags%tg_option/=tg_emc2014spc) return

    if(grid%hydro) then
       call wrf_error_fatal('Tornado genesis products require non-hydrostatic integration.')
    endif

    CALL get_ijk_from_grid (  grid ,      &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         ips, ipe, jps, jpe, kps, kpe    )
    call init_tg_vars(grid,config_flags, &
         ids, ide, jds, jde, kds, kde,   &
         ims, ime, jms, jme, kms, kme,   &
         ips, ipe, jps, jpe, kps, kpe    )
    call update_tg_time(grid,.true.)
  end subroutine init_tornado_genesis

  subroutine reset_tornado_genesis(grid,config_flags,stream)
    ! Called after writing output for a given stream.  Resets all
    ! min/max information for all fields if the stream is the
    ! tg_reset_stream.  
    use module_state_description, only: tg_emc2014spc
    use module_domain, only: domain, get_ijk_from_grid
    use module_configure, only : grid_config_rec_type
    use module_io_domain, only: first_history
    type(domain), intent(inout) :: grid
    type(grid_config_rec_type), intent(in) :: config_flags
    integer, intent(in) :: stream
    integer :: IDS,IDE,JDS,JDE,KDS,KDE
    integer :: IMS,IME,JMS,JME,KMS,KME
    integer :: IPS,IPE,JPS,JPE,KPS,KPE
    character*255 :: message
    integer :: histnum

    if(config_flags%tg_option/=tg_emc2014spc) return

    histnum=stream-first_history
    if(config_flags%tg_reset_stream == histnum) then
3012   format('Grid ',I2,': resetting tornado genesis data after stream ',I0,' output')
       write(message,3012) grid%id,histnum
       call wrf_message(trim(message))

       CALL get_ijk_from_grid (  grid ,      &
            ids, ide, jds, jde, kds, kde,    &
            ims, ime, jms, jme, kms, kme,    &
            ips, ipe, jps, jpe, kps, kpe    )
       call init_tg_vars(grid,config_flags, &
            ids, ide, jds, jde, kds, kde,   &
            ims, ime, jms, jme, kms, kme,   &
            ips, ipe, jps, jpe, kps, kpe    )

       ! Previous interval end time is now this interval's start time
       ! since we're entering the next interval:
       grid%tg_interval_start=grid%tg_interval_end
       grid%tg_duration=0.0 ! end=start, so interval length is 0
    else
3013   format('Grid ',I2,': NOT resetting tornado genesis data: stream ',I0,' is not reset stream ',I0)
       write(message,3013) grid%id,histnum,config_flags%tg_reset_stream
       call wrf_debug(1,trim(message))
    endif
  end subroutine reset_tornado_genesis

  subroutine calc_tornado_genesis(grid,config_flags)
    ! Updates max/min information for tornado genesis wind fields from
    ! grid data at the current time.  The tg_total_precip is handled
    ! in module_PHYSICS_CALLS instead.
    use module_comm_dm, only: HALO_NMM_C_sub
    use module_state_description, only: tg_emc2014spc
    use module_domain, only: domain, get_ijk_from_grid
    use module_configure, only : grid_config_rec_type
#ifdef DM_PARALLEL
    use module_dm, only: wrf_dm_maxval_real, wrf_dm_minval_real, &
         ntasks_x, ntasks_y, mytask, ntasks, local_communicator
#endif
    type(domain), intent(inout) :: grid
    type(grid_config_rec_type), intent(in) :: config_flags
    integer :: IDS,IDE,JDS,JDE,KDS,KDE
    integer :: IMS,IME,JMS,JME,KMS,KME
    integer :: IPS,IPE,JPS,JPE,KPS,KPE
    integer :: i,j,k, istart,iend, jstart,jend, a, imin,imax
    real :: dudy, dvdx, w, zhel, maxmaxwind, minminw, maxmaxw, sec
    real :: height, height1, height2, height0, maxmaxzhel, minminzhel
    character*255 :: message

    if(config_flags%tg_option/=tg_emc2014spc) return
    if(grid%hydro) then
       call wrf_error_fatal('Tornado genesis products require non-hydrostatic integration.')
    endif

    CALL get_ijk_from_grid (  grid ,      &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         ips, ipe, jps, jpe, kps, kpe    )
    jstart=max(jps,jds+1)
    jend=min(jpe,jde-2)
    istart=max(ips,ids+1)
    iend=min(ipe,ide-2)
    imin=max(ips,ids)
    imax=min(ipe,ide-1)

#ifdef DM_PARALLEL
#    include "HALO_NMM_C.inc"
#endif

    ! Maximum 10m wind vector magnitude:
    maxmaxwind=0.0
    do j=jstart,jend
       do i=istart,iend
          grid%tg_max_m10wind(i,j)=max(grid%tg_max_m10wind(i,j), &
               sqrt(grid%u10(i,j)*grid%u10(i,j) + grid%v10(i,j)*grid%v10(i,j)))
          maxmaxwind=max(maxmaxwind,grid%tg_max_m10wind(i,j))
       enddo
    enddo

#ifdef DM_PARALLEL
    call wrf_dm_maxval_real(maxmaxwind,i,j)
#endif

    ! Min/max vertical wind below 400mbar:
    minminw=0.0
    maxmaxw=0.0
    do j=jstart,jend
       do i=istart,iend
          kloop: do k=kds+1,kde-1
             if(grid%pint(i,j,k)<wwind_cutoff) exit kloop
             w=grid%w(i,j,k)
             grid%tg_min_wwind(i,j)=min(grid%tg_min_wwind(i,j),w)
             minminw=min(minminw,grid%tg_min_wwind(i,j))
             grid%tg_max_wwind(i,j)=max(grid%tg_max_wwind(i,j),w)
             maxmaxw=max(maxmaxw,grid%tg_max_wwind(i,j))
          enddo kloop
       enddo
    enddo

#ifdef DM_PARALLEL
    call wrf_dm_maxval_real(maxmaxw,i,j)
    call wrf_dm_minval_real(minminw,i,j)
#endif

    ! Min/max helicity for 0-3km layer and 2-5km layer.  Note this is
    ! X km above ground (lowest interface level geopotential height),
    ! not above sea level.
    minminzhel=0.0
    maxmaxzhel=0.0
    do j=jstart,jend
       a=mod(j,2)
       do i=istart,iend
          k=kds
          height0=grid%Z(i,j,k)   ! height0=lowest interface level height
          height1=0.0             ! height1=lower height bound for layer
          do while(k<kde-1 .and. height1<5000.0)
             height2=grid%Z(i,j,k+1)-height0 ! height2=layer upper height bound

             dvdx = (grid%v(i+1-a,j,k)-grid%v(i-a,j,k))/(2.*grid%dx_nmm(i,j))
             dudy = (grid%u(i,j+1,k)-grid%u(i,j-1,k))/(2.*grid%dy_nmm)
             zhel = (dvdx-dudy) * (grid%w(i,j,k) + grid%w(i,j,k+1))/2
             
             if(height1<3000.0) then
                grid%tg_max_zhel_03(i,j)=max(grid%tg_max_zhel_03(i,j),zhel)
                grid%tg_min_zhel_03(i,j)=min(grid%tg_min_zhel_03(i,j),zhel)
                minminzhel=min(grid%tg_min_zhel_03(i,j),minminzhel)
                maxmaxzhel=max(grid%tg_max_zhel_03(i,j),maxmaxzhel)
             endif
             if(height2>2000.0) then
                grid%tg_max_zhel_25(i,j)=max(grid%tg_max_zhel_25(i,j),zhel)
                grid%tg_min_zhel_25(i,j)=min(grid%tg_min_zhel_25(i,j),zhel)
                minminzhel=min(grid%tg_min_zhel_25(i,j),minminzhel)
                maxmaxzhel=max(grid%tg_max_zhel_25(i,j),maxmaxzhel)
             endif

             k=k+1
             height1=height2
          enddo
       enddo
    enddo

#ifdef DM_PARALLEL
    call wrf_dm_maxval_real(maxmaxzhel,i,j)
    call wrf_dm_minval_real(minminzhel,i,j)
#endif

    ! I boundaries copy from nearest point that has data, excluding corner points:
    if(ips<=ids) then
       grid%tg_max_zhel_25(ids,jstart:jend)=grid%tg_max_zhel_25(ids+1,jstart:jend)
       grid%tg_max_zhel_03(ids,jstart:jend)=grid%tg_max_zhel_03(ids+1,jstart:jend)
       grid%tg_min_zhel_25(ids,jstart:jend)=grid%tg_min_zhel_25(ids+1,jstart:jend)
       grid%tg_min_zhel_03(ids,jstart:jend)=grid%tg_min_zhel_03(ids+1,jstart:jend)
       grid%tg_max_wwind(ids,jstart:jend)=grid%tg_max_wwind(ids+1,jstart:jend)
       grid%tg_max_wwind(ids,jstart:jend)=grid%tg_max_wwind(ids+1,jstart:jend)
       grid%tg_max_m10wind(ids,jstart:jend)=grid%tg_max_m10wind(ids+1,jstart:jend)
    endif

    if(ipe>=ide-2) then
       grid%tg_max_zhel_25(ide-1,jstart:jend)=grid%tg_max_zhel_25(ide-2,jstart:jend)
       grid%tg_max_zhel_03(ide-1,jstart:jend)=grid%tg_max_zhel_03(ide-2,jstart:jend)
       grid%tg_min_zhel_25(ide-1,jstart:jend)=grid%tg_min_zhel_25(ide-2,jstart:jend)
       grid%tg_min_zhel_03(ide-1,jstart:jend)=grid%tg_min_zhel_03(ide-2,jstart:jend)
       grid%tg_max_wwind(ide-1,jstart:jend)=grid%tg_max_wwind(ide-2,jstart:jend)
       grid%tg_max_wwind(ide-1,jstart:jend)=grid%tg_max_wwind(ide-2,jstart:jend)
       grid%tg_max_m10wind(ide-1,jstart:jend)=grid%tg_max_m10wind(ide-2,jstart:jend)
    endif

    ! J boundaries: copy from nearest point that has data.  We use
    ! imin:imax instead of istart:iend to get the corner points.
    if(jps<=jds) then
       grid%tg_max_zhel_25(imin:imax,jds)=grid%tg_max_zhel_25(imin:imax,jds+1)
       grid%tg_max_zhel_03(imin:imax,jds)=grid%tg_max_zhel_03(imin:imax,jds+1)
       grid%tg_min_zhel_25(imin:imax,jds)=grid%tg_min_zhel_25(imin:imax,jds+1)
       grid%tg_min_zhel_03(imin:imax,jds)=grid%tg_min_zhel_03(imin:imax,jds+1)
       grid%tg_max_wwind(imin:imax,jds)=grid%tg_max_wwind(imin:imax,jds+1)
       grid%tg_max_wwind(imin:imax,jds)=grid%tg_max_wwind(imin:imax,jds+1)
       grid%tg_max_m10wind(imin:imax,jds)=grid%tg_max_m10wind(imin:imax,jds+1)
    endif

    if(jpe>=jde-2) then
       grid%tg_max_zhel_25(imin:imax,jde-1)=grid%tg_max_zhel_25(imin:imax,jde-2)
       grid%tg_max_zhel_03(imin:imax,jde-1)=grid%tg_max_zhel_03(imin:imax,jde-2)
       grid%tg_min_zhel_25(imin:imax,jde-1)=grid%tg_min_zhel_25(imin:imax,jde-2)
       grid%tg_min_zhel_03(imin:imax,jde-1)=grid%tg_min_zhel_03(imin:imax,jde-2)
       grid%tg_max_wwind(imin:imax,jde-1)=grid%tg_max_wwind(imin:imax,jde-2)
       grid%tg_max_wwind(imin:imax,jde-1)=grid%tg_max_wwind(imin:imax,jde-2)
       grid%tg_max_m10wind(imin:imax,jde-1)=grid%tg_max_m10wind(imin:imax,jde-2)
    endif

    call update_tg_time(grid,.false.)

3313 format('TG extrema: max(wind)=',F0.2,' max(w)=',F0.2,' min(w)=',F0.2,'  max(zhel)=',F0.4,' min(zhel)=',F0.4)
    write(message,3313) maxmaxwind,maxmaxw,minminw,maxmaxzhel,minminzhel
    call wrf_debug(1,message)
  end subroutine calc_tornado_genesis

end module module_tornado_genesis

subroutine nmm_reset_tornado_genesis(grid,config_flags,stream)
  ! This subroutine is a wrapper kludge to work around the WRF build
  ! order and limitations of make.  The module_tornado_genesis module
  ! file does not exist when mediation_integrate is compiled, so
  ! med_hist_out has to call a non-module function instead.
  use module_domain, only: domain
  use module_configure, only : grid_config_rec_type
  use module_tornado_genesis
  implicit none
  integer, intent(in) :: stream
  type(domain), intent(inout) :: grid
  type(grid_config_rec_type), intent(in) :: config_flags
  call reset_tornado_genesis(grid,config_flags,stream)
end subroutine nmm_reset_tornado_genesis
