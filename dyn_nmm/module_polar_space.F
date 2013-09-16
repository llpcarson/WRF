module module_polar_space
  implicit none
  
  type polar_space
     real, pointer :: radius(:,:), angle(:,:)
     integer :: refs
  end type polar_space

  type polar_var
     type(polar_space), pointer :: space
     integer :: Nr, Na
     integer, pointer :: have(Nr,Na)
     real, pointer :: data(Nr,Na)
  end type polar_var

contains

  subroutine polar_space_init(space, grid, clon, clat, rotated, &
       ids,ide, jds,jde, kds,kde, &
       ims,ime, jms,jme, kms,kme, &
       its,ite, jts,jte, kts,kte)

    logical, intent(in) :: rotated
    type(polar_space), intent(inout) :: space
    integer, intent(in) :: Nr, Na

    ! Earth equatorial radius and Requator/(Requator-Rpole)
    real, parameter :: Requator=6378137., flattening_inv=298.247
    real, parameter :: pi180=0.01745329251 ! pi/180

    ! Earth (geographic) radius and angles:
    real :: Rearth1,Rearth2,Rearth,xlon1,xlon2,ylat1,ylat2

    allocate(space%radius(ims:ime,jms:jme))
    allocate(space%angle(ims:ime,jms:jme))

    ! Clean up the center lat & lon:
    xlon1=(mod(clon+3600.+180.,360.)-180.)
    ylat1=(mod(clat+3600.+180.,360.)-180.)
    if(clat>90.) then
       ylat1=180.-ylat1
       xlon1=mod(xlon1+360.,360.)-180.
    elseif(clat<-90.) then
       ylat1=-180. - ylat1
       xlon1=mod(xlon1+360.,360.)-180.
    endif
    xlon1=xlon1*pi180
    ylat1=ylat1*pi180

    ! Earth radius at the polar coordinate system center:
    Rearth1=Requator*(1-sin(ylat1)**2/flattening_inv)

    ! Loop over all points in this tile, calculating radius and angle.
    ! We do not clean up the tile lats & lons; instead we assume
    ! (hope?) the model did that for us.
    do j=jts,min(jde-1,jte)
       do i=its,min(ide-1,ite)
          xlon2=grid%glon(i,j)
          ylat2=grid%glat(i,j)

          ! Earth radius at the point of interest
          Rearth2=Requator*(1-sin(ylat2)**2/flattening_inv)          

          ! Distance (radius from polar C.S. center) using haversine
          ! formula, and the mean of the earth radii at the two points:
          radius(i,j)=(Rearth1+Rearth2)*asin(min(1., &
            sqrt(sin((ylat1-ylat2)/2)**2+ &
            cos(ylat1)*cos(ylat2)*sin((xlon1-xlon2)/2)**2)))

          ! Mathematical (NOT COMPASS!) angle in radians.  This is
          ! angle=0 at due East, angle=pi/2 at due North.  Do
          ! (pi/2-angle)*180/pi) to get compass angle in degrees:
          angle(i,j)=atan2(xlon2-xlon1,ylon2-ylon1)

          ! Yes, atan2 is a Fortran 95 intrinsic.  It shouldn't be, so
          ! I'm using it anyway.  Deal with it.
       enddo
    enddo
  end subroutine init_polar_space

  subroutine polar_space_mean(space, grid, var, id, &
       ids,ide, jds,jde, kds,kde, &
       ims,ime, jms,jme, kms,kme, &
       its,ite, jts,jte, kts,kte)
    real, intent(in) :: var(ims:ime,jms:jme,kms:kme)
    integer, intent(in) :: id
    real(kind=8) :: sum(space%Nr,space%Na)

    space%have=0
    space%data=0
    sum=0

    

  end subroutine fill_polar_space
end module module_polar_space
