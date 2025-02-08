!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
! ETOPO
!
! Global Gridded Elevation Data
!
! by default (constants.h), it uses a smoothed ETOPO 4 dataset
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------
module etopo5
!--------------------------------------------------------------------------
  implicit none
  public :: load_etopo5_dos,filtre_topo,write_topo_filtre,init_etopo_filtre,etopo_filtre
  private
  integer, parameter ::  NBT_TOPO5=2160,NBP_TOPO5=4320,INITR_TOPO5=12
  integer, dimension(:,:), allocatable :: topo_start
  integer :: NBT,NBP,INITR
  real, parameter :: pi=3.141592653589793,deg2rad=pi/180.
  real :: drfiltre
!************************************************************************
  contains
!************************************************************************
!-----------------------------------------------------------------
  subroutine filtre_topo(tout,dr,NBTF,NBPF)
!dr: rayon du filtrage
!NBTF: nombre d'echantillons sur un meridien (de 0 a 180 + 1 echantillon)
!NBPF: nombre d'echantillons sur un paralelle (de 0 a 360)
!-----------------------------------------------------------------
    implicit none
    integer , intent(in) :: NBTF,NBPF
    real, intent(in) :: dr
    integer, dimension(:,:), intent(out) :: tout    
    integer :: i,j
    real :: theta,phi
    do j=1,NBTF
       print*,'filtering ',j,' ...'
       theta=real(j-1)*180/real(NBTF-1)
       do i=1,NBPF
          phi=real(i-1)*360/real(NBPF)
          tout(i,j)=gauss_filtre(topo_start,theta,phi,dr)
       enddo
    enddo
!-----------------------------------------------------------------
  end subroutine filtre_topo
!-----------------------------------------------------------------
!-----------------------------------------------------------------
  subroutine init_etopo_filtre(file,dr_)
!-----------------------------------------------------------------
    implicit none
    character(len=*), intent(in) :: file
    doubleprecision :: dr_
!
    open(52,file=file,status='old',action='read')
    call read_topo_filtre(52)
    drfiltre=dr_
    close(52)
!-----------------------------------------------------------------
  end subroutine init_etopo_filtre
!-----------------------------------------------------------------
!-----------------------------------------------------------------
  integer function etopo_filtre(theta,phi)
!theta phi en radians
!dr en degre
!reponse en metre
!-----------------------------------------------------------------
    implicit none
    doubleprecision :: theta,phi
    real::t,p
    t=theta/deg2rad
    p=phi  /deg2rad
    etopo_filtre=gauss_filtre(topo_start,t,p,drfiltre)
!-----------------------------------------------------------------
  end function etopo_filtre
!-----------------------------------------------------------------

!-----------------------------------------------------------------
  subroutine write_topo_filtre(unit,tin,dr,NBTF,NBPF)
!-----------------------------------------------------------------
    implicit none
    integer , intent(in) :: NBTF,NBPF,unit
    real, intent(in) :: dr
    integer, dimension(:,:), intent(in) :: tin
    integer :: i,j
!
    write(unit,*) NBPF,NBTF,dr
    do i=1,NBTF
       do j=1,NBPF
          write(unit,'(i5)') tin(j,i)
       enddo
    enddo
!-----------------------------------------------------------------
  end subroutine write_topo_filtre
!-----------------------------------------------------------------

!-----------------------------------------------------------------
  subroutine read_topo_filtre(unit)
!-----------------------------------------------------------------
    implicit none
    integer, intent(in) :: unit
    real :: dr
    integer :: i,j
!
    
    read(unit,*) NBP,NBT,dr
    if (dr>1.) STOP 'read_topo_filtre: dr muste be <=1'
    INITR=int(1./dr)
    allocate(topo_start(NBP,NBT))
    do i=1,NBT
       do j=1,NBP
          read(unit,'(i5)') topo_start(j,i)
       enddo
    enddo
    NBT=NBT-1 !pour etre coherent avec load_etopo5_dos (181-1=180)
!-----------------------------------------------------------------
  end subroutine read_topo_filtre
!-----------------------------------------------------------------
!-----------------------------------------------------------------
  subroutine load_etopo5_dos()
!pour lire etopo5.dos.
!attention doit etre complier en little_endian et 
!le recl unit doit etre le byte (et pas 4 bytes)
!-----------------------------------------------------------------
    implicit none
    integer*2 :: val
    integer :: unit,i,j,k
    unit=11
    k=0
    NBT=NBT_TOPO5
    NBP=NBP_TOPO5
    INITR=INITR_TOPO5
    allocate(topo_start(NBP_TOPO5,NBT_TOPO5+1))
    open(unit,file='ETOPO5.DOS',status='old',form='unformatted',access='direct',recl=2)
    do i=1,NBT ! 180*12
       if (mod(i,100)==0) print*,'Reading ETOPO5.DOS for i=',i ,' ...'
       do j=1,NBP ! 360*12
          k=k+1
          read(unit,rec=k,err=100) val
          if (i==1 .and. j==1 .and. val /= -4290 ) STOP 'ETOPO5.DOS reading problem: use little_endian!'
          topo_start(j,i)=val
       enddo
    enddo
    topo_start(:,NBT+1)=2810
    close(unit)
    return
100 continue
    print*,'fin de fichier prematuree pour ETOPO5.DOS'
!-----------------------------------------------------------------
  end subroutine load_etopo5_dos
!-----------------------------------------------------------------
!----------------------------------------------------------------------
  real function cos_cylindre(t0_,p0_,d0_,theta_,phi_)
!----------------------------------------------------------------------
    implicit none
    real :: t0,p0,d0,theta,phi, d_ang
    real :: t0_,p0_,d0_,theta_,phi_
    !
    ! sfrench
    real :: d_arg_clean
    ! sfrench
    !
!
    t0=t0_*deg2rad
    p0=p0_*deg2rad
    theta=theta_*deg2rad
    phi=phi_*deg2rad
    d0=d0_*deg2rad
!distance angulaire au centre du cylindre:
    !
    ! sfrench
    !d_ang=acos(cos(theta)*cos(t0)+sin(theta)*sin(t0)*cos(phi-p0))
    d_arg_clean = min( 1.0, max( 0.0, real( cos(theta)*cos(t0)+sin(theta)*sin(t0)*cos(phi-p0) ) ) )
    d_ang       = acos( d_arg_clean )
    ! sfrench
    !
    if (d_ang>d0) then
       cos_cylindre=0.d0
    else
       cos_cylindre=0.5d0*(1.d0+cos(PI*d_ang/d0))
    endif
!----------------------------------------------------------------------
  end function cos_cylindre
!----------------------------------------------------------------------
!-----------------------------------------------------------------
  integer function gauss_filtre(tin,theta,phi,dr)
!-----------------------------------------------------------------
    integer, dimension(:,:)  :: tin
    real :: theta,phi,dr,tmp,thetar,phir,tmpnorm,int
    integer :: i,ii,j,jj,LARG
!
    tmp=0.
    tmpnorm=0. 
    LARG=2
    do i=1,LARG*dr*INITR+1
       do j=1,LARG*dr*INITR+1
          call get_indexloc(phi,theta,i,j,dr,LARG,ii,jj,phir,thetar)
          int=cos_cylindre(theta,phi,dr,thetar,phir)*(1./real(INITR))**2*sin(thetar*deg2rad)
          tmp=tmp+tin(ii,jj)*int
          tmpnorm=tmpnorm+int
       enddo
    enddo
    gauss_filtre=tmp/tmpnorm
!-----------------------------------------------------------------
  end function gauss_filtre
!-----------------------------------------------------------------

!-----------------------------------------------------------------
  subroutine get_indexloc(phi,theta,i,j,dr,LARG,ii,jj,phir,thetar)
!-----------------------------------------------------------------
    implicit none
    real, intent(in) :: theta,phi,dr
    integer, intent(in) :: i,j,LARG
    doubleprecision :: t,p,eps
    real, intent(out) :: thetar,phir
    integer, intent(out) :: ii,jj
    eps=1.d-8
!       
    p  =phi+(i-1-LARG*dr*INITR/2)/real(INITR) 
    t  =theta+(j-1-LARG*dr*INITR/2)/real(INITR)
    if (p<0.d0-eps) p=p+360.d0
    if (p>=360.d0-eps) p=p-360.d0
    if (t>180.d0-eps) then
       t=t-180.d0
       p=360.d0-p
    else if (t< 0.d0-eps) then
       t=180.d0+t
       p=360.d0-p
    endif
    if (p<0.d0-eps) p=p+360.d0
    if (p>=360.d0-eps) p=p-360.d0
!
    ii=nint(p*INITR)+1
    if (ii>NBP) ii=NBP
    jj=nint(t*INITR)+1
    if (jj>NBT+1) jj=NBT+1
    thetar=t
    phir  =p
!-----------------------------------------------------------------
  end subroutine get_indexloc
!-----------------------------------------------------------------
!--------------------------------------------------------------------------
end module etopo5
!--------------------------------------------------------------------------

  subroutine model_topo_bathy_broadcast(myrank,ibathy_topo,LOCAL_PATH)

! standard routine to setup model

  use constants
  use etopo5

  implicit none

  integer :: myrank

  ! bathymetry and topography: use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  character(len=150) :: LOCAL_PATH

  double precision, parameter :: dr_ = 2.d0

  IF(PATHNAME_TOPO_FILE == 'DATA/berkeley_model/ETOPO5_1x1_filtre.dat')THEN       
    call init_etopo_filtre(PATHNAME_TOPO_FILE,dr_)
  else        
  if(myrank == 0) then
    ! user output
    write(IMAIN,*)
    write(IMAIN,*) 'incorporating topography'
    call flush_IMAIN()

    ! read/save topo file on master
    call read_topo_bathy_file(ibathy_topo)
    call save_topo_bathy_database(ibathy_topo,LOCAL_PATH)
  endif

  ! broadcast the information read on the master to the nodes
  call bcast_all_i(ibathy_topo,NX_BATHY*NY_BATHY)
  endif

  end subroutine model_topo_bathy_broadcast
!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_topo_bathy_file(ibathy_topo)
!
!---- read topography and bathymetry file once and for all
!
  use constants

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  ! local parameters
  real :: val
  integer :: ival
  integer :: itopo_x,itopo_y,ier
  character(len=150) :: topo_bathy_file
  integer,parameter :: TOPO_MINIMUM = - 10000 ! (depth in m )
  integer,parameter :: TOPO_MAXIMUM = + 10000 ! (height in m )

  call get_value_string(topo_bathy_file, 'model.topoBathy.PATHNAME_TOPO_FILE', PATHNAME_TOPO_FILE)

  ! reads in topography values from file
  open(unit=13,file=trim(topo_bathy_file),status='old',action='read',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening:',trim(topo_bathy_file)
    call exit_mpi(0,'error opening topography data file')
  endif

  ! reads in topography array
  do itopo_y=1,NY_BATHY
    do itopo_x=1,NX_BATHY
      read(13,*,iostat=ier) val

      ! checks
      if( ier /= 0 ) then
        print*,'error read topo_bathy: ix,iy = ',itopo_x,itopo_y,val
        print*,'topo_bathy dimension: nx,ny = ',NX_BATHY,NY_BATHY
        call exit_mpi(0,'error reading topo_bathy file')
      endif

      ! converts to integer
      ival = nint(val)

      ! checks values
      if( ival < TOPO_MINIMUM .or. ival > TOPO_MAXIMUM ) then
        print*,'error read topo_bathy: ival = ',ival,val,'ix,iy = ',itopo_x,itopo_y
        print*,'topo_bathy dimension: nx,ny = ',NX_BATHY,NY_BATHY
        call exit_mpi(0,'error reading topo_bathy file')
      endif

      ! stores in array
      ibathy_topo(itopo_x,itopo_y) = ival

    enddo
  enddo
  close(13)

  ! note: we check the limits after reading in the data. this seems to perform slightly faster
  !          however, reading ETOPO1.xyz will take ~ 2m 1.2s for a single process

  ! imposes limits
  if( USE_MAXIMUM_HEIGHT_TOPO .or. USE_MAXIMUM_DEPTH_OCEANS ) then
    do itopo_y=1,NY_BATHY
      do itopo_x=1,NX_BATHY

        ! impose maximum height of mountains, to suppress oscillations in Himalaya etc.
        if(USE_MAXIMUM_HEIGHT_TOPO .and. ibathy_topo(itopo_x,itopo_y) > MAXIMUM_HEIGHT_TOPO) &
          ibathy_topo(itopo_x,itopo_y) = MAXIMUM_HEIGHT_TOPO

        ! impose maximum depth of oceans, to suppress oscillations near deep trenches
        if(USE_MAXIMUM_DEPTH_OCEANS .and. ibathy_topo(itopo_x,itopo_y) < MAXIMUM_DEPTH_OCEANS) &
          ibathy_topo(itopo_x,itopo_y) = MAXIMUM_DEPTH_OCEANS

      enddo
    enddo

  endif

  ! user output
  write(IMAIN,*) "  topography/bathymetry: min/max = ",minval(ibathy_topo),maxval(ibathy_topo)

  end subroutine read_topo_bathy_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_topo_bathy_database(ibathy_topo,LOCAL_PATH)

  use constants

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo
  character(len=150) :: LOCAL_PATH

  ! local parameters
  character(len=150) :: prname
  integer :: ier

  ! create the name for the database of the current slide and region
  ! only master needs to save this
  call create_name_database(prname,0,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! saves topography and bathymetry file for solver
  open(unit=27,file=prname(1:len_trim(prname))//'topo.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) then
    ! inform about missing database topo file
    print*,'TOPOGRAPHY problem:'
    print*,'error opening file: ',prname(1:len_trim(prname))//'topo.bin'
    print*,'please check if path exists and rerun mesher'
    call exit_mpi(0,'error opening file for database topo')
  endif

  write(27) ibathy_topo
  close(27)

  end subroutine save_topo_bathy_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_topo_bathy_database(ibathy_topo,LOCAL_PATH)

  use constants
  use etopo5

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo
  character(len=150) :: LOCAL_PATH

  ! local parameters
  character(len=150) :: prname
  integer :: ier

  double precision, parameter :: dr_ = 2.d0

  IF(PATHNAME_TOPO_FILE == 'DATA/berkeley_model/ETOPO5_1x1_filtre.dat')THEN       
    print*,'reading berkeley topography...'
    call init_etopo_filtre(PATHNAME_TOPO_FILE,dr_)
  else   
  ! create the name for the database of the current slide and region
  ! only master needs to save this
  call create_name_database(prname,0,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! reads topography and bathymetry file from saved database file
  open(unit=27,file=prname(1:len_trim(prname))//'topo.bin', &
        status='unknown',form='unformatted',action='read',iostat=ier)
  if( ier /= 0 ) then
    ! inform user
    print*,'TOPOGRAPHY problem:'
    print*,'error opening file: ',prname(1:len_trim(prname))//'topo.bin'
    !print*,'please check if file exists and rerun solver'
    !call exit_mpi(0,'error opening file for database topo')

    ! read by original file
    print*,'trying original topography file...'
    call read_topo_bathy_file(ibathy_topo)

    ! saves database topo file for next time
    call save_topo_bathy_database(ibathy_topo,LOCAL_PATH)
  else
    ! database topo file exists
    read(27) ibathy_topo
    close(27)

    ! user output
    write(IMAIN,*) "  topography/bathymetry: min/max = ",minval(ibathy_topo),maxval(ibathy_topo)
    call flush_IMAIN()
  endif
  endif

  end subroutine read_topo_bathy_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_topo_bathy(xlat,xlon,value,ibathy_topo)

!
!---- get elevation or ocean depth in meters at a given latitude and longitude
!

  use constants
  use etopo5

  implicit none

  ! location latitude/longitude (in degree)
  double precision,intent(in):: xlat,xlon

  ! returns elevation (in meters)
  double precision,intent(out):: value

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY),intent(in) :: ibathy_topo

  ! local parameters
  integer :: iadd1,iel1
  double precision :: samples_per_degree_topo
  double precision :: xlo
  double precision :: lon_corner,lat_corner,ratio_lon,ratio_lat
  double precision :: theta,phi
  real, parameter :: deg2rad=pi/180.,rad2deg=180./pi
        

  IF(PATHNAME_TOPO_FILE == 'DATA/berkeley_model/ETOPO5_1x1_filtre.dat')THEN

    theta = (90-xlat)*deg2rad ! assumed lat range: [-90,90]
    xlo = xlon
    if(xlo <-180.d0) xlo = xlo + 360.d0
    if(xlo > 180.d0) xlo = xlo - 360.d0
    phi = xlo*deg2rad ! assumed lon range: [-180,180]
    value = dble(etopo_filtre(theta,phi))
  else   
  ! initializes elevation
  value = ZERO

  ! longitude within range [0,360] degrees
  xlo = xlon
  if(xlo < 0.d0) xlo = xlo + 360.d0
  if(xlo > 360.d0) xlo = xlo - 360.d0

  ! compute number of samples per degree
  samples_per_degree_topo = dble(RESOLUTION_TOPO_FILE) / 60.d0

  ! compute offset in data file and avoid edge effects
  iadd1 = 1 + int((90.d0-xlat)/samples_per_degree_topo)
  if(iadd1 < 1) iadd1 = 1
  if(iadd1 > NY_BATHY) iadd1 = NY_BATHY

  iel1 = int(xlo/samples_per_degree_topo)
  if(iel1 <= 0 .or. iel1 > NX_BATHY) iel1 = NX_BATHY

! Use bilinear interpolation rather nearest point interpolation

  ! convert integer value to double precision
  !  value = dble(ibathy_topo(iel1,iadd1))

  lon_corner=iel1*samples_per_degree_topo
  lat_corner=90.d0-iadd1*samples_per_degree_topo

  ratio_lon = (xlo-lon_corner)/samples_per_degree_topo
  ratio_lat = (xlat-lat_corner)/samples_per_degree_topo

  if(ratio_lon<0.0) ratio_lon=0.0
  if(ratio_lon>1.0) ratio_lon=1.0
  if(ratio_lat<0.0) ratio_lat=0.0
  if(ratio_lat>1.0) ratio_lat=1.0

  ! convert integer value to double precision
  if( iadd1 <= NY_BATHY-1 .and. iel1 <= NX_BATHY-1 ) then
    ! interpolates for points within boundaries
    value = dble(ibathy_topo(iel1,iadd1))*(1-ratio_lon)*(1.-ratio_lat) &
            + dble(ibathy_topo(iel1+1,iadd1))*ratio_lon*(1.-ratio_lat) &
            + dble(ibathy_topo(iel1+1,iadd1+1))*ratio_lon*ratio_lat &
            + dble(ibathy_topo(iel1,iadd1+1))*(1.-ratio_lon)*ratio_lat
  else if( iadd1 <= NY_BATHY-1 .and. iel1 == NX_BATHY ) then
    ! interpolates for points on longitude border
    value = dble(ibathy_topo(iel1,iadd1))*(1-ratio_lon)*(1.-ratio_lat) &
            + dble(ibathy_topo(1,iadd1))*ratio_lon*(1.-ratio_lat) &
            + dble(ibathy_topo(1,iadd1+1))*ratio_lon*ratio_lat &
            + dble(ibathy_topo(iel1,iadd1+1))*(1.-ratio_lon)*ratio_lat
  else
    ! for points on latitude boundaries
    value = dble(ibathy_topo(iel1,iadd1))
  endif
  endif

  end subroutine get_topo_bathy

! -------------------------------------------


