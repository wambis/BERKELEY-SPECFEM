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

!-------------------------------
!
! 1D Berkeley model 
!
! Add infos... 
!
!-------------------------------

module model_1dberkeley_par

  ! number of layers in model1Dberkeley.dat
  integer :: NR_REF_BERKELEY
  integer :: NR_inner_core_berk
  integer :: NR_outer_core_berk
  integer :: NR_water_berk
  integer :: ifanis_berk
  integer :: tref_berk
  integer :: ifdeck_berk

  ! model_1dberkeley_variables
  double precision, dimension(:), allocatable :: &
    Mref_V_radius_berkeley,                      &
    Mref_V_density_berkeley,                     &
    Mref_V_vpv_berkeley,                         &
    Mref_V_vph_berkeley,                         &
    Mref_V_vsv_berkeley,                         &
    Mref_V_vsh_berkeley,                         &
    Mref_V_eta_berkeley,                         &
    Mref_V_Qkappa_berkeley,                      &
    Mref_V_Qmu_berkeley

end module model_1dberkeley_par


  !--------------------------------------------
  subroutine model_1dberkeley_broadcast(myrank)
  !--------------------------------------------
  !
  ! reads and broadcasts berkeley 1D model
  !
  use constants
  use model_1dberkeley_par
  !
  implicit none
  !
  integer :: i, icode, ier, myrank
  integer, parameter :: lunit = 54
  character (len=100) :: filename, title
  !
  ! define the berkeley 1D model 
  !
  filename = 'DATA/berkeley_model/model1D.dat'
  !
  ! root nodes read header
  !
  if(myrank==0)then
    !
    open(lunit,file=trim(filename),status='old')
    !
    read(lunit,100,iostat=icode) title
    read(lunit,*  ,iostat=icode) ifanis_berk,&
                                 tref_berk,  &
                                 ifdeck_berk
    !                             
    read(lunit,*  ,iostat=icode) NR_REF_BERKELEY,    &
                                 NR_inner_core_berk, &
                                 NR_outer_core_berk, &
                                 NR_water_berk
    !
  endif
  !
  ! broadcast header values
  !
  !call MPI_BCAST(ifanis_berk       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(tref_berk         ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(ifdeck_berk       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(NR_REF_BERKELEY   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(NR_inner_core_berk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(NR_outer_core_berk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(NR_water_berk     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !
  call BCAST_ALL_SINGLEI(ifanis_berk       )
  call BCAST_ALL_SINGLEI(tref_berk         )
  call BCAST_ALL_SINGLEI(ifdeck_berk       )
  call BCAST_ALL_SINGLEI(NR_REF_BERKELEY   )
  call BCAST_ALL_SINGLEI(NR_inner_core_berk)
  call BCAST_ALL_SINGLEI(NR_outer_core_berk)
  call BCAST_ALL_SINGLEI(NR_water_berk     )
  !
  ! allocate arrays
  !
  allocate(Mref_V_radius_berkeley(NR_REF_BERKELEY),  &
           Mref_V_density_berkeley(NR_REF_BERKELEY), &
           Mref_V_vpv_berkeley(NR_REF_BERKELEY),     &
           Mref_V_vsv_berkeley(NR_REF_BERKELEY),     &
           Mref_V_Qkappa_berkeley(NR_REF_BERKELEY),  &
           Mref_V_Qmu_berkeley(NR_REF_BERKELEY),     &
           Mref_V_vph_berkeley(NR_REF_BERKELEY),     &
           Mref_V_vsh_berkeley(NR_REF_BERKELEY),     &
           Mref_V_eta_berkeley(NR_REF_BERKELEY)      )
  !
  ! root proc reads data
  !
  if(myrank==0)then
  !
  do i = 1,NR_REF_BERKELEY
    read(lunit,*)Mref_V_radius_berkeley(i),  &
                 Mref_V_density_berkeley(i), &
                 Mref_V_vpv_berkeley(i),     &
                 Mref_V_vsv_berkeley(i),     &
                 Mref_V_Qkappa_berkeley(i),  &
                 Mref_V_Qmu_berkeley(i),     &
                 Mref_V_vph_berkeley(i),     &
                 Mref_V_vsh_berkeley(i),     &
                 Mref_V_eta_berkeley(i)     
  enddo
  !
  close(lunit)
  !
  endif
  !
  ! broadcast data
  !
  call BCAST_ALL_DP(Mref_V_radius_berkeley ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_density_berkeley,NR_REF_BERKELEY) 
  call BCAST_ALL_DP(Mref_V_vpv_berkeley    ,NR_REF_BERKELEY)     
  call BCAST_ALL_DP(Mref_V_vph_berkeley    ,NR_REF_BERKELEY)     
  call BCAST_ALL_DP(Mref_V_vsv_berkeley    ,NR_REF_BERKELEY)     
  call BCAST_ALL_DP(Mref_V_vsh_berkeley    ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_eta_berkeley    ,NR_REF_BERKELEY)  
  call BCAST_ALL_DP(Mref_V_Qkappa_berkeley ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_Qmu_berkeley    ,NR_REF_BERKELEY)
  !
  ! reading formats
  !
100 format(a80)
105 format(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)   
  !
  !----------------------------------------
  end subroutine model_1dberkeley_broadcast
  !----------------------------------------
  



!
!----------------------------------------------
!

  subroutine model_1dberkeley(x,rho,vpv,vph,vsv,vsh,eta,Qkappa,Qmu,iregion_code,CRUSTAL)

  use constants
  use model_1dberkeley_par

  implicit none

! model_1dref_variables

! input:
! dimensionless radius x

! output: non-dimensionalized
!
! mass density             : rho
! compressional wave speed : vpv
! compressional wave speed : vph
! shear wave speed         : vsv
! shear wave speed         : vsh
! dimensionless parameter  : eta
! shear quality factor     : Qmu
! bulk quality factor      : Qkappa

  double precision :: x,rho,vpv,vph,vsv,vsh,eta,Qmu,Qkappa
  integer :: iregion_code
  logical :: CRUSTAL

  ! local parameters
  double precision :: r,frac,scaleval
  integer :: i
  logical, parameter :: mimic_native_specfem = .true. 

  ! compute real physical radius in meters
  r = x * R_EARTH

  i = 1
  do while(r >= Mref_V_radius_berkeley(i) .and. i /= NR_REF_BERKELEY)
    i = i + 1
  enddo

! make sure we stay in the right region
  if(mimic_native_specfem .and. iregion_code == IREGION_INNER_CORE .and. i > NR_inner_core_berk) i = NR_inner_core_berk

  if(mimic_native_specfem .and. iregion_code == IREGION_OUTER_CORE .and. i < NR_inner_core_berk+2) i = NR_inner_core_berk+2
  if(mimic_native_specfem .and. iregion_code == IREGION_OUTER_CORE .and. i > NR_outer_core_berk) i = NR_outer_core_berk

  if(mimic_native_specfem .and. iregion_code == IREGION_CRUST_MANTLE .and. i < NR_outer_core_berk+2) i = NR_outer_core_berk+2

  ! if crustal model is used, mantle gets expanded up to surface
  ! for any depth less than 24.4 km, values from mantle below moho are taken
  if(mimic_native_specfem .and. CRUSTAL .and. i > 713) i = 713 ! Warining : may need to be changed if file is modified !
  !
  if(i == 1) then
    ! first layer in inner core
    rho    = Mref_V_density_berkeley(i)
    vpv    = Mref_V_vpv_berkeley(i)
    vph    = Mref_V_vph_berkeley(i)
    vsv    = Mref_V_vsv_berkeley(i)
    vsh    = Mref_V_vsh_berkeley(i)
    eta    = Mref_V_eta_berkeley(i)
    Qkappa = Mref_V_Qkappa_berkeley(i)
    Qmu    = Mref_V_Qmu_berkeley(i)
  else
    ! interpolates between one layer below to actual radius layer,
    ! that is from radius_ref(i-1) to r using the values at i-1 and i
    frac = (r-Mref_V_radius_berkeley(i-1))/(Mref_V_radius_berkeley(i)-Mref_V_radius_berkeley(i-1))
    ! interpolated model parameters
    rho = Mref_V_density_berkeley(i-1)  + frac * (Mref_V_density_berkeley(i)- Mref_V_density_berkeley(i-1))
    vpv = Mref_V_vpv_berkeley(i-1)      + frac * (Mref_V_vpv_berkeley(i)    - Mref_V_vpv_berkeley(i-1)    )
    vph = Mref_V_vph_berkeley(i-1)      + frac * (Mref_V_vph_berkeley(i)    - Mref_V_vph_berkeley(i-1)    )
    vsv = Mref_V_vsv_berkeley(i-1)      + frac * (Mref_V_vsv_berkeley(i)    - Mref_V_vsv_berkeley(i-1)    )
    vsh = Mref_V_vsh_berkeley(i-1)      + frac * (Mref_V_vsh_berkeley(i)    - Mref_V_vsh_berkeley(i-1)    )
    eta = Mref_V_eta_berkeley(i-1)      + frac * (Mref_V_eta_berkeley(i)    - Mref_V_eta_berkeley(i-1)    )
    Qkappa = Mref_V_Qkappa_berkeley(i-1)+ frac * (Mref_V_Qkappa_berkeley(i) - Mref_V_Qkappa_berkeley(i-1) )
    Qmu = Mref_V_Qmu_berkeley(i-1)      + frac * (Mref_V_Qmu_berkeley(i)    - Mref_V_Qmu_berkeley(i-1)    )
  endif

  ! make sure Vs is zero in the outer core even if roundoff errors on depth
  ! also set fictitious attenuation to a very high value (attenuation is not used in the fluid)
  if(mimic_native_specfem .and. iregion_code == IREGION_OUTER_CORE) then
    vsv = 0.d0
    vsh = 0.d0
    Qkappa = 3000.d0
    Qmu = 3000.d0
  endif

  ! non-dimensionalize
  ! time scaling (s^{-1}) is done with scaleval
  scaleval=dsqrt(PI*GRAV*RHOAV)
  rho=rho/RHOAV
  vpv=vpv/(R_EARTH*scaleval)
  vph=vph/(R_EARTH*scaleval)
  vsv=vsv/(R_EARTH*scaleval)
  vsh=vsh/(R_EARTH*scaleval)

  end subroutine model_1dberkeley

