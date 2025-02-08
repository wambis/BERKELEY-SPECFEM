module mirror

use specfem_par
use specfem_par_crustmantle
use constants

! mirror receiver information
  integer nrec_mir,nrec_local_mir, irec_mir
  integer, dimension(:), allocatable :: islice_selected_rec_mir,ispec_selected_rec_mir,number_receiver_global_mir
  double precision, dimension(:), allocatable :: xi_receiver_mir,eta_receiver_mir,gamma_receiver_mir
  character(len=150) :: STATIONS_mir,rec_filename_mir
  double precision, dimension(:,:,:), allocatable :: nu_mir
  double precision, allocatable, dimension(:) :: stlat_mir,stlon_mir,stele_mir,stbur_mir
  character(len=MAX_LENGTH_STATION_NAME), dimension(:), allocatable  :: station_name_mir
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:), allocatable :: network_name_mir  
  character(len=150) :: dummystring
  integer :: ios
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: seismograms_mir
  character :: fname_mirror*100
  integer, parameter :: spln_order_mir = 3
  integer :: decim_fact_mir
  real, parameter :: period_mir = 5.e0
  integer :: rec_len_mir, count_rec_mir, j, rec_len_all_mir
  real, allocatable :: b_coef_mir(:),seismograms_tmp_mir(:,:,:),diag(:,:),ddiag(:)
  real, allocatable  :: seismograms_all_mir(:,:)
  integer, allocatable :: seismograms_sizes_mir(:),seismograms_displs_mir(:),nrec_slice_mir(:)

! Lagrange interpolators at receivers
  double precision, dimension(:,:), allocatable :: hxir_store_mir,hetar_store_mir,hgammar_store_mir

contains

subroutine setup_mirror()

implicit none

integer irec, ier, nrec_tot_found_mir, nrec_simulation_mir

! counts mirror receiver
  if (SIMULATION_TYPE == 1) then
    rec_filename_mir = 'DATA/MIRROR'
  ! get total number of receivers
  if(myrank == 0) then
    open(unit=IIN,file=rec_filename_mir,iostat=ios,status='old',action='read')
    nrec_mir = 0
    do while(ios == 0)
      read(IIN,"(a)",iostat=ios) dummystring
      if(ios == 0) nrec_mir = nrec_mir + 1
    enddo
    close(IIN)
  endif
  ! broadcast the information read on the master to the nodes
  call BCAST_ALL_SINGLEI(nrec_mir)

  if(nrec_mir < 1) call exit_MPI(myrank,trim(rec_filename_mir)//': need at least one receiver')
    ! allocate memory for receiver arrays
  allocate(islice_selected_rec_mir(nrec_mir), &
          ispec_selected_rec_mir(nrec_mir), &
          xi_receiver_mir(nrec_mir), &
          eta_receiver_mir(nrec_mir), &
          gamma_receiver_mir(nrec_mir), &
          station_name_mir(nrec_mir), &
          network_name_mir(nrec_mir), &
          stlat_mir(nrec_mir), &
          stlon_mir(nrec_mir), &
          stele_mir(nrec_mir), &
          stbur_mir(nrec_mir), &
          nu_mir(NDIM,NDIM,nrec_mir),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating mirror receiver arrays')
  endif

 !  receivers
  if(myrank == 0) then
    write(IMAIN,*)
    if (SIMULATION_TYPE == 1) then
      write(IMAIN,*) 'Total number of mirror receivers = ', nrec_mir
    endif
    write(IMAIN,*)
  endif

  ! locate receivers in the crust in the mesh
  call locate_receivers_mirror(myrank,DT,NSTEP,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,ibool_crust_mantle, &
                               xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                               xigll,yigll,zigll,trim(rec_filename_mir), &
                               nrec_mir,islice_selected_rec_mir,ispec_selected_rec_mir, &
                               xi_receiver_mir,eta_receiver_mir,gamma_receiver_mir,station_name_mir,network_name_mir, &
                               stlat_mir,stlon_mir,stele_mir,stbur_mir,nu_mir, &
                               yr_SAC,jda_SAC,ho_SAC,mi_SAC,sec_SAC,NPROCTOT_VAL,ELLIPTICITY_VAL,TOPOGRAPHY, &
                               theta_source(1),phi_source(1),rspl,espl,espl2,nspl, &
                               ibathy_topo,RECEIVERS_CAN_BE_BURIED,NCHUNKS_VAL)


  !call synchronize_all()
  !write(IMAIN,*) 'exit locate_receiver_mirror'       
  !call flush_IMAIN()
  
  ! count number of receivers located in this slice
  nrec_local_mir = 0
  if (SIMULATION_TYPE == 1) then
    nrec_simulation_mir = nrec_mir
    do irec = 1,nrec_mir
      if(myrank == islice_selected_rec_mir(irec)) nrec_local_mir = nrec_local_mir + 1
    enddo
  endif


  ! check that the sum of the number of receivers in each slice is nrec_mir
  call SUM_ALL_I(nrec_local_mir,nrec_tot_found_mir)
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'found a total of ',nrec_tot_found_mir,' mirror receivers in all slices'
    if(nrec_tot_found_mir /= nrec_simulation_mir) then
      call exit_MPI(myrank,'problem when dispatching the mirror receivers')
    else
      write(IMAIN,*) 'this total is okay'
    endif
  endif       

  !call synchronize_all()
  !write(IMAIN,*) 'exit locate_receiver_mirror 1',myrank       
  !call flush_IMAIN()

  ! allocates receiver interpolators
  if (nrec_local_mir > 0) then
    ! allocate Lagrange interpolators for receivers
    allocate(hxir_store_mir(nrec_local_mir,NGLLX), &
            hetar_store_mir(nrec_local_mir,NGLLY), &
            hgammar_store_mir(nrec_local_mir,NGLLZ),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating mirror receiver interpolators')
    ! define local to global receiver numbering mapping
    allocate(number_receiver_global_mir(nrec_local_mir),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating global mirror receiver numbering')
    ! define and store Lagrange interpolators at all the receivers

  !call synchronize_all()
  !write(IMAIN,*) 'exit locate_receiver_mirror 2',myrank
  !call flush_IMAIN()
    ! stores interpolators for receiver positions
    call setup_sources_receivers_intp_mir(myrank, &
                      xigll,yigll,zigll, &
                      SIMULATION_TYPE,nrec_mir,nrec_local_mir, &
                      islice_selected_rec_mir,number_receiver_global_mir, &
                      xi_receiver_mir,eta_receiver_mir,gamma_receiver_mir, &
                      hxir_store_mir,hetar_store_mir,hgammar_store_mir)

  !write(IMAIN,*) 'exit locate_receiver_mirror 3',myrank
  !call flush_IMAIN()
    ! allocate seismogram array
    if (SIMULATION_TYPE == 1) then
      allocate(seismograms_mir(NDIM,nrec_local_mir),stat=ier)
      if(ier /= 0) stop 'error while allocating mirror seismograms'
    endif
    ! initialize seismograms
    seismograms_mir(:,:) = 0._CUSTOM_REAL
  else
    ! allocate dummy array since we need it to pass as argument e.g. in write_seismograms() routine
    ! note: nrec_local_mir is zero, fortran 90/95 should allow zero-sized array allocation...
    allocate(seismograms_mir(NDIM,nrec_local_mir),stat=ier)
    if( ier /= 0) stop 'error while allocating zero mirror seismograms'
    allocate(number_receiver_global_mir(nrec_local_mir),stat=ier)
    if( ier /= 0) stop 'error while allocating zero mirror number_receiver_global'
  endif
  !write(IMAIN,*) 'exit locate_receiver_mirror 4',myrank       
  !call flush_IMAIN()

return

end subroutine setup_mirror

subroutine write_mirror()

implicit none

integer :: ier, i, j, k

   if(it==it_begin)then

     decim_fact_mir = (period_mir)/dt
     if(myrank==0)write(imain,*)'the decimation factor for a period of ',period_mir,'is : ',decim_fact_mir,dt,scale_t

        allocate(b_coef_mir((spln_order_mir+1)*decim_fact_mir+1),stat=ier)
         if( ier /= 0 )call exit_MPI(myrank,'error allocating mirror arrays')
        call bmn(b_coef_mir,decim_fact_mir,spln_order_mir)     
        count_rec_mir = 0 ! init decimate counter  
        if(decim_fact_mir>1)then
         allocate(seismograms_tmp_mir(NDIM,nrec_local_mir,spln_order_mir+1),stat=ier)
         if( ier /= 0 )call exit_MPI(myrank,'error allocating mirror arrays')
         seismograms_tmp_mir(:,:,:) = 0.d0
        endif
        inquire(iolength=rec_len_mir)seismograms_tmp_mir(:,:,1)
        write (fname_mirror,"(a,I3.3)") " MIRROR_FILES/mirror.",myrank
        if (rec_len_mir/=0) open(33,file=fname_mirror,access='direct',form="unformatted",recl=rec_len_mir,status='replace')
    endif 

   call compute_seismograms_mirror(nrec_local_mir,nrec_mir,displ_crust_mantle, &
                      nu_mir,hxir_store_mir,hetar_store_mir,hgammar_store_mir, &
                                               scale_displ,ibool_crust_mantle, &
                            ispec_selected_rec_mir,number_receiver_global_mir, &
                            seismograms_mir)               
           !
!
! filter (convolve with B-spline)
!
if(decim_fact_mir>1)then
 do j = 1,spln_order_mir+1
  i = mod(it-it_begin,decim_fact_mir)+(spln_order_mir+1-j)*decim_fact_mir+1
  seismograms_tmp_mir(:,:,j) = seismograms_tmp_mir(:,:,j)+b_coef_mir(i)*seismograms_mir(:,:) ! HERE !
 enddo
endif
!
! Decimate/write
!
if(decim_fact_mir==1)then
 if (rec_len_mir/=0) write(33,rec=it-it_begin+1)seismograms_tmp_mir(:,:,1)
else
 if(mod(it-it_begin,decim_fact_mir)==decim_fact_mir-1)then
    count_rec_mir = count_rec_mir+1  
    if (rec_len_mir/=0) write(33,rec=count_rec_mir)seismograms_tmp_mir(:,:,1)
    do j = 1,spln_order_mir
     seismograms_tmp_mir(:,:,j) = seismograms_tmp_mir(:,:,j+1)
    enddo
     seismograms_tmp_mir(:,:,spln_order_mir+1) = 0
 endif
endif
   call synchronize_all()
!
if (it==it_end)then
 ! dump buffer
    do j = 1,spln_order_mir+1
     count_rec_mir = count_rec_mir+1  
     if(rec_len_mir/=0)write (33,rec=count_rec_mir)seismograms_tmp_mir(:,:,j)
    enddo
   call synchronize_all()
 ! post processing : obtain B-spline coefs
   deallocate(seismograms_tmp_mir) 
   call synchronize_all()
   allocate(diag(spln_order_mir+1,count_rec_mir))
   allocate(ddiag(count_rec_mir))
   allocate(seismograms_tmp_mir(NDIM,nrec_local_mir,0:spln_order_mir))
   call synchronize_all()
    call fill_spmat_diag(decim_fact_mir,spln_order_mir,count_rec_mir,it_end-it_begin+1,diag)
    call bchfac ( diag, spln_order_mir+1, count_rec_mir, ddiag )
    if(rec_len_mir/=0)call bchslv (seismograms_tmp_mir, NDIM,nrec_local_mir, spln_order_mir, diag, spln_order_mir+1, count_rec_mir, 33)
   deallocate(diag)
   deallocate(ddiag)
   deallocate(seismograms_tmp_mir) 
   call synchronize_all()
 if(rec_len_mir/=0)close (33)
 deallocate(b_coef_mir)

 ! 
 ! Merge files
 !
 
 ! get sismo sizes

   call synchronize_all()
  if(myrank==0)then
   allocate(seismograms_sizes_mir(NPROCTOT_VAL))
   allocate(seismograms_displs_mir(NPROCTOT_VAL))
   allocate(nrec_slice_mir(NPROCTOT_VAL))
  endif  
   call synchronize_all()

   call gather_all_singlei(nrec_local_mir,seismograms_sizes_mir,NPROCTOT_VAL)

   IF(MYRANK==0)then
   do i = 1,NPROCTOT_VAL
        write(imain,*)i,seismograms_sizes_mir(i)
   enddo
   endif

 ! setup displacement

if(myrank==0)then
  if(sum(seismograms_sizes_mir)/=nrec_mir)write(IMAIN,*)'Problem number of mirror receiver doesn"t match'
  seismograms_sizes_mir(:) = seismograms_sizes_mir(:)*3
  seismograms_displs_mir(:) = seismograms_sizes_mir(:)
  do i = 2,NPROCTOT_VAL
    seismograms_displs_mir(i) = seismograms_displs_mir(i)+seismograms_displs_mir(i-1)
  enddo
  do i = NPROCTOT_VAL,2,-1
    seismograms_displs_mir(i) = seismograms_displs_mir(i-1)
  enddo
  seismograms_displs_mir(1) = 0
  allocate(seismograms_all_mir(NDIM,nrec_mir))
endif    
 
 IF(MYRANK==0)then
   do i = 1,NPROCTOT_VAL
        write(imain,*)i,seismograms_sizes_mir(i),seismograms_displs_mir(i)
   enddo
   endif
  ! open output file (main only)

  if(myrank == 0)then
    inquire(iolength=rec_len_all_mir)seismograms_all_mir(:,:)
    open(35,file='mirror_src.dat',access='stream',form="unformatted",status='replace')
    write(35)int(spln_order_mir),sngl(t0-(spln_order_mir+1)*dt*decim_fact_mir/2.),sngl(dt*decim_fact_mir),int(nrec_mir),int(count_rec_mir)
  endif 

 ! allocate array to store all seismo

 if (rec_len_mir/=0) open(33,file=fname_mirror,access='direct',form="unformatted",recl=rec_len_mir,status='old')
 do j = 1,count_rec_mir
    if (rec_len_mir/=0) read(33,rec=j)seismograms_mir(:,:) ! read seismo
    ! gather seismo
    !call mpi_gatherv(seismograms_mir,nrec_local_mir*3,mpi_real,seismograms_all_mir,seismograms_sizes_mir,seismograms_displs_mir,mpi_real,0,mpi_comm_world,ier)
    call gatherv_all_r(seismograms_mir,nrec_local_mir*3,seismograms_all_mir,seismograms_sizes_mir,seismograms_displs_mir,NDIM*nrec_mir,NPROCTOT_VAL)
    ! write seismo
!    if(myrank==0)write(35)seismograms_all_mir(:,:)
     if(myrank==0)then
       nrec_slice_mir(:) = seismograms_sizes_mir(:)/3
       do i = 2,NPROCTOT_VAL
         nrec_slice_mir(i) = nrec_slice_mir(i)+nrec_slice_mir(i-1)
       enddo
       do i = NPROCTOT_VAL,2,-1
         nrec_slice_mir(i) = nrec_slice_mir(i-1)
       enddo
       nrec_slice_mir(1) = 0
       do i = 1,nrec_mir
          irec_mir = islice_selected_rec_mir(i)+1    
          nrec_slice_mir(irec_mir) = nrec_slice_mir(irec_mir) + 1
          write(35)seismograms_all_mir(:,nrec_slice_mir(irec_mir))
       enddo
       write(imain,*)'Merging mirror files... writing record:',j,' to file'
     !
     ! check compare sizes and nrec_slice
     !
       do i = 1,NPROCTOT_VAL
!         if(nrec_slice_mir(i)-seismograms_displs_mir(i)/=seismograms_sizes_mir(i))write(IMAIN,*)'sizes doesn"t match....'
       enddo
     endif
 enddo
 if (rec_len_mir/=0) close(33)
  ! close output file (main only)

  if(myrank == 0)then
    close(35)
    deallocate(seismograms_sizes_mir)
    deallocate(seismograms_displs_mir)
    deallocate(seismograms_all_mir)
    deallocate(nrec_slice_mir)
  endif 
endif


return

end subroutine write_mirror

end module mirror 

  subroutine compute_seismograms_mirror(nrec_local,nrec,displ_crust_mantle, &
                                nu,hxir_store,hetar_store,hgammar_store, &
                                scale_displ,ibool_crust_mantle, &
                                ispec_selected_rec,number_receiver_global, &
                                seismograms)

  implicit none
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer nrec_local,nrec
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
    displ_crust_mantle

  double precision, dimension(NDIM,NDIM,nrec) :: nu

  double precision, dimension(nrec_local,NGLLX) :: hxir_store
  double precision, dimension(nrec_local,NGLLY) :: hetar_store
  double precision, dimension(nrec_local,NGLLZ) :: hgammar_store

  double precision scale_displ

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  integer, dimension(nrec) :: ispec_selected_rec
  integer, dimension(nrec_local) :: number_receiver_global

  integer :: seismo_current
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS

  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local) :: &
    seismograms

  ! local parameters
  double precision :: uxd,uyd,uzd,hlagrange
  integer :: i,j,k,iglob,irec_local,irec

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ! perform the general interpolation using Lagrange polynomials
    uxd = ZERO
    uyd = ZERO
    uzd = ZERO

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          iglob = ibool_crust_mantle(i,j,k,ispec_selected_rec(irec))

          hlagrange = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)

          uxd = uxd + dble(displ_crust_mantle(1,iglob))*hlagrange
          uyd = uyd + dble(displ_crust_mantle(2,iglob))*hlagrange
          uzd = uzd + dble(displ_crust_mantle(3,iglob))*hlagrange

        enddo
      enddo
    enddo
    ! store North, East and Vertical components

    ! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
     ! seismograms(:,irec_local) = sngl(scale_displ*(nu(:,1,irec)*uxd + &
     !            nu(:,2,irec)*uyd + nu(:,3,irec)*uzd))
     seismograms(1,irec_local) = sngl(scale_displ*uxd)
     seismograms(2,irec_local) = sngl(scale_displ*uyd)
     seismograms(3,irec_local) = sngl(scale_displ*uzd)
    else
     ! seismograms(:,irec_local) = scale_displ*(nu(:,1,irec)*uxd + &
     !            nu(:,2,irec)*uyd + nu(:,3,irec)*uzd)
     seismograms(1,irec_local) = scale_displ*uxd
     seismograms(2,irec_local) = scale_displ*uyd
     seismograms(3,irec_local) = scale_displ*uzd
    endif

  enddo

  end subroutine compute_seismograms_mirror


 subroutine locate_receivers_mirror(myrank,DT,NSTEP,nspec,nglob,ibool, &
                                     xstore,ystore,zstore,xigll,yigll,zigll,rec_filename, &
                                     nrec,islice_selected_rec,ispec_selected_rec, &
                                     xi_receiver,eta_receiver,gamma_receiver,station_name,network_name, &
                                     stlat,stlon,stele,stbur,nu, &
                                     yr,jda,ho,mi,sec,NPROCTOT,ELLIPTICITY,TOPOGRAPHY, &
                                     theta_source,phi_source,rspl,espl,espl2,nspl, &
                                     ibathy_topo,RECEIVERS_CAN_BE_BURIED,NCHUNKS)

  use constants

  implicit none
  
  integer NPROCTOT,NCHUNKS

  logical ELLIPTICITY,TOPOGRAPHY,RECEIVERS_CAN_BE_BURIED

  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

  integer nspec,nglob,nrec,myrank,nrec_found

  integer yr,jda,ho,mi
  double precision sec

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)
  integer NSTEP
  double precision DT

! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ)

  character(len=*)  rec_filename

! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  integer, allocatable, dimension(:) :: ix_initial_guess,iy_initial_guess,iz_initial_guess

  integer iorientation
  integer iprocloop
  double precision stazi,stdip

  double precision, allocatable, dimension(:) :: x_target,y_target,z_target
  double precision, allocatable, dimension(:) :: epidist
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision, allocatable, dimension(:,:) :: x_found_all,y_found_all,z_found_all

  integer irec
  integer i,j,k,ispec,iglob
  integer ier

  double precision ell
  double precision elevation
  double precision n(3)
  double precision thetan,phin
  double precision sint,cost,sinp,cosp
  double precision r0,p20
  double precision theta,phi
  double precision theta_source,phi_source
  double precision dist
  double precision xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma

! topology of the control points of the surface element
  integer iax,iay,iaz
  integer iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer iter_loop,ispec_iterate

  integer ia
  double precision x,y,z
  double precision xix,xiy,xiz
  double precision etax,etay,etaz
  double precision gammax,gammay,gammaz

! timer MPI
  double precision time_start,tCPU

! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision, dimension(:,:), allocatable :: final_distance_all
  double precision distmin,final_distance_max

! receiver information
! timing information for the stations
! station information for writing the seismograms
  integer nsamp
  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(3,3,nrec) :: nu
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer, dimension(nrec) :: islice_selected_rec_found,ispec_selected_rec_found
  double precision, dimension(nrec) :: xi_receiver_found,eta_receiver_found,gamma_receiver_found
  double precision, dimension(3,3,nrec) :: nu_found
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name_found
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name_found
  double precision, dimension(nrec) :: stlat_found,stlon_found,stele_found,stbur_found,epidist_found
  character(len=150) STATIONS

  integer, allocatable, dimension(:,:) :: ispec_selected_rec_all
  double precision, dimension(nrec) :: stlat,stlon,stele,stbur
  double precision, allocatable, dimension(:,:) :: xi_receiver_all,eta_receiver_all,gamma_receiver_all

  character(len=150) OUTPUT_FILES
  character(len=2) bic

! added by yder masson

integer :: ix_initial_guess_selected(27),iy_initial_guess_selected(27),iz_initial_guess_selected(27),ispec_selected(27)
integer :: nspec_selected,iglob_selected,iselected
double precision, dimension(27) :: xi_receiver_selected,eta_receiver_selected,gamma_receiver_selected
double precision, dimension(27) :: x_found_selected,y_found_selected,z_found_selected,final_distance_selected
double precision :: distmin27(27)
integer :: ispec27,ispec_selected_27(27),idistmin
logical :: inside
logical, allocatable :: found(:)
integer, allocatable :: inside_all(:)
double precision :: x_rec_min,x_rec_max,y_rec_min,y_rec_max,z_rec_min,z_rec_max
double precision :: x_elm_min,x_elm_max,y_elm_min,y_elm_max,z_elm_min,z_elm_max
double precision :: l_elm_x,l_elm_y,l_elm_z,distance
integer :: n_rec_in_elm,n_rec_out_elm
real :: r_max,r_min
double precision :: wtime


! **************
! get MPI starting time
  time_start = wtime()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) ' locating receivers mirror'
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
  endif

! define topology of the control element
  call hex_nodes(iaddx,iaddy,iaddr)

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) 'reading mirror information'
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
  endif

! allocate memory for arrays using number of stations
  allocate(epidist(nrec), &
          inside_all(nrec), &
          ix_initial_guess(nrec), &
          iy_initial_guess(nrec), &
          iz_initial_guess(nrec), &
          x_target(nrec), &
          y_target(nrec), &
          z_target(nrec), &
          x_found(nrec), &
          y_found(nrec), &
          z_found(nrec), &
          final_distance(nrec), &
          ispec_selected_rec_all(nrec,0:NPROCTOT-1), &
          xi_receiver_all(nrec,0:NPROCTOT-1), &
          eta_receiver_all(nrec,0:NPROCTOT-1), &
          gamma_receiver_all(nrec,0:NPROCTOT-1), &
          x_found_all(nrec,0:NPROCTOT-1), &
          y_found_all(nrec,0:NPROCTOT-1), &
          z_found_all(nrec,0:NPROCTOT-1), &
          found(nrec), &
          final_distance_all(nrec,0:NPROCTOT-1),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating temporary mirror receiver arrays')

  ! read that STATIONS file on the master
  if(myrank == 0) then
    open(unit=1,file=rec_filename,status='old',action='read',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening MIRROR file')

    ! loop on all the stations to read station information
    do irec = 1,nrec
      read(1,*,iostat=ier)x_target(irec), y_target(irec), z_target(irec)! station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
      if( ier /= 0 ) then
        write(IMAIN,*) 'error reading in mirror station ',irec
        call exit_MPI(myrank,'error reading in mirror station in MIRROR file')
      endif
    enddo
    ! close receiver file
    close(1)
!
  endif

! broadcast the information read on the master to the nodes
!
  call BCAST_ALL_DP(x_target,nrec)
  call BCAST_ALL_DP(y_target,nrec)
  call BCAST_ALL_DP(z_target,nrec)
!
! find mirror box extent
!
x_target(:) = x_target(:)/R_EARTH
y_target(:) = y_target(:)/R_EARTH
z_target(:) = z_target(:)/R_EARTH
!
x_rec_min = minval(x_target) 
x_rec_max = maxval(x_target)
y_rec_min = minval(y_target)
y_rec_max = maxval(y_target)
z_rec_min = minval(z_target)
z_rec_max = maxval(z_target)
!
! Find elements that are inside the mirror box
!
final_distance(:) = HUGEVAL
found(:) = .false.
ispec_selected_rec(:) = 0
xi_receiver(:)    = 2
eta_receiver(:)   = 2
gamma_receiver(:) = 2
x_found(:)        = 0
y_found(:)        = 0
z_found(:)        = 0
!

do ispec=1,nspec
     do ia=1,NGNOD

      if(iaddx(ia) == 0) then
        iax = 1
      else if(iaddx(ia) == 1) then
        iax = (NGLLX+1)/2
      else if(iaddx(ia) == 2) then
        iax = NGLLX
      else
        call exit_MPI(myrank,'incorrect value of iaddx')
      endif

      if(iaddy(ia) == 0) then
        iay = 1
      else if(iaddy(ia) == 1) then
        iay = (NGLLY+1)/2
      else if(iaddy(ia) == 2) then
        iay = NGLLY
      else
        call exit_MPI(myrank,'incorrect value of iaddy')
      endif

      if(iaddr(ia) == 0) then
        iaz = 1
      else if(iaddr(ia) == 1) then
        iaz = (NGLLZ+1)/2
      else if(iaddr(ia) == 2) then
        iaz = NGLLZ
      else
        call exit_MPI(myrank,'incorrect value of iaddr')
      endif

      iglob = ibool(iax,iay,iaz,ispec)
      xelm(ia) = dble(xstore(iglob))
      yelm(ia) = dble(ystore(iglob))
      zelm(ia) = dble(zstore(iglob))

    enddo
!    
    x_elm_min = minval(xelm(1:NGNOD))
    x_elm_max = maxval(xelm(1:NGNOD))
    y_elm_min = minval(yelm(1:NGNOD))
    y_elm_max = maxval(yelm(1:NGNOD))
    z_elm_min = minval(zelm(1:NGNOD))
    z_elm_max = maxval(zelm(1:NGNOD))
!
    l_elm_x = x_elm_max-x_elm_min 
    l_elm_y = y_elm_max-y_elm_min 
    l_elm_z = z_elm_max-z_elm_min
!
    x_elm_min = x_elm_min-0.25*l_elm_x
    x_elm_max = x_elm_max+0.25*l_elm_x
    y_elm_min = y_elm_min-0.25*l_elm_y
    y_elm_max = y_elm_max+0.25*l_elm_y
    z_elm_min = z_elm_min-0.25*l_elm_z
    z_elm_max = z_elm_max+0.25*l_elm_z
!
! check if element box is inside or overlap with mirror box
!
    if(.not.(x_elm_max<=x_rec_min.or.x_elm_min>=x_rec_max.OR.&
             y_elm_max<=y_rec_min.or.y_elm_min>=y_rec_max.OR.&
             z_elm_max<=z_rec_min.or.z_elm_min>=Z_rec_max))then
!
      do irec = 1,nrec ! loop over receiver
      ! check if receiver is element box
        if(x_target(irec)<=x_elm_max.and.x_target(irec)>=x_elm_min.AND.&
          y_target(irec)<=y_elm_max.and.y_target(irec)>=y_elm_min.AND.&
          z_target(irec)<=z_elm_max.and.z_target(irec)>=z_elm_min)then
!
! use initial guess in xi and eta
    xi    = 0
    eta   = 0
    gamma = 0
! iterate to solve the non linear system
    do iter_loop = 1,NUM_ITER
! recompute jacobian for the new point
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
! compute distance to target location
      dx = - (x - x_target(irec))
      dy = - (y - y_target(irec))
      dz = - (z - z_target(irec))
! compute increments
      dxi    = xix   *dx + xiy   *dy + xiz   *dz
      deta   = etax  *dx + etay  *dy + etaz  *dz
      dgamma = gammax*dx + gammay*dy + gammaz*dz
! update values
      xi    = xi    + dxi
      eta   = eta   + deta
      gamma = gamma + dgamma
! impose that we stay in that element
! (useful if user gives a receiver outside the mesh for instance)
! we can go slightly outside the [1,1] segment since with finite elements
! the polynomial solution is defined everywhere
! can be useful for convergence of iterative scheme with distorted elements
      if (xi    >  1.10d0) xi    =  1.10d0
      if (xi    < -1.10d0) xi    = -1.10d0
      if (eta   >  1.10d0) eta   =  1.10d0
      if (eta   < -1.10d0) eta   = -1.10d0
      if (gamma >  1.10d0) gamma =  1.10d0
      if (gamma < -1.10d0) gamma = -1.10d0
! end of non linear iterations
    enddo
! compute final coordinates of point found
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
         xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
! compute final distance between asked and found (converted to km)
         distance = dsqrt((x_target(irec)-x)**2&
                        +(y_target(irec)-y)**2&
                        +(z_target(irec)-z)**2)*R_EARTH/1000.d0

        if(      xi>=-1.and.   xi<=1&
        .and.  eta>=-1.and.  eta<=1&
        .and.gamma>=-1.and.gamma<=1)then
! store xi,eta and x,y,z of point found
         found(irec) = .true.
         xi_receiver(irec)    = xi
         eta_receiver(irec)   = eta
         gamma_receiver(irec) = gamma
         x_found(irec)        = x
         y_found(irec)        = y
         z_found(irec)        = z
         final_distance(irec) = distance
         ispec_selected_rec(irec) = ispec
        endif
        if((.not.found(irec)).and.(distance<final_distance(irec)))then
! store xi,eta and x,y,z of point found
         xi_receiver(irec)    = xi
         eta_receiver(irec)   = eta
         gamma_receiver(irec) = gamma
         x_found(irec)        = x
         y_found(irec)        = y
         z_found(irec)        = z
         final_distance(irec) = distance
         ispec_selected_rec(irec) = ispec
        endif

       endif
      enddo
!
    endif
enddo

! for MPI version, gather information from all the nodes
  ispec_selected_rec_all(:,:) = -1
  call GATHER_ALL_I(ispec_selected_rec,nrec,ispec_selected_rec_all,nrec, NPROCTOT)

  call GATHER_ALL_DP(xi_receiver   , nrec, xi_receiver_all   , nrec, NPROCTOT)
  call GATHER_ALL_DP(eta_receiver  , nrec, eta_receiver_all  , nrec, NPROCTOT)
  call GATHER_ALL_DP(gamma_receiver, nrec, gamma_receiver_all, nrec, NPROCTOT)
  call GATHER_ALL_DP(final_distance, nrec, final_distance_all, nrec, NPROCTOT)
  call GATHER_ALL_DP(x_found       , nrec, x_found_all       , nrec, NPROCTOT)
  call GATHER_ALL_DP(y_found       , nrec, y_found_all       , nrec, NPROCTOT)
  call GATHER_ALL_DP(z_found       , nrec, z_found_all       , nrec, NPROCTOT)

! this is executed by main process only
  if(myrank == 0) then

    ! check that the gather operation went well
    if(any(ispec_selected_rec_all(:,:) == -1)) call exit_MPI(myrank,'gather operation failed for mirror receivers')

    ! MPI loop on all the results to determine the best slice
    islice_selected_rec(:) = -1
    n_rec_in_elm = 0
    n_rec_out_elm = 0
    r_max = 0
    r_min = hugeval
    do irec = 1,nrec
      distmin = HUGEVAL
      inside = .false.
      do iprocloop = 0,NPROCTOT-1
        if(      xi_receiver_all(irec,iprocloop)>=-1.and.   xi_receiver_all(irec,iprocloop)<=1&
         .and.  eta_receiver_all(irec,iprocloop)>=-1.and.  eta_receiver_all(irec,iprocloop)<=1&
         .and.gamma_receiver_all(irec,iprocloop)>=-1.and.gamma_receiver_all(irec,iprocloop)<=1)then
          inside=.true.
          distmin = final_distance_all(irec,iprocloop)
          islice_selected_rec(irec) = iprocloop
          ispec_selected_rec(irec) = ispec_selected_rec_all(irec,iprocloop)
          xi_receiver(irec) = xi_receiver_all(irec,iprocloop)
          eta_receiver(irec) = eta_receiver_all(irec,iprocloop)
          gamma_receiver(irec) = gamma_receiver_all(irec,iprocloop)
          x_found(irec) = x_found_all(irec,iprocloop)
          y_found(irec) = y_found_all(irec,iprocloop)
          z_found(irec) = z_found_all(irec,iprocloop)
        endif
        if((.not.inside).and.(final_distance_all(irec,iprocloop) < distmin)) then
          distmin = final_distance_all(irec,iprocloop)
          islice_selected_rec(irec) = iprocloop
          ispec_selected_rec(irec) = ispec_selected_rec_all(irec,iprocloop)
          xi_receiver(irec) = xi_receiver_all(irec,iprocloop)
          eta_receiver(irec) = eta_receiver_all(irec,iprocloop)
          gamma_receiver(irec) = gamma_receiver_all(irec,iprocloop)
          x_found(irec) = x_found_all(irec,iprocloop)
          y_found(irec) = y_found_all(irec,iprocloop)
          z_found(irec) = z_found_all(irec,iprocloop)
        endif
      enddo
        if(inside)then
            n_rec_in_elm = n_rec_in_elm+1
        else
            n_rec_out_elm = n_rec_out_elm+1
            r_max = max(r_max,sqrt((x_found(irec))**2.+(y_found(irec))**2.+(z_found(irec))**2.))
            r_min = min(r_max,sqrt((x_found(irec))**2.+(y_found(irec))**2.+(z_found(irec))**2.))
        endif
      final_distance(irec) = distmin
!      write(IMAIN,*)irec,n_rec_in_elm,n_rec_out_elm,distmin,islice_selected_rec(irec),ispec_selected_rec(irec)
!      write(IMAIN,*)xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec)
      if(final_distance(irec) == HUGEVAL)then
        print*,irec,nrec,ispec
        call exit_MPI(myrank,'error locating mirror receiver')
      endif   
    enddo
!    
    final_distance_max = maxval(final_distance(:))

    write(IMAIN,*)
    write(IMAIN,*)n_rec_in_elm,'receivers were found to be inside an element'
    write(IMAIN,*)
    write(IMAIN,*)n_rec_out_elm,'receivers were found to be outside an element'
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of all the mirror receivers: ',sngl(final_distance_max),' km'

    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for mirror receivers detection in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of mirror receivers detection - done'
    write(IMAIN,*)


  endif    ! end of section executed by main process only

! main process broadcasts the results to all the slices
  call BCAST_ALL_SINGLEI(nrec)
  !
  call synchronize_all()
  !
  call BCAST_ALL_I (islice_selected_rec, nrec)
  call BCAST_ALL_I (ispec_selected_rec , nrec)
  call BCAST_ALL_DP(xi_receiver        , nrec)
  call BCAST_ALL_DP(eta_receiver       , nrec)
  call BCAST_ALL_DP(gamma_receiver     , nrec)


    call flush_IMAIN()  

  ! deallocate arrays
  deallocate(inside_all)
  deallocate(epidist)
  deallocate(ix_initial_guess)
  deallocate(iy_initial_guess)
  deallocate(iz_initial_guess)
  deallocate(x_target)
  deallocate(y_target)
  deallocate(z_target)
  deallocate(x_found)
  deallocate(y_found)
  deallocate(z_found)
  deallocate(final_distance)
  deallocate(ispec_selected_rec_all)
  deallocate(xi_receiver_all)
  deallocate(eta_receiver_all)
  deallocate(gamma_receiver_all)
  deallocate(x_found_all)
  deallocate(y_found_all)
  deallocate(z_found_all)
  deallocate(final_distance_all)

  end subroutine locate_receivers_mirror

   subroutine setup_sources_receivers_intp_mir(myrank, &
                      xigll,yigll,zigll, &
                      SIMULATION_TYPE,nrec,nrec_local, &
                      islice_selected_rec,number_receiver_global, &
                      xi_receiver,eta_receiver,gamma_receiver, &
                      hxir_store,hetar_store,hgammar_store)

  use constants

  implicit none

  integer myrank

  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll


  integer SIMULATION_TYPE

  integer nrec,nrec_local
  integer, dimension(nrec) :: islice_selected_rec
  integer, dimension(nrec_local) :: number_receiver_global
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver

  double precision, dimension(nrec_local,NGLLX) :: hxir_store
  double precision, dimension(nrec_local,NGLLY) :: hetar_store
  double precision, dimension(nrec_local,NGLLZ) :: hgammar_store

  ! local parameters
  integer :: isource,irec,irec_local
  double precision, dimension(NGLLX) :: hxir,hpxir
  double precision, dimension(NGLLY) :: hpetar,hetar
  double precision, dimension(NGLLZ) :: hgammar,hpgammar


  ! select local receivers

  ! define local to global receiver numbering mapping
  irec_local = 0
  if (SIMULATION_TYPE == 1) then
    do irec = 1,nrec
      if(myrank == islice_selected_rec(irec)) then
        irec_local = irec_local + 1
        number_receiver_global(irec_local) = irec
      endif
    enddo
  endif

  ! define and store Lagrange interpolators at all the receivers
  if (SIMULATION_TYPE == 1) then
    do irec_local = 1,nrec_local
      irec = number_receiver_global(irec_local)
      call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
      call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
      hxir_store(irec_local,:) = hxir(:)
      hetar_store(irec_local,:) = hetar(:)
      hgammar_store(irec_local,:) = hgammar(:)
    enddo
  endif

  end subroutine setup_sources_receivers_intp_mir

!========================================
recursive function bspln(i,k,x) result(b)
!========================================

implicit none

real :: x,b
integer :: k,i

b=0.d0
if(k+1>0)then
   if(k+1==1)then
      if(x>=i.and.x<i+1)b=1.d0   
   else
      b=(x-i)*bspln(i,k-1,x)/(k+1-1)+(i+k+1-x)*bspln(i+1,k-1,x)/(k+1-1)
      if(k==0)b=0
   endif
endif

!=================
end function bspln
!=================

!====================
subroutine bmn(b,m,n)
!====================
  !
  implicit none
  !
  integer :: m,n,i
  real(4) :: b((n+1)*m+1),bspln
  !
  do i = 1,(n+1)*m+1
     b(i) = bspln(0,n,real(real(i-1)/real((n+1)*m)*real(n+1)))
  enddo
  !
return
!=================
end subroutine bmn
!=================


!===========================================
subroutine fill_spmat_diag(m,n,nsp,nt,spmat)
!===========================================
!
implicit none
!
integer :: m,n,nsp,nt
real :: spmat(n+1,nsp),b((n+1)*m+1)
!
integer :: i,j,i1,i2,i3,i4
!
call bmn(b,m,n)
!
do j = 1,nsp
   do i = 1,n+1
!
      i1 = max(1+(i-1)*m,(n+1-j)*m+1)
      i2 = min(1+(n+1)*m,(n-j+1)*m+nt)
      i3 = max(1,(n+2-j-i)*m+1)
      i4 = min(1+(n+2-i)*m,(n+2-j-i)*m+nt)
!
      if(j+i-1<=nsp)then
         spmat(i,j) = dot_product(b(i1:i2),b(i3:i4))
      else
         spmat(i,j) = 0.d0
      endif
!
   enddo
enddo
!
return
!=============================
end subroutine fill_spmat_diag
!=============================

!==========================================
subroutine bchfac ( w, nbands, nrow, diag )
!==========================================

!*****************************************************************************80
!
!! BCHFAC constructs a Cholesky factorization of a matrix.
!
!  Discussion:
!
!    The factorization has the form
!
!      C = L * D * L'
!  
!    with L unit lower triangular and D diagonal, for a given matrix C of 
!    order NROW, where C is symmetric positive semidefinite and banded, 
!    having NBANDS diagonals at and below the main diagonal.
! 
!    Gauss elimination is used, adapted to the symmetry and bandedness of C.
! 
!    Near-zero pivots are handled in a special way.  The diagonal 
!    element C(N,N) = W(1,N) is saved initially in DIAG(N), all N. 
! 
!    At the N-th elimination step, the current pivot element, W(1,N), 
!    is compared with its original value, DIAG(N).  If, as the result 
!    of prior elimination steps, this element has been reduced by about 
!    a word length, that is, if W(1,N) + DIAG(N) <= DIAG(N), then the pivot 
!    is declared to be zero, and the entire N-th row is declared to
!    be linearly dependent on the preceding rows.  This has the effect 
!    of producing X(N) = 0 when solving C * X = B for X, regardless of B.
! 
!    Justification for this is as follows.  In contemplated applications 
!    of this program, the given equations are the normal equations for 
!    some least-squares approximation problem, DIAG(N) = C(N,N) gives 
!    the norm-square of the N-th basis function, and, at this point, 
!    W(1,N) contains the norm-square of the error in the least-squares 
!    approximation to the N-th basis function by linear combinations 
!    of the first N-1.  
!
!    Having W(1,N)+DIAG(N) <= DIAG(N) signifies that the N-th function 
!    is linearly dependent to machine accuracy on the first N-1 
!    functions, therefore can safely be left out from the basis of 
!    approximating functions.
!
!    The solution of a linear system C * X = B is effected by the 
!    succession of the following two calls:
! 
!      call bchfac ( w, nbands, nrow, diag )
!
!      call bchslv ( w, nbands, nrow, b, x )
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) W(NBANDS,NROW).
!    On input, W contains the NBANDS diagonals in its rows, 
!    with the main diagonal in row 1.  Precisely, W(I,J) 
!    contains C(I+J-1,J), I=1,...,NBANDS, J=1,...,NROW.
!    For example, the interesting entries of a seven diagonal
!    symmetric matrix C of order 9 would be stored in W as
! 
!      11 22 33 44 55 66 77 88 99
!      21 32 43 54 65 76 87 98  *
!      31 42 53 64 75 86 97  *  *
!      41 52 63 74 85 96  *  *  *
!
!    Entries of the array not associated with an
!    entry of C are never referenced.
!    On output, W contains the Cholesky factorization 
!    C = L*D*L', with W(1,I) containing 1/D(I,I) and W(I,J) 
!    containing L(I-1+J,J), I=2,...,NBANDS.
!
!    Input, integer ( kind = 4 ) NBANDS, indicates the bandwidth of the
!    matrix C, that is, C(I,J) = 0 for NBANDS < abs(I-J).
! 
!    Input, integer ( kind = 4 ) NROW, is the order of the matrix C.
! 
!    Work array, real ( kind = 8 ) DIAG(NROW).
!
  implicit none

  integer nbands
  integer nrow

  real diag(nrow)
  integer i
  integer imax
  integer j
  integer jmax
  integer n
  real ratio
  real w(nbands,nrow)

  if ( nrow <= 1 ) then
    if ( 0.0D+00 < w(1,1) ) then
      w(1,1) = 1.0D+00 / w(1,1)
    end if
    return
  end if
!
!  Store the diagonal.
!
  diag(1:nrow) = w(1,1:nrow)
!
!  Factorization.
!
  do n = 1, nrow
 
    if ( w(1,n) + diag(n) <= diag(n) ) then
      w(1:nbands,n) = 0.0D+00
    else
 
      w(1,n) = 1.0D+00 / w(1,n)
 
      imax = min ( nbands - 1, nrow - n )
 
      jmax = imax
 
      do i = 1, imax
 
        ratio = w(i+1,n) * w(1,n)
 
        do j = 1, jmax
          w(j,n+i) = w(j,n+i) - w(j+i,n) * ratio
        end do
 
        jmax = jmax-1
        w(i+1,n) = ratio
 
      end do
 
    end if
 
  end do
 
  return

!==
end
!==


!==========================================================
subroutine bchslv (tmp, n1, n2, n3, w, nbands, nrow, lunit)
!==========================================================

!*****************************************************************************
!
!! BCHSLV solves a banded symmetric positive definite system.
!
!  Discussion:
!
!    The system is of the form:
!
!      C * X = B 
!  
!    and the Cholesky factorization of C has been constructed 
!    by BCHFAC.
! 
!    With the factorization 
!
!      C = L * D * L'
!
!    available, where L is unit lower triangular and D is diagonal, 
!    the triangular system 
!
!      L * Y = B 
!
!    is solved for Y (forward substitution), Y is stored in B, the 
!    vector D**(-1)*Y is computed and stored in B, then the 
!    triangular system L'*X = D**(-1)*Y is solved for X 
!    (back substitution).
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W(NBANDS,NROW), the Cholesky factorization for C, 
!    as computed by BCHFAC.
! 
!    Input, integer ( kind = 4 ) NBANDS, the bandwidth of C.
!
!    Input, integer ( kind = 4 ) NROW, the order of the matrix C.
! 
!    Input/output, real ( kind = 8 ) B(NROW).
!    On input, the right hand side.
!    On output, the solution.
!
  
  implicit none

  integer :: j
  integer :: n,n1,n2,n3
  integer :: nrow,nbands,lunit
  real :: w(nbands,nrow)
  character :: opt*5      
  real :: tmp(n1,n2,0:n3)   
! 
!nbands = mir%spln_order+1, nrow =  mir%count_rec
!
!  allocate(mir%tmp(0:mir%recl_mirror-1,0:mir%spln_order))
!
  if ( nrow <= 1 ) then
    write(*,*)'warning few time steps after decimation'
    return
  end if
!
!  Forward substitution. 
!  Solve L*Y = B.
!
  do n = 1, nrow-nbands+1
     if(n==1)then
        do j = 0, nbands - 1
           read(lunit,rec=j+n)tmp(:,:,j)
        enddo
     else
        do j = 0,nbands-2
           tmp(:,:,j) = tmp(:,:,j+1)
        enddo
         read(lunit,rec=n+nbands-1)tmp(:,:,nbands-1)
     endif
     do j = 1, nbands - 1
        tmp(:,:,j) = tmp(:,:,j) - w(j+1,n) * tmp(:,:,0)
     end do
     
     write(lunit,rec=n)tmp(:,:,0)
    
 end do

  do n = nrow-nbands+2, nrow

     do j = 0,nbands-2
        tmp(:,:,j) = tmp(:,:,j+1)
     enddo

     do j = 1, nrow - n
        tmp(:,:,j) = tmp(:,:,j) - w(j+1,n) * tmp(:,:,0)
     end do

      write(lunit,rec=n)tmp(:,:,0)
     
  end do
!
!  Back substitution. 
!  Solve L'*X = D**(-1)*Y.
!
  do n = nrow, nrow-nbands+2, -1

     do j = nrow-n,1,-1
        tmp(:,:,j) = tmp(:,:,j-1)
     enddo
     read(lunit,rec=n)tmp(:,:,0)
     
    tmp(:,:,0) = tmp(:,:,0) * w(1,n)

    do j = 1, nrow - n
      tmp(:,:,0) = tmp(:,:,0) - w(j+1,n) * tmp(:,:,j)
    end do

     write(lunit,rec=n)tmp(:,:,0)

  end do

  do n = nrow-nbands+1, 1, -1

     do j = nbands-1,1,-1
        tmp(:,:,j) = tmp(:,:,j-1)
     enddo
     read(lunit,rec=n)tmp(:,:,0)
     
    tmp(:,:,0) = tmp(:,:,0) * w(1,n)

    do j = 1, nbands - 1
      tmp(:,:,0) = tmp(:,:,0) - w(j+1,n) * tmp(:,:,j)
    end do

     write(lunit,rec=n)tmp(:,:,0)

  end do

!  deallocate(mir%tmp)
  return

!====================
end subroutine bchslv
!====================

