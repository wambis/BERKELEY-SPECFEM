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

!----
!---- assemble the contributions between slices and chunks using MPI
!----

! non-blocking routines

  subroutine assemble_MPI_vector_s(NPROC,nglob, &
                                           array_val, &
                                           buffer_send_vector,buffer_recv_vector, &
                                           num_interfaces,max_nibool_interfaces, &
                                           nibool_interfaces,ibool_interfaces, &
                                           my_neighbours, &
                                           request_send_vector,request_recv_vector)

! sends data

  use constants

  implicit none

  integer :: NPROC
  integer :: nglob

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: array_val

  integer :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces) :: &
       buffer_send_vector,buffer_recv_vector

  integer, dimension(num_interfaces) :: nibool_interfaces,my_neighbours
  integer, dimension(max_nibool_interfaces,num_interfaces) :: ibool_interfaces
  integer, dimension(num_interfaces) :: request_send_vector,request_recv_vector

  integer ipoin,iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then

     ! partition border copy into the buffer
     do iinterface = 1, num_interfaces
        do ipoin = 1, nibool_interfaces(iinterface)
           buffer_send_vector(:,ipoin,iinterface) = &
                array_val(:,ibool_interfaces(ipoin,iinterface))
        enddo
     enddo

     ! send messages
     do iinterface = 1, num_interfaces
        call isend_cr(buffer_send_vector(1,1,iinterface), &
                      NDIM*nibool_interfaces(iinterface), &
                      my_neighbours(iinterface), &
                      itag, &
                      request_send_vector(iinterface))

        call irecv_cr(buffer_recv_vector(1,1,iinterface), &
                      NDIM*nibool_interfaces(iinterface), &
                      my_neighbours(iinterface), &
                      itag, &
                      request_recv_vector(iinterface))
     enddo

  endif

  end subroutine assemble_MPI_vector_s

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_w(NPROC,nglob, &
                                           array_val, &
                                           buffer_recv_vector, &
                                           num_interfaces,max_nibool_interfaces, &
                                           nibool_interfaces,ibool_interfaces, &
                                           request_send_vector,request_recv_vector)

! waits for data to receive and assembles

  use constants

  implicit none

  integer :: NPROC
  integer :: nglob

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: array_val

  integer :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces) :: &
       buffer_recv_vector

  integer, dimension(num_interfaces) :: nibool_interfaces
  integer, dimension(max_nibool_interfaces,num_interfaces) :: ibool_interfaces
  integer, dimension(num_interfaces) :: request_send_vector,request_recv_vector

  integer ipoin,iinterface

! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! adding contributions of neighbours
    do iinterface = 1, num_interfaces
      do ipoin = 1, nibool_interfaces(iinterface)
        array_val(:,ibool_interfaces(ipoin,iinterface)) = &
             array_val(:,ibool_interfaces(ipoin,iinterface)) &
             + buffer_recv_vector(:,ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces
      call wait_req(request_send_vector(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_w


!-------------------------------------------------------------------------------------------------
!
! CUDA routines
!
!-------------------------------------------------------------------------------------------------


! interrupt might improve MPI performance
! see: https://computing.llnl.gov/tutorials/mpi_performance/#Sender-ReceiverSync
!
! check: MP_CSS_INTERRUPT environment variable on IBM systems


  subroutine assemble_MPI_vector_send_cuda(NPROC, &
                                          buffer_send_vector,buffer_recv_vector, &
                                          num_interfaces,max_nibool_interfaces, &
                                          nibool_interfaces, &
                                          my_neighbours, &
                                          request_send_vector,request_recv_vector)

  ! sends data
  ! note: array to assemble already filled into buffer_send_vector array
  use constants

  implicit none

  integer :: NPROC

  integer :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces) :: &
       buffer_send_vector,buffer_recv_vector

  integer, dimension(num_interfaces) :: nibool_interfaces,my_neighbours
  integer, dimension(num_interfaces) :: request_send_vector,request_recv_vector

  ! local parameters
  integer :: iinterface

  ! note: preparation of the contribution between partitions using MPI
  !          already done in transfer_boun_from_device() routines

  ! send only if more than one partition
  if(NPROC > 1) then

    ! send messages
    do iinterface = 1, num_interfaces
      call isend_cr(buffer_send_vector(1,1,iinterface), &
                    NDIM*nibool_interfaces(iinterface), &
                    my_neighbours(iinterface), &
                    itag, &
                    request_send_vector(iinterface))

      call irecv_cr(buffer_recv_vector(1,1,iinterface), &
                    NDIM*nibool_interfaces(iinterface), &
                    my_neighbours(iinterface), &
                    itag, &
                    request_recv_vector(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_send_cuda

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_write_cuda(Mesh_pointer,NPROC, &
                                            buffer_recv_vector, &
                                            num_interfaces,max_nibool_interfaces, &
                                            request_send_vector,request_recv_vector, &
                                            IREGION,FORWARD_OR_ADJOINT )

! waits for data to receive and assembles

  use constants

  implicit none

  integer(kind=8) :: Mesh_pointer

  integer :: NPROC

  integer :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces) :: buffer_recv_vector
  integer, dimension(num_interfaces) :: request_send_vector,request_recv_vector

  integer :: IREGION
  integer :: FORWARD_OR_ADJOINT

  ! local parameters

  integer :: iinterface

! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! adding contributions of neighbours
    call transfer_asmbl_accel_to_device(Mesh_pointer, &
                                        buffer_recv_vector, &
                                        IREGION,FORWARD_OR_ADJOINT)

    ! This step is done via previous function transfer_and_assemble...
    ! do iinterface = 1, num_interfaces
    !   do ipoin = 1, nibool_interfaces(iinterface)
    !     array_val(:,ibool_interfaces(ipoin,iinterface)) = &
    !          array_val(:,ibool_interfaces(ipoin,iinterface)) + buffer_recv_vector(:,ipoin,iinterface)
    !   enddo
    ! enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces
      call wait_req(request_send_vector(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_write_cuda


!
!-------------------------------------------------------------------------------------------------
!

! with cuda functions...

  subroutine transfer_boundary_to_device(Mesh_pointer, NPROC, &
                                            buffer_recv_vector, &
                                            num_interfaces,max_nibool_interfaces,&
                                            request_recv_vector, &
                                            IREGION,FORWARD_OR_ADJOINT)

  use constants

  implicit none

  integer(kind=8) :: Mesh_pointer

  integer :: NPROC

  ! array to assemble
  integer :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces) :: &
       buffer_recv_vector

  integer, dimension(num_interfaces) :: request_recv_vector

  integer :: IREGION
  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then

    ! waits for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! sends contributions to GPU
    call transfer_buffer_to_device_async(Mesh_pointer, &
                                         buffer_recv_vector, &
                                         IREGION,FORWARD_OR_ADJOINT)
  endif

  ! This step is done via previous function transfer_and_assemble...
  ! do iinterface = 1, num_interfaces
  !   do ipoin = 1, nibool_interfaces(iinterface)
  !     array_val(:,ibool_interfaces(ipoin,iinterface)) = &
  !          array_val(:,ibool_interfaces(ipoin,iinterface)) + buffer_recv_vector(:,ipoin,iinterface)
  !   enddo
  ! enddo

  end subroutine transfer_boundary_to_device

