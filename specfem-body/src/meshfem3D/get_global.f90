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

  subroutine get_global(nspec,xp,yp,zp,iglob,locval,ifseg,nglob,npointot)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by Paul F. Fischer

! leave sorting subroutines in same source file to allow for inlining

  use constants

  implicit none

  ! input parameters
  integer, intent(in) :: npointot,nspec

  double precision, dimension(npointot), intent(in) :: xp,yp,zp

  integer, dimension(npointot), intent(out) :: iglob,locval
  logical, dimension(npointot), intent(out) :: ifseg
  integer, intent(out) :: nglob

  ! local variables
  double precision, dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: ind,ninseg,iwork
  integer :: ispec,i,j,ier
  integer :: ieoff,ilocnum,nseg,ioff,iseg,ig

  ! dynamically allocate arrays
  allocate(ind(npointot), &
          ninseg(npointot), &
          iwork(npointot), &
          work(npointot), &
          stat=ier)
  if( ier /= 0 ) stop 'error allocating local array in get_global'

  ! establish initial pointers
  do ispec=1,nspec
    ieoff=NGLLX * NGLLY * NGLLZ * (ispec-1)
    do ilocnum=1,NGLLX * NGLLY * NGLLZ
      locval(ilocnum+ieoff)=ilocnum+ieoff
    enddo
  enddo

  ifseg(:) = .false.

  nseg = 1
  ifseg(1) = .true.
  ninseg(1) = npointot

  do j=1,NDIM
    ! sort within each segment
    ioff=1
    do iseg=1,nseg
      if(j == 1) then
        call rank(xp(ioff),ind,ninseg(iseg))
      else if(j == 2) then
        call rank(yp(ioff),ind,ninseg(iseg))
      else
        call rank(zp(ioff),ind,ninseg(iseg))
      endif

      call swap_all(locval(ioff),xp(ioff),yp(ioff),zp(ioff),iwork,work,ind,ninseg(iseg))

      ioff=ioff+ninseg(iseg)
    enddo

    ! check for jumps in current coordinate
    ! compare the coordinates of the points within a small tolerance
    if(j == 1) then
      do i=2,npointot
        if(dabs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
      enddo
    else if(j == 2) then
      do i=2,npointot
        if(dabs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
      enddo
    else
      do i=2,npointot
        if(dabs(zp(i)-zp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
      enddo
    endif

    ! count up number of different segments
    nseg=0
    do i=1,npointot
      if(ifseg(i)) then
        nseg=nseg+1
        ninseg(nseg)=1
      else
        ninseg(nseg)=ninseg(nseg)+1
      endif
    enddo
  enddo

  ! assign global node numbers (now sorted lexicographically)
  ig=0
  do i=1,npointot
    if(ifseg(i)) ig=ig+1
    iglob(locval(i))=ig
  enddo

  nglob=ig

  ! deallocate arrays
  deallocate(ind,ninseg,iwork,work)

  end subroutine get_global

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_global_indirect_addressing(nspec,nglob,ibool)

!
!- we can create a new indirect addressing to reduce cache misses
! (put into this subroutine but compiler keeps on complaining that it can't vectorize loops...)

  use constants

  implicit none

  integer,intent(in) :: nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! local parameters
  ! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber
  integer:: i,j,k,ispec,ier

  ! copies original array
  allocate(copy_ibool_ori(NGLLX,NGLLY,NGLLZ,nspec), &
          mask_ibool(nglob), &
          stat=ier)
  if(ier /= 0) stop 'error allocating local arrays in get_global_indirect_addressing'

  ! initializes arrays
  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:,:) = ibool(:,:,:,:)

  ! reduces misses
  inumber = 0
  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          if(mask_ibool(copy_ibool_ori(i,j,k,ispec)) == -1) then
            ! creates a new point
            inumber = inumber + 1
            ibool(i,j,k,ispec) = inumber
            mask_ibool(copy_ibool_ori(i,j,k,ispec)) = inumber
          else
            ! uses an existing point created previously
            ibool(i,j,k,ispec) = mask_ibool(copy_ibool_ori(i,j,k,ispec))
          endif
        enddo
      enddo
    enddo
  enddo

  ! cleanup
  deallocate(copy_ibool_ori,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(mask_ibool,stat=ier); if(ier /= 0) stop 'error in deallocate'

  end subroutine get_global_indirect_addressing

!
!-------------------------------------------------------------------------------------------------
!

! sorting routines put in same file to allow for inlining

  subroutine rank(A,IND,N)
!
! Use Heap Sort (Numerical Recipes)
!
  implicit none

  integer :: n
  double precision,dimension(n) :: A
  integer,dimension(n) :: IND

  ! local parameters
  integer :: i,j,l,ir,indx
  double precision :: q

  do j=1,n
   IND(j)=j
  enddo

  if (n == 1) return

  L = n/2 + 1
  ir = n

  do while( .true. )

    IF ( l > 1 ) THEN
      l = l-1
      indx = ind(l)
      q = a(indx)
    ELSE
      indx = ind(ir)
      q = a(indx)
      ind(ir) = ind(1)
      ir = ir-1

      ! checks exit criterion
      if (ir == 1) then
         ind(1) = indx
         return
      endif
    endif

    i = l
    j = l+l

    do while( J <= IR )
      IF ( J < IR ) THEN
        IF ( A(IND(j)) < A(IND(j+1)) ) j=j+1
      endif
      IF ( q < A(IND(j)) ) THEN
        IND(I) = IND(J)
        I = J
        J = J+J
      ELSE
        J = IR+1
      endif
    enddo

    IND(I)=INDX
  enddo

  end subroutine rank

!
!-------------------------------------------------------------------------------------------------
!

  subroutine swap_all(IA,A,B,C,IW,W,ind,n)
!
! swap arrays IA, A, B and C according to addressing in array IND
!
  implicit none

  integer :: n
  integer,dimension(n) :: IND
  integer,dimension(n) :: IA,IW
  double precision,dimension(n) :: A,B,C,W

  ! local parameter
  integer :: i

  IW(:) = IA(:)
  W(:) = A(:)

  do i=1,n
    IA(i) = IW(ind(i))
    A(i) = W(ind(i))
  enddo

  W(:) = B(:)

  do i=1,n
    B(i) = W(ind(i))
  enddo

  W(:) = C(:)

  do i=1,n
    C(i) = W(ind(i))
  enddo

  end subroutine swap_all

