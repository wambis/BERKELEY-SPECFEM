program convert_topo

integer, parameter :: nx = 5400
integer,parameter :: ny = 2700
integer :: mat(nx,ny)

open(1,file='DATA/topo_bathy/topo_bathy_etopo4_smoothed_window_7.dat',status='old')
open(2,file='new_topo.dat',status='unknown')
do i = 1,nx
      do j = 1,ny
        read(1,*)mat(i,j)
      enddo
      enddo

close(1)
close(2)
stop
end             
