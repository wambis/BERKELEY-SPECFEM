program read_patch

integer :: io
character*1000 :: line
character*1 :: c1

open(1,file='patch',status='old')
open(2,file='patch.txt',status='unknown')

do 
read(1,'(a1000)',iostat=io)line
if(io<0)exit
c1 = line
if(line(1:1)=='+')then
if(line(2:2)=='*')line(2:2) = '!'
write(2,'(a200)')line(2:)
endif
enddo

close(2)
close(1)

end               
