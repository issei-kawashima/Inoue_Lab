call system_clock(t0)
open(999, file="LU_o.csv")
 do i=0,N
   write(999,*) LU(-1,i),",",LU(0,i),",",LU(1,i)
 end do
close(999)
call system_clock(t1,tr)
write(*,'(f10.3,A)')10000*(t1-t0)/dble(tr),"[s]"
stop
