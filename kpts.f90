subroutine kpts(kp,q,wt,swt)
use global
implicit none
integer(4) :: i, j
integer(4), intent(in) :: kp
real(8), intent(out) :: q(kp,3), swt, wt(kp)
swt = 0.0d0; q(:,:) = 0.0d0; wt(:) = 0.0d0
open(7,file='cpaweights.dat')
read(7,*)
do i = 1, kp
  read(7,1000)(q(i,j),j=1,3),wt(i)
  swt = swt + wt(i)
end do
close(7)
write(6,1001) swt
return
1000 format(3F10.6,1X,1F10.6)
1001 format(' sum of weights is',e16.7)
end subroutine kpts
