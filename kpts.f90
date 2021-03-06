subroutine kpts(kp,q,wt,swt)
!--------------------------------------------------------------------------
! Reads in the k-point information; Direction vector, weight at given 
! k-point and total weight at all k-points
!--------------------------------------------------------------------------
! Variables:
! kp - Number of k-points
! q(kp,3) - k-point vector
! wt(kp) - Weights for each k-point
! swt - Total weight
!--------------------------------------------------------------------------
use global
implicit none
integer(kind=4) :: i, j
integer(kind=4), intent(in) :: kp
real(kind=8), intent(out) :: q(kp,3), swt, wt(kp)
if (verbose) print 1000
swt = 0.0d0; q(:,:) = 0.0d0; wt(:) = 0.0d0
open(10,file='cpaweights.dat')
read(10,*)
do i = 1, kp
  read(10,1001)(q(i,j),j=1,3),wt(i)
  swt = swt + wt(i)
end do
if (verbose.and.vlvl.ge.1) then
  do i = 1, kp
    print 1001, (q(i,j),j=1,3),wt(i)
  end do
end if
close(10)
write(7,1002) swt
if (verbose.and.vlvl.ge.1) print 1002, swt
if (verbose) print 2000
return
1000 format(/,'Begin subroutine kpts')
1001 format(3F10.6,1X,1F10.6)
1002 format(' sum of weights is',e16.7)
2000 format('End kpts',/)
end subroutine kpts
