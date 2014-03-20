subroutine trace(A,n,t,verbose)
!--------------------------------------------------------------------------
! Finds the trace of a square complex matrix.
!--------------------------------------------------------------------------
! Variables:
! A(n,n) - Complex matrix in which the trace needs to be calculated
! n - Rank of matrix A
! t - The trace of the matrix A
!--------------------------------------------------------------------------
implicit none
integer(4) :: i
integer(4), intent(in) :: n
real(8) :: tr, ti
complex(8), intent(in) :: A(n,n)
complex(8), intent(out) :: t
logical :: verbose
character(len=100) :: f1000
if (verbose) then
  print 1000
  write(f1000,"(A,I1,A)") "(",n,"(2(F10.5,1X)))"
  print *, 'A: '
  print f1000, transpose(A)
end if
do i = 1, n
  tr = tr + dble(A(i,i))
  ti = ti + aimag(A(i,i))
end do
t = cmplx(tr,ti,8)
return
1000 format(/,'Begin subroutine trace')
1001 format('End trace',/)
end subroutine
