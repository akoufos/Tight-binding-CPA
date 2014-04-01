subroutine trace(A, n, t, verbose, vlvl)
!--------------------------------------------------------------------------
! Finds the trace of a square complex matrix.
!--------------------------------------------------------------------------
! Variables:
! A(n,n) - Complex matrix in which the trace needs to be calculated
! n - Rank of matrix A
! t - The trace of the matrix A
!--------------------------------------------------------------------------
implicit none
integer(kind=4) :: i
integer(kind=4), intent(in) :: n, vlvl
real(kind=8) :: tr, ti
complex(kind=8), intent(in) :: A(n,n)
complex(kind=8), intent(out) :: t
logical :: verbose
character(len=100) :: f1000
if (verbose) print 1000
if (verbose.and.vlvl.gt.3) then
  write(f1000,"(A,I1,A)") "(",n,"(2(F10.5,1X)))"
  print *, 'A: '
  print f1000, transpose(A)
end if
do i = 1, n
  tr = tr + dble(A(i,i))
  ti = ti + aimag(A(i,i))
end do
t = cmplx(tr,ti,8)
if (verbose) print 2000
return
1000 format(/,'Begin subroutine trace')
2000 format('End trace',/)
end subroutine
