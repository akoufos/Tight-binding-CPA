subroutine ccrlu(A,n,p)
!--------------------------------------------------------------------------
! LU decompostion subroutine with pivots for square matrix using
! Crout's algorithm (from rosetta code)
!--------------------------------------------------------------------------
! Variables:
! n - Rank of matrix A
! A(n,n) - Complex matrix of size nxn to be decomposed
! p(n) - Pivot vector
!--------------------------------------------------------------------------
implicit none
integer(4) :: i, j, k, piv
integer(4), intent(in) :: n
integer(4), intent(out) :: p(n)
real(8), parameter :: tiny = 1.0d-15
complex(8),intent(inout) :: A(n,n)
p = (/ (i, i=1, n) /)
do k = 1, n-1
  piv = k - 1 + maxloc(abs(a(p(k:),k)),1)
  if (piv.ne.k) then
    p((/k, piv/)) = p((/piv, k/))
  end if
  a(p(k+1:),k) = a(p(k+1:),k)/a(p(k),k)
! if (a(p(k+1:),k).le.tiny) a(p(k+1:),k) = 0.0d0
  forall (j = k+1:n)
    a(p(k+1:),j) = a(p(k+1:),j)-a(p(k+1:),k)*a(p(k),j)
!   if (a(p(k+1:),j).le.tiny) a(p(k+1:),j) = 0.0d0
  end forall
end do
return
end subroutine ccrlu
