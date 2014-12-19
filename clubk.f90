subroutine clubk(A, n, piv, b, verbose, vlvl)
!--------------------------------------------------------------------------
!> Solves A*x = b with partial pivoting from Crout's LU 
!!   [\f$(P*L)*U = P*A\f$]
!--------------------------------------------------------------------------
! Variables:
!> @param n - Rank of matrix A
!> @param piv(n) - Pivot vector from subroutine ccrlu
!> @param A(n,n) - Complex matrix that is decomposed by subroutine ccrlu
!> @param L(n,n) - Lower decomposition of matrix a
!> @param U(n,n) - Upper decomposition of matrix a
!> @param P(n,n) - Pivot matrix
!> @param y(n) - (U*x = y)
!> @param x(n) - Working vector (b will equal x at the end)
!> @param b(n) - Input and solution of A*x = b [\f$ (P*L)*y = b \f$]
!> @param verbose - Logical for debugging flags (.true. = debug info on)
!> @param vlvl - Level of debugging verboseness
!--------------------------------------------------------------------------
implicit none
integer(kind=4) :: i, j, P(n,n)
integer(kind=4), intent(in) :: n, piv(n), vlvl
logical, intent(in) :: verbose
complex(kind=8) :: L(n,n), U(n,n), y(n), x(n)
complex(kind=8), intent(in) :: A(n,n)
complex(kind=8), intent(inout) :: b(n)
character(len=100) :: f1000
if (verbose) print 1000
P = 0; L = 0.0d0; U = 0.0d0; x = 0.0d0; y = 0.0d0
do i = 1, n
  P(i,i) = 1
  L(i,i) = 1.0d0
end do
do i = 1, n
  L(i,:i-1) = A(piv(i),:i-1)
end do
do i = 1, n
  U(i,i:) = A(piv(i),i:)
end do
P(:,piv) = P
b = matmul(dble(P),b)
if (verbose.and.vlvl.ge.4) then
  write(f1000,"(A,I1,A)") "(",n,"(2(F10.5,1X)))"
  print *, 'P'
  print f1000, transpose(dble(P))
  print *, 'L'
  print f1000, transpose(L)
  print *, 'U'
  print f1000, transpose(U)
  print *, 'b'
  print 1001, b
end if
y(1) = b(1)/L(1,1)
do i = 2, n
  y(i) = b(i)/L(i,i)
  do j = 1, i-1
    y(i) = y(i) - L(i,j)*y(j)/L(i,i)
  end do
end do
if (verbose.and.vlvl.ge.3) then
  print *, 'y'
  print 1001, y 
end if
x(n) = y(n)/U(n,n)
do i = n-1, 1, -1
  x(i) = y(i)/U(i,i)
  do j = i+1, n
    x(i) = x(i) - U(i,j)*x(j)/U(i,i)
  end do
end do
if (verbose.and.vlvl.ge.2) then
  print *, 'x'
  print 1001, x
end if
b = x
if (verbose) print 2000
1000 format(/,'Begin subroutine clubk')
1001 format(2F10.5)
2000 format('End clubk',/)
end subroutine clubk
