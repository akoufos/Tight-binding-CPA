subroutine cmplxInv(A,n,verbose)
! Solves for the inverse of a complex square matrix of size nxn. This is
! done by using Crout's LU decomposition with partial pivoting.
! n - rank of matrix A
! A(n,n) - complex matrix to be inversed
! B(n,n) - Identity matrix
! p(n) - pivots for Crout's algorithm
! f**** - formatting strings
implicit none
integer(4) :: i, p(n)
integer(4), intent(in):: n
logical, intent(in) :: verbose
complex(8) :: AA(n,n), B(n,n)
complex(8), intent(inout) :: A(n,n)
character(len=100) :: f1000
write(f1000,"(A,I1,A)") "(",n,"((2F10.5,1X)))"
AA = A; B = 0.0d0
do i = 1, n
  B(i,i) = cmplx(1.0d0,0.0d0,8)
end do
if (verbose) then
  print *, 'A (real + imaginary)'
  print f1000, transpose(AA)
  print *, 'B'
  print f1000, transpose(B)
end if
! Do LU decomposition using Crout's algorithm with partial pivoting on A
call ccrlu(AA,n,p)
if (verbose) then
  print *, 'LU (real + imaginary'
  print f1000, transpose(AA)
end if
! Now use the PLU decomposition to solve for the inverse of A by columns
do i = 1, n
  call clubk(AA,n,p,b(i,:),verbose)
  if (verbose) then
    print *, 'b'
    print 1001, b(i,:)
  end if
end do
! Write matrix B (inverse of A) to matrix A
A = B
if (verbose) then
  print *, 'Inverse of A'
  print f1000, A
end if
return
1001 format(2F10.5)
1002 format(18X,1F15.10)
end subroutine cmplxInv
