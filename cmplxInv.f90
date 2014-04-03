subroutine cmplxInv(A,n,verbose,vlvl)
!--------------------------------------------------------------------------
! Solves for the inverse of a complex square matrix of size nxn. This is
! done by using Crout's LU decomposition with partial pivoting.
!--------------------------------------------------------------------------
! Variables:
! n - rank of matrix A
! A(n,n) - complex matrix to be inversed
! B(n,n) - Identity matrix
! p(n) - pivots for Crout's algorithm
! f**** - formatting strings
! verbose - Logical for debugging flags (.true. = debug info on)
! vlvl - Level of debugging verboseness
!--------------------------------------------------------------------------
implicit none
integer(kind=4) :: i, p(n)
integer(kind=4), intent(in):: n, vlvl
logical, intent(in) :: verbose
complex(kind=8) :: AA(n,n), B(n,n)
complex(kind=8), intent(inout) :: A(n,n)
character(len=100) :: f1000
if (verbose) print 1000
write(f1000,"(A,I1,A)") "(",n,"((2F10.5,1X)))"
AA = A; B = 0.0d0
do i = 1, n
  B(i,i) = cmplx(1.0d0,0.0d0,8)
end do
if (verbose.and.vlvl.ge.3) then
  print *, 'A (real + imaginary)'
  print f1000, transpose(AA)
  print *, 'B'
  print f1000, transpose(B)
end if
! Do LU decomposition using Crout's algorithm with partial pivoting on A
call ccrlu(AA,n,p)
if (verbose.and.vlvl.ge.3) then
  print *, 'LU (real + imaginary'
  print f1000, transpose(AA)
end if
! Now use the PLU decomposition to solve for the inverse of A by columns
!$OMP PARALLEL &
!$OMP SHARED(AA, n, p, b) &
!$OMP PRIVATE(i)
!$OMP DO
do i = 1, n
  call clubk(AA,n,p,b(i,:),verbose,vlvl)
  if (verbose.and.vlvl.ge.2) then
    print *, 'b'
    print 1001, b(i,:)
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
! Write matrix B (inverse of A) to matrix A
A = B
if (verbose.and.vlvl.ge.2) then
  print *, 'Inverse of A'
  print f1000, A
end if
if (verbose) print 2000
return
1000 format(/,'Begin subroutine cmplxInv')
1001 format(2F10.5)
1002 format(18X,1F15.10)
2000 format('End cmplxInv',/)
end subroutine cmplxInv
