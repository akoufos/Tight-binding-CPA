subroutine cmplxLin(A, n, b, verbose, vlvl)
!--------------------------------------------------------------------------
! Solves a set of n linear complex equations (A*x=b). This is
! done by using Crout's LU decomposition with partial pivoting.
!--------------------------------------------------------------------------
! Variables:
! n - Rank of matrix A
! A(n,n) - Complex matrix to be inversed
! b(n) - Right hand side of equation
! p(n) - Pivots for Crout's algorithm
! f**** - Formatting strings
! verbose - Logical for debugging flags (.true. = debug info on)
! vlvl - Level of debugging verboseness
!--------------------------------------------------------------------------
implicit none
integer(4) :: p(n)
integer(4), intent(in):: n, vlvl
logical, intent(in) :: verbose
complex(8) :: AA(n,n) 
complex(8), intent(inout) :: A(n,n), b(n)
character(len=100) :: f1000
if (verbose) print 1000
write(f1000,"(A,I1,A)") "(",n,"((2F10.5,1X)))"
AA = A
if (verbose.and.vlvl.ge.3) then
  print *, 'A (real + imaginary)'
  print f1000, transpose(AA)
  print *, 'b (real + imaginary)'
  print 1001, b
end if
! Do LU decomposition using Crout's algorithm with partial pivoting on A
call ccrlu(AA,n,p)
if (verbose.and.vlvl.ge.2) then
  print *, 'LU (real + imaginary'
  print f1000, transpose(AA)
end if
! Now use the PLU decomposition to solve for the inverse of A by columns
call clubk(AA,n,p,b,verbose)
if (verbose.and.vlvl.ge.1) then
  print *, 'Solution of A*x=b'
  print 1001, b
end if
if (verbose) print 1002
return
1000 format(/,'Begin subroutine cmplxLin')
1001 format(2F10.5)
1002 format('End cmplxLin',/)
end subroutine cmplxLin
