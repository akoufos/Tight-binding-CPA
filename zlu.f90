subroutine zlu(n,m,a)
use global
implicit none
integer(4) :: i, j, k, p
integer(4), intent(in) :: m, n ! row, column
real(8), allocatable :: piv(:,:) ! pivot matrix
complex(8), allocatable :: mul(:,:) ! multiplier matrix
complex(8),intent(inout) :: a(n,m)
if(n.le.0.or.m.le.0) then
  write(6,85)
  write(*,85)
  write(*,*)'n: ',n,'m: ',m
  return
end if
allocate(piv(n,m),mul(n,m))
piv(:,:) = zero
! LU factorization
do k = 1, n - 1
  p = maxloc(real(a(k,:)))
  select case(p)
    case(p.eq.k)
      continue
    case(p.ne.k)
      piv(k,p) = 1.0d0
  end select
  if(real(a(k,k).eq.zero)) cycle
  do i = k+1, n
    real(mul(k,i)) = real(a(k,i))/real(a(k,k))
    if(aimag(a(k,k).eq.zero) then
      aimag(mul(k,i)) = zero
    else
      aimag(mul(k,i)) = aimag(a(k,i))/aimag(a(k,k))
    end if
  end do
  do j = k+1, n
    do i = k+1, m
      a(j,i) = a(j,i) - mul(j,i)*a(j,i)
    end do
  end do
end do
80 format('system near singular --- calculation continuing')
85 format('dimensioning problem in zlu')
return
end subroutine zlu
