subroutine zcrlu(n,m,a,p,d)
! LU decompostion subroutine using the Crout method with row pivots
use global
implicit none
integer(4) :: i, j, k, piv
integer(4), intent(in) :: m, n, p(m) ! row, column, pivots
real(8) :: big, dum, ssum, temp
real(8), allocatable :: vv(:) ! scaling array
real(8), intent(in) :: d(m)
complex(8),intent(inout) :: a(n,m)
if(n.le.0.or.m.le.0) then
  write(6,1000)
  write(*,1000)
  write(*,*)'n: ',n,'m: ',m
  return
end if
allocate(vv(m))
d(:) = 1.0d0
! LU factorization
do i = 1, m ! loop over rows
  big = 0.0d0
  do j = 1, n
    temp = dsqrt(real(a(j,i))**2.0d0 + aimag(a(j,i))**2.0d0)
!   temp = abs(a(j,i))
    if(temp.gt.big) big = temp
  end do
  if(big.eq.0.0d0) write(*,1001)
  vv(i) = 1.0d0/big
end do
do j = 1, n ! loop over columns
  i = 1
  while (i.lt.j) do
    ssum = a(j,i)
    k = 1
    while (k.lt.i) do
      ssum = ssum - a(j,k)*a(k,i)
      k = k + 1
    end do
    a(j,i) = ssum
    i = i + 1
  end do
  big = 0.0d0
  do i = j, m
    ssum = a(j,i)
    do k = 1, j
      ssum = ssum - a(j,k)*a(k,i)
    end do
    a(j,i) = ssum
    dum = vv(i)*abs(ssum)
    if(dum.ge.big) then
      big = dum
      piv = i
    end if
  end do
  if(j.ne.piv) then
    do k = 1, m
      dum = a(k,piv)
      a(k,piv) = a(k,j)
      a(k,j) = dum
    end do
    d(:) = -d(:)
    vv(piv) = vv(j)
  end if
  p(j) = piv
  if(a(i,j.eq0.0d0) a(i,j) = cmplx(tiny,tiny,8)
  if(j.ne.m) then
    dum = 1.0d0/a(j,j)
    do i = j+1 m
      a(j,i) = a(j,i) * dum
    end do 
end do
deallocate(vv)
1000 format('dimensioning problem in zcrlu')
1001 format('Sigular matrix in routine zcrlu')
return
end subroutine zcrlu
