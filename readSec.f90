subroutine readSec(h,g)
!--------------------------------------------------------------------------
! Reads in the secular equations of the system given by the output of the
! Static code.
!--------------------------------------------------------------------------
! Variables:
! h(jsz,sec,sec) - Real part of the Hamiltonian for each k-point
! j(jsz,sec,sec) - Imaginary part of the Hamiltonian for each k-point
!--------------------------------------------------------------------------
use global
implicit none
integer(4) :: i, j, k
real(8), intent(out) :: h(jsz,sec,sec), g(jsz,sec,sec)
g(:,:,:) = 0.0d0; h(:,:,:) = 0.0d0
open(16,file='cpamat.dat',blank='zero')
do i = 1, jsz
  do j = 1, 4
    read(16,*)
  end do
  do j = 1, sec
    do k = 1, j
      read(16,1000) h(i,j,k), g(i,j,k)
    end do
  end do
!$OMP PARALLEL DO DEFAULT(SHARED)
  do j = 1, sec
    do k = 1, j
      g(i,j,k) = -g(i,k,j)
      h(i,j,k) = h(i,k,j)
    end do
  end do
end do
!$OMP END PARALLEL DO
close(16)
return
1000 format(18X,2F15.10)
end subroutine
