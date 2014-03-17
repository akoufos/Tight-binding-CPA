subroutine cpaDOS(dos,ham,w,tot)
!--------------------------------------------------------------------------
! Solves for the DOS of the system.
!--------------------------------------------------------------------------
! Variables:
! w(jsz) - Weights at each k-point
! tot - Total weight
! dos(sec+1) - Total and decomposed DOS of the system
! H(jsz,sec,sec) - Working Hamiltonian matrix
! ham(jsz,sec,sec) - Hamiltonian matrix
!--------------------------------------------------------------------------
use global
implicit none
integer(4) :: i, l
real(8), intent(in) :: w(jsz), tot
real(8), intent(out) :: dos(sec+1) ! total & decomposed DOS
complex(8) :: H(jsz,sec,sec)
complex(8), intent(in) :: ham(jsz,sec,sec)
dos(:) = 0.0d0; H(:,:,:) = ham(:,:,:)
do i = 1, jsz
  call cmplxInv(H(i,:,:),sec,verbose)
  do l = 1, sec
    dos(l) = dos(l) + aimag(H(i,l,l))*w(i)*pp/tot
  end do
end do
dos(sec+1) = sum(dos(1:sec))
end subroutine cpaDOS
