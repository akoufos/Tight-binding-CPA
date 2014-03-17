module global
! Currently this global file is setup for calculations of P4/nmm Fe2Se/Te2
!-------------------------------------------------------------------------
! jsz - Number of kpoints
! ntype - Number of different atom types
! natom(ntype) - Number of different atoms (# of atoms in should 
  ! be given alphabetical order; i.e. UPd2Al3 should be (3,2,1))
! nse - Number of self energies
! sec - Number of secular equations
! verbose - Used for debug information to screen
!-------------------------------------------------------------------------
implicit none
integer(4), parameter :: jsz = 904 ! From TB band calculation
integer(4), parameter :: ntype = 2 ! Fe & Se
integer(4), parameter :: natom(ntype) = (/ 2, 2 /) ! Fe = 2 , Se/Te = 2
integer(4), parameter :: nse = 4*natom(2) ! s(1) + p(3) for each Se/Te
integer(4), parameter :: sec = 9*sum(natom) ! 9 (s(1) + p(3) + d(5))
logical :: verbose
real(8), parameter :: pi = 4.0d0*datan(1.0d0), pp = 1.0d0/pi
end module global
