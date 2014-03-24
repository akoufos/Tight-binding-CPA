module global
!--------------------------------------------------------------------------
! Currently this global file is setup for calculations of P4/nmm Fe2Se/Te2
!-------------------------------------------------------------------------
! Variables:
! jsz - Number of kpoints
! ntype - Number of different atom types
! natom(ntype) - Number of different atoms (# of atoms in should 
  ! be given alphabetical order; i.e. UPd2Al3 should be (3,2,1))
! nse - Number of self energies
! sec - Number of secular equations
! verbose - Logical for debugging flags (.true. = debug info on)
! vlvl - Level of debugging verboseness
!-------------------------------------------------------------------------
! Notes for future release:
!
! Some of these parmeters should not be explicitly set in the global 
! module and just declared to be set from one of the reading subroutines.
! This will allow for more flexibility in the future and make it easier 
! to modify the code for more than just the FeSe/Te system.
!--------------------------------------------------------------------------
implicit none
integer(4) :: vlvl
integer(4), parameter :: jsz = 75 ! From TB band calculation
integer(4), parameter :: ntype = 2 ! Fe & Se
integer(4), parameter :: natom(ntype) = (/ 2, 2 /) ! Fe = 2 , Se/Te = 2
integer(4), parameter :: nse = 4*natom(2) ! s(1) + p(3) for each Se/Te
integer(4), parameter :: sec = 9*2*2 !sum(natom) ! 9 (s(1) + p(3) + d(5))
logical :: verbose
real(8), parameter :: pi = 4.0d0*datan(1.0d0), pp = 1.0d0/pi
real(8), parameter :: small = 1.0d-20
end module global
