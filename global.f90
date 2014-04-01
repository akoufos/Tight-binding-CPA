module global
!--------------------------------------------------------------------------
! Currently this global file is setup for calculations of P4/nmm Fe2Se/Te2
!--------------------------------------------------------------------------
! Variables:
! jsz - Number of kpoints
! ntype - Number of different atom types
! natom(ntype) - Number of different atoms (# of atoms in should 
  ! be given alphabetical order; i.e. UPd2Al3 should be (3,2,1))
! nse - Number of self energies
! sec - Number of secular equations
! verbose - Logical for debugging flags (.true. = debug info on)
! vlvl - Level of debugging verboseness
!--------------------------------------------------------------------------
! Notes for future release:
!
! Some of these parmeters should not be explicitly set in the global 
! module and just declared to be set from one of the reading subroutines.
! This will allow for more flexibility in the future and make it easier 
! to modify the code for more than just the FeSe/Te system.
!--------------------------------------------------------------------------
implicit none
integer(kind=4) :: vlvl
integer(kind=4), parameter :: jsz = 904 ! From TB band calculation
integer(kind=4), parameter :: ntype = 2 ! Fe & Se
integer(kind=4), parameter :: natom(ntype) = (/ 2, 2 /) ! Fe = 2 , Se/Te = 2
integer(kind=4), parameter :: nse = 4*natom(2) ! s(1) + p(3) for each Se/Te
integer(kind=4), parameter :: sec = 9*2*2 !sum(natom) ! 9 (s(1) + p(3) + d(5))
logical :: verbose
real(kind=8), parameter :: pi = 4.0d0*datan(1.0d0), pp = 1.0d0/pi
real(kind=8), parameter :: small = 1.0d-20
save
end module global

!--------------------------------------------------------------------------
module converge
!--------------------------------------------------------------------------
! Module for the convergence criterion related variables
!--------------------------------------------------------------------------
! Variables:
! cr - Convergence criterion for the real part of the Newton-Raphson 
  ! procedure
! ci - Convergence criterion for the imaginary part of the Newton-Raphson 
  ! procedure
!--------------------------------------------------------------------------
implicit none
real(kind=8) :: cr 
real(kind=8) :: ci
save
end module converge

!--------------------------------------------------------------------------
module concentration
!--------------------------------------------------------------------------
! Module for the concentration related variables of the CPA program
!--------------------------------------------------------------------------
! Variables:
! con - The concentration of the first atom type (e.g. FeSe_[con]Te-[1-con]
!--------------------------------------------------------------------------
implicit none
real(kind=8) :: con
save
end module

!--------------------------------------------------------------------------
module onsites
!--------------------------------------------------------------------------
! Module for the onsites parameters and related variables
!--------------------------------------------------------------------------
! Variables:
! ons(natom(2),sec) - The onsite parameters of the systems
! ons_bar(sec) - Average onsite parameters between the two systems
!--------------------------------------------------------------------------
use global
implicit none
real(kind=8) :: ons(natom(2),sec)
real(kind=8) :: ons_bar(sec)
save
end module onsites

!--------------------------------------------------------------------------
module hamiltonians
!--------------------------------------------------------------------------
! Module for hamiltonians of the two systems and related variables
!--------------------------------------------------------------------------
! Variables:
! hma(jsz,sec,sec) - Real part of initial Hamiltonian of system A (e.g.
  ! FeSe)
! vsa(jsz,sec,sec) - Imaginary part of initial Hamiltonian of system A
! hmb(jsz,sec,sec) - Real part of initial Hamiltonian of system B (e.g.
  ! FeTe)
! vsb(jsz,sec,sec) - Imaginary part of initial Hamiltonian of system B
!--------------------------------------------------------------------------
use global
implicit none
real(kind=8) :: hma(jsz,sec,sec)
real(kind=8) :: hmb(jsz,sec,sec)
real(kind=8) :: vsa(jsz,sec,sec)
real(kind=8) :: vsb(jsz,sec,sec)
save
end module hamiltonians
