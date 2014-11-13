module global
!--------------------------------------------------------------------------
!> Currently this global file is setup for calculations of P4/nmm Fe2Se/Te2
!--------------------------------------------------------------------------
! Variables:
!> @param jsz - Number of kpoints
!> @param mode - Used to decide how the program is run
!!  1 - Perform full CPA prgoram, including GG calculations
!!  2 - Perform only GG calculations. Must have run mode 1 or 3 previously
!!  3 - Perform VCA with DOS and GG calculations. Essentially skips 
!!      setting and calculating self-energies.
!> @param ntype - Number of different atom types
!> @param natom(ntype) - Number of different atoms (# of atoms in should 
!!   be given alphabetical order; i.e. UPd2Al3 should be (3,2,1))
!> @param nse - Number of self energies
!> @param pi Mathematical value of pi
!> @param sec - Number of secular equations
!> @param small Real value considered to be "small"
!> @param title - Title from cpaper.in to be used in all output files
!> @param verbose - Logical for debugging flags (.true. = debug info on)
!> @param vlvl - Level of debugging verboseness
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
integer(kind=4) :: mode
integer(kind=4), parameter :: jsz = 196 !904 ! From TB band calculation
integer(kind=4), parameter :: ntype = 2 ! Fe & Se
integer(kind=4), parameter :: natom(ntype) = (/ 2, 2 /) ! Fe = 2 , Se/Te = 2
integer(kind=4), parameter :: nse = 4*natom(2) ! s(1) + p(3) for each Se/Te
integer(kind=4), parameter :: sec = 9*2*2 !sum(natom) ! 9 (s(1) + p(3) + d(5))
logical :: verbose
real(kind=8), parameter :: pi = 4.0d0*datan(1.0d0), pp = 1.0d0/pi
real(kind=8), parameter :: small = 1.0d-20
character(len=75) :: title
save
end module global

!--------------------------------------------------------------------------
module converge
!--------------------------------------------------------------------------
!> Module for the convergence criterion related variables
!--------------------------------------------------------------------------
! Variables:
!> @param cr - Convergence criterion for the real part of the 
!!   Newton-Raphson procedure
!> @param ci - Convergence criterion for the imaginary part of the 
!!   Newton-Raphson procedure
!--------------------------------------------------------------------------
implicit none
real(kind=8) :: cr 
real(kind=8) :: ci
save
end module converge

!--------------------------------------------------------------------------
module concentration
!--------------------------------------------------------------------------
!> Module for the concentration related variables of the CPA program
!--------------------------------------------------------------------------
! Variables:
!> @param con - The concentration of the first atom type (e.g. 
!!   \f$FeSe_{con}Te_{1-con}\f$)
!--------------------------------------------------------------------------
implicit none
real(kind=8) :: con
save
end module

!--------------------------------------------------------------------------
module onsites
!--------------------------------------------------------------------------
!> Module for the onsites parameters and related variables
!--------------------------------------------------------------------------
! Variables:
!> @param ons(natom(2),sec) - The onsite parameters of the systems
!> @param ons_bar(sec) - Average onsite parameters between the two systems
!--------------------------------------------------------------------------
use global
implicit none
real(kind=8) :: ons(natom(2),sec)
real(kind=8) :: ons_bar(sec)
complex(kind=8) :: onsA(nse,nse)
complex(kind=8) :: onsB(nse,nse)
complex(kind=8) :: onsAvg(nse,nse)
save
end module onsites

!--------------------------------------------------------------------------
module hamiltonians
!--------------------------------------------------------------------------
!> Module for hamiltonians of the two systems and related variables
!--------------------------------------------------------------------------
! Variables:
!> @param hma(jsz,sec,sec) - Real part of initial Hamiltonian of system A 
!!   (e.g. FeSe)
!> @param vsa(jsz,sec,sec) - Imaginary part of initial Hamiltonian of 
!!   system A
!> @param hmb(jsz,sec,sec) - Real part of initial Hamiltonian of system B
!!   (e.g. FeTe)
!> @param vsb(jsz,sec,sec) - Imaginary part of initial Hamiltonian of 
!!   system B
!> @param ham(jsz,sec,sec) - Average Hamiltonian between the two systems
!--------------------------------------------------------------------------
use global
implicit none
real(kind=8) :: hma(jsz,sec,sec)
real(kind=8) :: hmb(jsz,sec,sec)
real(kind=8) :: vsa(jsz,sec,sec)
real(kind=8) :: vsb(jsz,sec,sec)
complex(kind=8) :: ham(jsz,sec,sec)
save
end module hamiltonians

!--------------------------------------------------------------------------
module green
!--------------------------------------------------------------------------
!> Module for greens function matrix and related variables
!--------------------------------------------------------------------------
! Variables:
!> @param grn(sec,sec) - Green's function
!--------------------------------------------------------------------------
use global
implicit none
complex(kind=8) :: grn(sec,sec)
save
end module green

!--------------------------------------------------------------------------
module sigma
!--------------------------------------------------------------------------
!> Module for the self-energies and related variables
!--------------------------------------------------------------------------
! Variables:
!> @param sig(nse) - Self energies of the system
!--------------------------------------------------------------------------
use global
implicit none
complex(kind=8) :: sig(nse)
save
end module sigma

!--------------------------------------------------------------------------
module unitconvert
!--------------------------------------------------------------------------
!> Module for all unit conversions
!--------------------------------------------------------------------------
! Variables:
!> @param ang2m - Angstroms to si units, meters
!> @param bohr2ang - Atomic radius, bohrs, to angstroms
!> @param bohr2m - Atomic radius, bohrs, to si units, meters
!> @param eV2J - Electron volts to si units, Joules (kgm^2/s^2 or Nm)
!> @param K2eV - Kelvin to electron volts
!> @param u2kg - Atomic mass to si units, kilograms
!> @param u2eVc2 - Atomic mass to electron volts per speed of light squared
!--------------------------------------------------------------------------
use global
implicit none
real(kind=8), parameter :: ang2m = 1.0d-10
real(kind=8), parameter :: bohr2ang = 0.52917721092d0
real(kind=8), parameter :: bohr2m = 0.52917721092d-11
real(kind=8), parameter :: eV_AA2kg_ss = 16.0217656d0
real(kind=8), parameter :: eV2Hz = 2.417989348d14
real(kind=8), parameter :: eV2J = 1.602176565d-18
real(kind=8), parameter :: K2eV = 8.621738d-5
real(kind=8), parameter :: K2meV = 8.621738d-2
real(kind=8), parameter :: kg_ss2eV_AA = 6.2415094d-2
real(kind=8), parameter :: meV2Hz = 2.417989348d11
real(kind=8), parameter :: Ry2ev = 13.60569253d0
real(kind=8), parameter :: u2kg = 1.660538921d-27
real(kind=8), parameter :: u2eV_c2 = 931.494061d6
real(kind=8), parameter :: uKK2eV_AA = 1.776386273d-6
real(kind=8), parameter :: uKK2eV_AA2 = 4.504401531323d-8
save
end module unitconvert
