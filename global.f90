module global
implicit none
! jsz = # of kpoints
! nse = # of self energies [2*states (2 for real and imag)]
! sec = # of secular equations
integer(4), parameter :: jsz = 904, nse = 8, sec = 36
logical :: verbose
real(8), parameter :: pi = 4.0d0*datan(1.0d0), pp = 1.0d0/pi
end module global
