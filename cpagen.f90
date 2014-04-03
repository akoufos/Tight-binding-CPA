program cpagen
!--------------------------------------------------------------------------
! This program is for tetragoanl pervoskites with disorder in the s & p
! orbitals of the second atom i.e. FeSe/Te, where Se/Te is the second atom.
!--------------------------------------------------------------------------
use global
implicit none
vlvl = 0
verbose = .false.
call mainn()
write(*,1000)
stop
1000 format(/,'Program finished running. Check cpaper.out for final &
  results',/)
end program cpagen
