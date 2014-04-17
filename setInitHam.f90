subroutine setInitHam
!--------------------------------------------------------------------------
! Sets up the initial Hamiltonian to be modified at each specific energy
!--------------------------------------------------------------------------
! Variables:
! con - The concentration of the first atom type (e.g. FeSe_[con]Te-[1-con]
! ham(jsz,sec,sec) - Hamiltonian of the system
! H1(jsz,sec,sec) - Hamiltonian of system A
! H2(jsz,sec,sec) - Hamiltonian of system B
! f1000 - Formatting string for printing matrices verbosely
!--------------------------------------------------------------------------
use global
use concentration
use hamiltonians
implicit none
integer(kind=4) :: i, kpts
complex(kind=8) :: H1(jsz,sec,sec), H2(jsz,sec,sec)
character(len=100) :: f1000
if (verbose) print 1000
kpts = 1
H1(:,:,:) = cmplx(hma(:,:,:),vsa(:,:,:),8)
H2(:,:,:) = cmplx(hmb(:,:,:),vsb(:,:,:),8)
ham(:,:,:) = (con*H1(:,:,:) + (1.0d0-con)*H2(:,:,:))
write(f1000,'(A,I1,A)') "(",sec,"(2(F10.6,1x)))"
if (verbose.and.vlvl.ge.2) then
  if (vlvl.ge.3) kpts = jsz
  do i = 1, kpts
    print 1001
    print f1000, transpose(H1(i,:,:))
    print 1002
    print f1000, transpose(H2(i,:,:))
    print 1003
    print f1000, transpose(ham(i,:,:))
  end do
end if
if (verbose) print 2000
return
1000 format(/,'Begin subroutine setInitHam')
1001 format(/,'Hamiltonian for system A')
1002 format(/,'Hamiltonian for system B')
1003 format(/,'Averaged Hamiltonian before complex energy added and &
  onsite energies replaced with self-energies')
2000 format('End setInitHam',/)
end subroutine setInitHam
