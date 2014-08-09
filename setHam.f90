subroutine setHam(H,e,eps)
!--------------------------------------------------------------------------
! Sets up the Hiltonian to be used in the greens and cpaDOS
! subroutines
!--------------------------------------------------------------------------
! Variables:
! H(jsz,sec,sec) - Hamiltonian of the system at specific energy level
! ons_bar(sec) - Averaged onsite energies
! sig(nse) - Self-energies for disordered states
! f1000 - Formatting string for printing matrices verbosely
!--------------------------------------------------------------------------
use global
use concentration
use onsites, only : ons_bar
use sigma
implicit none
integer(kind=4) :: i, l, kpts
real(kind=8), intent(in) :: e, eps
complex(kind=8), intent(inout) :: H(jsz,sec,sec)
character(len=100) :: f1000
if (verbose) print 1000
do l = 1, sec
  H(:,l,l) = H(:,l,l) + cmplx(e,eps,8)
  if (mode.eq.1) then
    if (l.ge.19.and.l.le.22) then
      H(:,l,l) = H(:,l,l) - sig(l-18) + ons_bar(l)
    else if (l.ge.28.and.l.le.31) then
      H(:,l,l) = H(:,l,l) - sig(l-23) + ons_bar(l)
    end if
  end if
end do
write(f1000,'(A,I1,A)') "(",sec,"(2(F10.6,1x)))"
if (verbose.and.vlvl.ge.2) then
  if (vlvl.ge.3) kpts = jsz
  do i = 1, kpts
    print 1001
    print f1000, transpose(H(i,:,:))
  end do
end if
if (verbose) print 2000
return
1000 format(/,'Begin subroutine setHam')
1001 format(/,'Averaged Hamiltonian with complex energy added and &
  onsite energies replaced with self-energies')
2000 format('End setHam',/)
end subroutine setHam
