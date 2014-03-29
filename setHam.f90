subroutine setHam(ham,e,eps)
!--------------------------------------------------------------------------
! Sets up the hamiltonian to be used in the greens and cpaDOS
! subroutines
!--------------------------------------------------------------------------
! Variables:
! ham(jsz,sec,sec) - Hamiltonian of the system
! ons_bar(sec) - Averaged onsite energies
! hma(jsz,sec,sec) - Real part of initial Hamiltonian of system A (e.g.
  ! FeSe)
! vsa(jsz,sec,sec) - Imaginary part of initial Hamiltonian of system A
! hmb(jsz,sec,sec) - Real part of initial Hamiltonian of system B (e.g.
  ! FeTe)
! vsb(jsz,sec,sec) - Imaginary part of initial Hamiltonian of system B
! sig(nse) - Self-energies for disordered states
! grn(sec,sec) - Green's function
! f1000 - Formatting string for printing matrices verbosely
!--------------------------------------------------------------------------
use global
implicit none
common /d1/ con
common /d3/ sig, grn
common /d4/ hma, vsa, hmb, vsb
common /d6/ ons
integer(4) :: i, l, l1, kpts
real(8) :: hma(jsz,sec,sec), vsa(jsz,sec,sec), hmb(jsz,sec,sec), &
  vsb(jsz,sec,sec), ons(natom(2),sec), ons_bar(sec), sig(nse), &
  grn(sec,sec), con
real(8), intent(in) :: e, eps
complex(8) :: H1(jsz,sec,sec), H2(jsz,sec,sec)
complex(8), intent(out) :: ham(jsz,sec,sec)
character(len=100) :: f1000
if (verbose) print 1000
kpts = 1
H1(:,:,:) = cmplx(-hma(:,:,:),-vsa(:,:,:),8)
H2(:,:,:) = cmplx(-hmb(:,:,:),-vsb(:,:,:),8)
ham(:,:,:) = (con*H1(:,:,:) + (1.0d0-con)*H2(:,:,:))
ons_bar(:) = (con*ons(1,:) + (1.0d0-con)*ons(2,:))
do l = 1, sec
  ham(:,l,l) = ham(:,l,l) + cmplx(e,eps,8)
  if (l.ge.19.and.l.le.22) then
    l1 = l + 9
    ham(:,l,l) = ham(:,l,l) - sig(l-18) + ons_bar(l)
    ham(:,l1,l1) = ham(:,l1,l1) - sig(l1-23) + ons_bar(l1)
  end if
end do
write(f1000,'(A,I1,A)') "(",sec,"(2(F10.6,1x)))"
if (verbose.and.vlvl.ge.2) then
  print 1005
  print f1000, ons_bar(:)
  if (vlvl.ge.3) kpts = jsz
  do i = 1, kpts
    print 1002
    print f1000, transpose(H1(i,:,:))
    print 1003
    print f1000, transpose(H2(i,:,:))
    print 1004
    print f1000, transpose(ham(i,:,:))
  end do
end if
if (verbose) print 2000
return
1000 format(/,'Begin subroutine setHam')
1002 format(/,'Hamiltonian for system A')
1003 format(/,'Hamiltonian for system B')
1004 format(/,'Averaged Hamiltonian with complex energy added and &
  onsite energies replaced with self-energies')
1005 format(/,'Averaged onsite parameters')
2000 format('End setHam',/)
end subroutine setHam
