subroutine cpaDOS(dos,w,tot,e,eps)
!--------------------------------------------------------------------------
! Solves for the DOS of the system.
!--------------------------------------------------------------------------
! Variables:
! w(jsz) - Weights at each k-point
! tot - Total weight
! dos(sec+1) - Total and decomposed DOS of the system
! H(jsz,sec,sec) - Working Hamiltonian matrix
! ham(jsz,sec,sec) - Hamiltonian matrix
! ons_bar(sec) - Averaged onsite energies
! f1000 - Formatting string for printing matrices verbosely
!--------------------------------------------------------------------------
use global
implicit none
common /d1/ con
common /d3/ sig, grn
common /d4/ hma, vsa, hmb, vsb
integer(4) :: i, l, l1
real(8) :: gr(sec,sec), gi(sec,sec), ons_bar(sec)
real(8), intent(in) :: w(jsz), tot, e, eps
real(8), intent(out) :: dos(sec+1)
complex(8) :: H(jsz,sec,sec), grn(sec,sec), sig(nse)
character(len=100) :: f1000
if (verbose) print 1000
dos(:) = 0.0d0
H(:,:,:) = (0.0d0,0.0d0)
ons_bar(:) = 0.0d0
call setHam(H,ons_bar,e,eps)
if (verbose.and.vlvl.ge.3) then
  write(f1000,'(A,I1,A)') "(",sec,"(2(F10.6,1x)))"
  print*, 'ham:'
  print f1000, transpose(H(1,:,:))
end if
! this loop is only for Se/Te s & p disorder
do l = 19, 22
  l1 = l + 9
  ham(:,l,l) = ham(:,l,l) - sig(l-18) + ons_bar(l)
  ham(:,l1,l1) = ham(:,l1,l1) - sig(l1-23) + ons_bar(l1)
end do
if (verbose.and.vlvl.ge.3) then
  print*, 'H:'
  print f1000, transpose(H(1,:,:))
end if
do i = 1, jsz
  call cmplxInv(H(i,:,:),sec,verbose,vlvl)
  gr(:,:) = dble(H(i,:,:))
  gi(:,:) = aimag(H(i,:,:))
  do l = 1, sec
    dos(l) = dos(l) - pp*w(i)/tot*gi(l,l)
!    dos(sec+1) = dos(sec+1) - gi(l,l)*w(i)*pp/tot
  end do
end do
dos(sec+1) = sum(dos(1:sec))
if (verbose.and.vlvl.ge.1) then
  print 1001, dos(:)
  print 1001, sum(dos(:sec))
end if
if (verbose) print 1002
1000 format(/,'Start subroutine cpaDOS')
1001 format('DOS:',/,37F10.6)
1002 format('End cpaDOS',/)
end subroutine cpaDOS
