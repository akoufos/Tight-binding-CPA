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
!--------------------------------------------------------------------------
use global
implicit none
common /d3/ sig, grn
common /d4/ hmz, vst
common /d6/ ons
integer(4) :: i, l, l1
real(8) :: ons(natom(2),sec), gr(sec,sec), gi(sec,sec), &
  hmz(jsz,sec,sec), vst(jsz,sec,sec), gsig(sec,sec)
real(8), intent(in) :: w(jsz), tot, e, eps
real(8), intent(out) :: dos(sec+1)
complex(8) :: H(jsz,sec,sec), Se(sec,sec), Te(sec,sec), grn(sec,sec), &
  sig(nse), Sigma(sec,sec), ham(jsz,sec,sec), tr
if (verbose) print 1000
dos(:) = 0.0d0
Se(:,:) = (0.0d0,0.0d0)
Te(:,:) = (0.0d0,0.0d0)
Sigma(:,:) = (0.0d0,0.0d0)
ham(:,:,:) = cmplx(-hmz(:,:,:),-vst(:,:,:),8)
if (verbose.and.vlvl.ge.3) then
  print*, 'ham:'
  print 1003, transpose(ham(1,:,:))
end if
do l = 1, sec
  ham(:,l,l) = ham(:,l,l) + cmplx(e,eps)
end do
H(:,:,:) = ham(:,:,:)
! this loop is only for Se/Te s & p disorder
do l = 19, 22
  l1 = l + 9
  ham(:,l,l) = ham(:,l,l) - sig(l-18) + ons(1,l)
  ham(:,l1,l1) = ham(:,l1,l1) - sig(l1-23) + ons(1,l1)
  Sigma(l,l) = sig(l-18)
  Sigma(l1,l1) = sig(l1-23)
end do
Sigma(:,:) = Sigma(:,:) - Te(:,:)
!H(:,:,:) = ham(:,:,:)
if (verbose.and.vlvl.ge.3) then
  print*, 'H:'
  print 1003, transpose(H(1,:,:))
end if
do l = 1, sec
  Se(l,l) = cmplx(ons(1,l),0.0d0,8)
  Te(l,l) = cmplx(ons(2,l),0.0d0,8)
end do
do i = 1, jsz
  call cmplxInv(H(i,:,:),sec,verbose,vlvl)
  gr(:,:) = dble(H(i,:,:))
  gi(:,:) = aimag(H(i,:,:))
  gsig = matmul(gi,Sigma)
  do l = 1, sec
!    dos(l) = dos(l) - pp*w(i)/tot*(gr(l,l)+gsig(l,l))/(Te(l,l)-Se(l,l))
    dos(l) = dos(l) - pp*w(i)/tot*gi(l,l)
    dos(sec+1) = dos(sec+1) - gi(l,l)*w(i)*pp/tot
  end do
! call trace(H(i,:,:),sec,tr,verbose)
! dos(sec+1) = dos(sec+1) - pp*w(i)/tot*aimag(tr)
end do
!dos(sec+1) = sum(dos(1:sec))
if (verbose.and.vlvl.ge.1) then
  print 1001, dos(:)
  print 1001, sum(dos(:sec))
end if
if (verbose) print 1002
1000 format(/,'Start subroutine cpaDOS')
1001 format('DOS:',/,37F10.6)
1002 format('End cpaDOS',/)
1003 format(36(2(F10.6,1X)))
end subroutine cpaDOS
