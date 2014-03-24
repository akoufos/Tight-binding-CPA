subroutine crete(dels,delp)
!--------------------------------------------------------------------------
! Sets up the matrices to be inverted for the elements that have a self-
! energy
!--------------------------------------------------------------------------
! Variables:
! dels(2) - Change in the s states
! delp(6) - Change in the p states
! ge(nse,nse) - Diagonally reduced Green's matrix
! usi(nse,nse) - Diagonally reduced self-energy matrix
! up(nse,nse) - Matrix of ge*usi
! su(nse,nse) - Matrix to be inverted
! term(nse,nse) - Final matrix with concentrations. Results are the changes
  ! in states (del**)
! ons(2,sec) - All onsite parameters for Se and Te
! Se/Te(nse,nse) - Complex versions of onsite matices for calculations
!--------------------------------------------------------------------------
use global
implicit none
common /d1/ con
common /d3/ sig, grn
common /d6/ ons
integer(4) :: l, l1, l2, l3
real(8) :: con, ons(natom(2),sec)
complex(8) :: up(nse,nse), ge(nse,nse), su(nse,nse), term(nse,nse), &
  usi(nse,nse), grn(sec,sec), sig(nse), Se(nse,nse), Te(nse,nse), &
  eps(nse,nse)
complex(8), intent(out) :: dels(2), delp(6)
if (verbose) print 1000
dels(:) = (0.0d0,0.0d0)
delp(:) = (0.0d0,0.0d0)
usi(:,:) = (0.0d0,0.0d0)
ge(:,:) = (0.0d0,0.0d0)
su(:,:) = (0.0d0,0.0d0)
Se(:,:) = (0.0d0,0.0d0)
Te(:,:) = (0.0d0,0.0d0)
eps(:,:) = (0.0d0,0.0d0)
do l = 19, 22
  l1 = l + 9; l2 = l - 18; l3 = l - 14
  Se(l2,l2) = cmplx(ons(1,l),0.0d0,8)
  Te(l2,l2) = cmplx(ons(2,l),0.0d0,8)
  Se(l3,l3) = cmplx(ons(1,l1),0.0d0,8)
  Te(l3,l3) = cmplx(ons(2,l1),0.0d0,8)
end do
eps(:,:) = con*(Se(:,:)-Te(:,:)) + Te(:,:)
if(verbose.and.vlvl.ge.1) print 1001, sig
ge(1,1) = grn(19,19)
ge(2,2) = grn(20,20)
ge(3,3) = grn(21,21)
ge(4,4) = grn(22,22)
ge(5,5) = grn(28,28)
ge(6,6) = grn(29,29)
ge(7,7) = grn(30,30)
ge(8,8) = grn(31,31)
usi(1,1) = sig(1)
usi(2,2) = sig(2)
usi(3,3) = sig(3)
usi(4,4) = sig(4)
usi(5,5) = sig(5)
usi(6,6) = sig(6)
usi(7,7) = sig(7)
usi(8,8) = sig(8)
Se(:,:) = Se(:,:) - usi(:,:)
Te(:,:) = Te(:,:) - usi(:,:)
up = matmul(ge,Se)
su = matmul(up,Te)
term(:,:) = eps(:,:) - su(:,:)
dels(1) = sig(1) - term(1,1)
dels(2) = sig(5) - term(5,5)
delp(1) = sig(2) - term(2,2)
delp(2) = sig(3) - term(3,3)
delp(3) = sig(4) - term(4,4)
delp(4) = sig(6) - term(6,6)
delp(5) = sig(7) - term(7,7)
delp(6) = sig(8) - term(8,8)
if(verbose.and.vlvl.ge.1) print 1002, dels, delp
do l = 1, nse
  sig(l) = term(l,l)
end do
if(verbose.and.vlvl.ge.1) print 1001, sig
if(verbose) print 1003
return
1000 format(/,'Begin subroutine crete')
1001 format('Sigma: ',4(2(F10.6,1X)))
1002 format('Delta sig:',4(2(F10.6,1X)))
1003 format('End crete',/)
end subroutine crete
