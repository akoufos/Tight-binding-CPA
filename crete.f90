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
!--------------------------------------------------------------------------
use global
implicit none
common /d1/ con
common /d3/ sig, grn
integer(4) :: l
real(8) :: con
complex(8) :: up(nse,nse), ge(nse,nse), su(nse,nse), term(nse,nse), &
  usi(nse,nse), grn(sec,sec), sig(nse)
complex(8), intent(out) :: dels(2), delp(6)
usi(:,:) = (0.0d0,0.0d0)
ge(:,:) = (0.0d0,0.0d0)
su(:,:) = (0.0d0,0.0d0)
dels(:) = (0.0d0,0.0d0)
delp(:) = (0.0d0,0.0d0)
do l = 1, nse
  su(l,l) = cmplx(1.0d0,0.0d0,8)
end do
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
!call herakl(ge,usi,su,ans)
up = matmul(ge,usi)
su(:,:) = su(:,:) - up(:,:)
call cmplxInv(su,nse,verbose)
term(:,:) = con*matmul(usi,su)
dels(1) = term(1,1)
dels(2) = term(5,5)
delp(1) = term(2,2)
delp(2) = term(3,3)
delp(3) = term(4,4)
delp(4) = term(6,6)
delp(5) = term(7,7)
delp(6) = term(8,8)
return
end subroutine crete
