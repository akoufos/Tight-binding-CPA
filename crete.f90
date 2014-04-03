subroutine crete(G,dels,delp)
!--------------------------------------------------------------------------
! Sets up the matrices to be inverted for the elements that have a self-
! energy
!--------------------------------------------------------------------------
! Variables:
! dels(2) - Change in the s states
! delp(6) - Change in the p states
! G(sec,sec) - Green's function calculated from greens subroutine
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
use concentration
use onsites
use sigma
implicit none
integer(kind=4) :: l
complex(kind=8) :: up(nse,nse), ge(nse,nse), su(nse,nse), term(nse,nse), &
  usi(nse,nse), Se(nse,nse), Te(nse,nse)
complex(kind=8), intent(in) :: G(sec,sec)
complex(kind=8), intent(out) :: dels(2), delp(6)
if (verbose) print 1000
dels(:) = (0.0d0,0.0d0)
delp(:) = (0.0d0,0.0d0)
usi(:,:) = (0.0d0,0.0d0)
ge(:,:) = (0.0d0,0.0d0)
su(:,:) = (0.0d0,0.0d0)
if(verbose.and.vlvl.ge.1) print 1001, sig
ge(1,1) = G(19,19)
ge(2,2) = G(20,20)
ge(3,3) = G(21,21)
ge(4,4) = G(22,22)
ge(5,5) = G(28,28)
ge(6,6) = G(29,29)
ge(7,7) = G(30,30)
ge(8,8) = G(31,31)
usi(1,1) = sig(1)
usi(2,2) = sig(2)
usi(3,3) = sig(3)
usi(4,4) = sig(4)
usi(5,5) = sig(5)
usi(6,6) = sig(6)
usi(7,7) = sig(7)
usi(8,8) = sig(8)
Se(:,:) = onsA(:,:) - usi(:,:)
Te(:,:) = onsB(:,:) - usi(:,:)
up = matmul(ge,Se)
su = matmul(up,Te)
term(:,:) = onsAvg(:,:) - su(:,:)
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
if(verbose) print 2000
return
1000 format(/,'Begin subroutine crete')
1001 format('Sigma: ',4(2(F10.6,1X)))
1002 format('Delta sig:',4(2(F10.6,1X)))
2000 format('End crete',/)
end subroutine crete
