subroutine newton(G,sig,del)
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
implicit none
integer(kind=4) :: l
complex(kind=8) :: ge(nse), term(nse), F(nse), dF(nse), alpha(nse), &
  beta(nse)
complex(kind=8), intent(in) :: G(sec,sec)
complex(kind=8), intent(out) :: del(nse)
complex(kind=8), intent(inout) :: sig(nse)
if (verbose) print 1000
del(:) = (0.0d0,0.0d0)
!delp(:) = (0.0d0,0.0d0)
if(verbose.and.vlvl.ge.1) print 1001, sig
ge(1) = G(19,19)
ge(2) = G(20,20)
ge(3) = G(21,21)
ge(4) = G(22,22)
ge(5) = G(28,28)
ge(6) = G(29,29)
ge(7) = G(30,30)
ge(8) = G(31,31)
do l = 1, nse
  alpha(l) = 1.0d0 - onsB(l,l)*ge(l) - onsA(l,l)*ge(l)
  beta(l) = ge(l)*onsA(l,l)*onsB(l,l) - onsAvg(l,l)
  F(l) = ge(l)*sig(l)**2.0d0 + alpha(l)*sig(l) + beta(l)
  dF(l) = alpha(l) + 2.0d0*ge(l)*sig(l)
  term(l) = sig(l) - F(l)/dF(l)
  del(l) = sig(l) - term(l)
end do
!dels(1) = sig(1) - term(1)
!dels(2) = sig(5) - term(5)
!delp(1) = sig(2) - term(2)
!delp(2) = sig(3) - term(3)
!delp(3) = sig(4) - term(4)
!delp(4) = sig(6) - term(6)
!delp(5) = sig(7) - term(7)
!delp(6) = sig(8) - term(8)
!verbose = .true.; vlvl = 3
if(verbose.and.vlvl.ge.1) then
  if(vlvl.ge.2) then
    print 1004
    print 1003, F(:)
    print 1003, dF(:)
    print 1003, term(:)
    print 1003, -alpha(:)/(2.0d0*ge(:))
  end if
  print 1002, del
end if
!verbose = .false.; vlvl = 3
do l = 1, nse
  sig(l) = term(l)
end do
if(verbose.and.vlvl.ge.1) print 1001, sig
if(verbose) print 2000
return
1000 format(/,'Begin subroutine newton')
1001 format('Sigma: ',/,4(2(F10.6,1X)))
1002 format('Delta sig:',/,4(2(E12.5,2X)))
1003 format(/,8(2F10.6,1x),/)
1004 format(//,'F, dF, term')
2000 format('End newton',/)
end subroutine newton
