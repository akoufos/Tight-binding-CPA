subroutine fixpt(G,sig,del)
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
complex(kind=8) :: ge(nse), term(nse), alpha(nse), beta(nse)
complex(kind=8), intent(in) :: G(sec,sec)
complex(kind=8), intent(out) :: del(nse)
complex(kind=8), intent(inout) :: sig(2,nse)
if (verbose) print 1000
del(:) = (0.0d0,0.0d0)
if(verbose.and.vlvl.ge.1) print 1001, sig(1,:)
if(verbose.and.vlvl.ge.1) print 1001, sig(2,:)
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
  beta(l) = onsAvg(l,l) - onsA(l,l)*onsB(l,l)*ge(l) - ge(l)*sig(1,l)**2.0d0
  term(l) = beta(l)/alpha(l)
  del(l) = sig(1,l) - term(l)
end do
if(verbose.and.vlvl.ge.1) then
  if(vlvl.ge.2) then
    print 1004
    print 1003, beta(:)
    print 1003, alpha(:)
    print 1003, term(:)
  end if
  print 1002, del(:)
end if
do l = 1, nse
  sig(1,l) = term(l)
end do
if(verbose.and.vlvl.ge.1) print 1001, sig(1,:)
if(verbose) print 2000
return
1000 format(/,'Begin subroutine fixpt')
1001 format('Sigma: ',/,4(2(F10.6,1X)))
1002 format('Delta sig:',/,4(2(E12.5,2X)))
1003 format(/,8(2F10.6,1x),/)
1004 format(//,'Numerator, Denominator, Self-energies')
2000 format('End fixpt',/)
end subroutine fixpt
