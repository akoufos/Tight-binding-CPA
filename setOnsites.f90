subroutine setOnsites
!--------------------------------------------------------------------------
! Sets up the global onsite energies
!--------------------------------------------------------------------------
! Variables:
! ons_bar(sec) - Averaged onsite energies
! f1000 - Formatting string for printing matrices verbosely
!--------------------------------------------------------------------------
use global
use concentration
use onsites
implicit none
integer(kind=4) :: l, l1
character(len=100) :: f1000
if (verbose) print 1000
onsA(:,:) = (0.0d0,0.0d0)
onsB(:,:) = (0.0d0,0.0d0)
onsAvg(:,:) = (0.0d0,0.0d0)
ons_bar(:) = (con*ons(1,:) + (1.0d0-con)*ons(2,:))
do l = 1, sec
  if (l.ge.19.and.l.le.22) then
    l1 = l - 18
    onsA(l1,l1) = cmplx(ons(1,l),0.0d0,8)
    onsB(l1,l1) = cmplx(ons(2,l),0.0d0,8)
    onsAvg(l1,l1) = cmplx(ons_bar(l),0.0d0,8)
  else if (l.ge.28.and.l.le.31) then
    l1 = l - 23
    onsA(l1,l1) = cmplx(ons(1,l),0.0d0,8)
    onsB(l1,l1) = cmplx(ons(2,l),0.0d0,8)
    onsAvg(l1,l1) = cmplx(ons_bar(l),0.0d0,8)
  end if
end do
write(f1000,'(A,I1,A)') "(",nse,"(2(F10.6,1x)))"
if (verbose.and.vlvl.ge.2) then
  print 1001
  print f1000, onsA(:,:)
  print 1002
  print f1000, onsB(:,:)
  print 1003
  print f1000, onsAvg(:,:)
end if
if (verbose) print 2000
return
1000 format(/,'Begin subroutine setOnsites')
1001 format(/,'Onsite energies for system A that will be replaced')
1002 format(/,'Onsite energies for system B that will be repleaced')
1003 format(/,'Averaged onsite energies to be replaced by self-energies')
2000 format('End setOnsites',/)
end subroutine setOnsites
