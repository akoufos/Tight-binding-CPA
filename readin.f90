subroutine readin(nd, ne, sagr1, sagi1)
!--------------------------------------------------------------------------
! Reads inputs from cpaper.dat and writes some information to output file.
!--------------------------------------------------------------------------
! Variables:
! ci/cr - Convergence criterion for real and imaginary parts of Green's 
!     function, respectively
! del - Energy shift for site ?
! es/ep - Initial onsite s & p energies for Se ?
! nd - # of dimension (# of secular equations; probably should remove)
! ne - Number of electrons
! ons(natom(2),nse) - Onsite parameters of Se and Te respectively
!--------------------------------------------------------------------------
use global
use converge
use concentration
use onsites
implicit none
common /d2/ emin, emax, eps, del, epiv
integer(kind=4) :: i
integer(kind=4), intent(out) :: nd
real(kind=8) :: del, epiv, eps, es, ep, emin, emax
real(kind=8), intent(out) :: ne, sagi1(nse), sagr1(nse)
character(len=1) :: a(50)
if (verbose) print 1000
open(5,file='cpaper.dat',blank='zero')
read(5,1001) (a(i),i=1,50)
write(6,1001) (a(i),i=1,50)
read (5,*) nd,ne,cr,ci
write(6,1002) nd,ne,cr,ci
do i = 1, sec
  read(5,*)ons(1,i), ons(2,i)
end do
write(6,1008)
write(6,1009)ons(1,:)
write(6,1010)
write(6,1009)ons(2,:)
read(5,*) es,ep
write(6,1004) es,ep
read(5,*) con, emin, emax, del, epiv
read(5,*) (sagr1(i), sagi1(i), i=1,nse)
write(6,1005) con, emin, emax, del, epiv, (sagr1(i), sagi1(i), i=1,nse)
read(5,*) eps
close(5)
write(6,1006) eps, jsz
write(6,1007) cr, ci
if (verbose.and.vlvl.ge.2) then 
  print 1001, (a(i),i=1,50)
  print 1002, nd,ne,cr,ci
  print 1008
  print 1009, ons(1,:)
  print 1007
  print 1009, ons(2,:)
  print 1004, es,ep
  print 1005, con, emin, emax, del, epiv, (sagr1(i), sagi1(i), i=1,nse)
  print 1006, eps, jsz
  print 1007, cr, ci
end if
if (con.lt.0.0d0.or.con.gt.1.0d0) then
  print 1011
  stop
end if
if (verbose) print 2000
return
1000 format(/,'Begin subroutine readin')
1001 format(50A1)
1002 format(I5,5X,3F10.7)
1003 format(2F9.5)
1004 format(2F15.10)
1005 format(5F10.4/2(4(2F10.4),/))
1006 format(//'eps= ',F15.6,10X,'k-points: ',i5//)
1007 format(//,'Convergence criterion (real and imaginary)'/,2F15.9,//)
1008 format(//,'Onsite parameters of Selenium',/)
1009 format(4((9F12.8,1X),/))
1010 format(//,'Onsite parameters of Tellurium',/)
1011 format(//,'Concentration is not an element of [0:1]',/,'Check &
  cpaper.dat and rerun',//)
2000 format('End readin',/)
end subroutine readin
