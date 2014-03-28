subroutine readin(nd, ne)
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
implicit none
common /d1/ con
common /d2/ emin, emax, eps, del, epiv, sagsr1, sagsi1, sagpr1, sagpi1
common /d6/ ons
common /d7/ cr, ci
integer(4) :: i
integer(4), intent(out) :: nd
real(8) :: con, del, epiv, eps, es, ep, emin, emax, sagsr1, sagsi1, &
  sagpr1, sagpi1, ons(natom(2),sec), cr, ci
real(8), intent(out) :: ne
character(len=1) :: a(50)
if (verbose) print 1020
open(5,file='cpaper.dat',blank='zero')
read(5,1000) (a(i),i=1,50)
write(6,1000) (a(i),i=1,50)
read (5,*) nd,ne,cr,ci
write(6,1001) nd,ne,cr,ci
do i = 1, sec
  read(5,*)ons(1,i), ons(2,i)
end do
write(6,1007)
write(6,1008)ons(1,:)
write(6,1009)
write(6,1008)ons(2,:)
read(5,*) es,ep
write(6,1003) es,ep
read(5,*) con, emin, emax, del, epiv, sagsr1, sagsi1, sagpr1, sagpi1
write(6,1004) con, emin, emax, del, epiv, sagsr1, sagsi1, sagpr1, sagpi1
read(5,*) eps
close(5)
write(6,1005) eps, jsz
write(6,1006) cr, ci
if (verbose.and.vlvl.ge.2) then 
  print 1000, (a(i),i=1,50)
  print 1001, nd,ne,cr,ci
  print 1007
  print 1008, ons(1,:)
  print 1009
  print 1008, ons(2,:)
  print 1003, es,ep
  print 1004, con, emin, emax, del, epiv, sagsr1, sagsi1, sagpr1, sagpi1
  print 1005, eps, jsz
  print 1006, cr, ci
end if
if (con.lt.0.0d0.or.con.gt.1.0d0) then
  print 1010
  stop
end if
if (verbose) print 1021
return
1000 format(50A1)
1001 format(I5,5X,3F10.7)
1002 format(2F9.5)
1003 format(2F15.10)
1004 format(5F10.4/4F10.4,/)
1005 format(//'eps= ',F15.6,10X,'k-points: ',i5//)
1006 format(//,'Convergence criterion (real and imaginary)'/,2F15.9,//)
1007 format(//,'Onsite parameters of Selenium',/)
1008 format(4((9F12.8,1X),/))
1009 format(//,'Onsite parameters of Tellurium',/)
1010 format(//,'Concentration is not an element of [0:1]',/,'Check &
  cpaper.dat and rerun',//)
1020 format(/,'Begin subroutine readin')
1021 format('End readin',/)
end subroutine readin
