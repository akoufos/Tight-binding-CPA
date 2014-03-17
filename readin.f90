subroutine readin(nd, ne, cr, ci)
!--------------------------------------------------------------------------
! Reads inputs from cpaper.dat and writes some information to output file.
!--------------------------------------------------------------------------
! Variables:
! ci/cr = convergence criterion for real and imaginary parts of Green's 
!     function, respectively
! del = energy shift for site ?
! es/ep = initial onsite s & p energies for Se ?
! nd = # of dimension (# of secular equations; probably should remove)
! ne = # of electrons
! ni = # of ? (probably should remove)
! ns = # of disordered sites (probably should remove)
!--------------------------------------------------------------------------
use global
implicit none
common /d1/ con
common /d2/ emin, emax, eps, del, epiv, sagsr1, sagsi1, sagpr1, sagpi1
integer(4) :: i
integer(4), intent(out) :: nd
real(8) :: con, del, epiv, eps, es, ep, emin, emax, sagsr1, sagsi1, &
  sagpr1, sagpi1
real(8), intent(out) :: cr, ci, ne
character(len=1) :: a(20)
open(5,file='cpaper.dat',blank='zero')
read(5,1000) (a(i),i=1,20)
write(6,1000) (a(i),i=1,20)
read (5,*) nd,ne,cr,ci
write(6,1001) nd,ne,cr,ci
read(5,*) es,ep
write(6,1003) es,ep
read(5,*) con, emin, emax, del, epiv, sagsr1, sagsi1, sagpr1, sagpi1
write(6,1004) con, emin, emax, del, epiv, sagsr1, sagsi1, sagpr1, sagpi1
write(6,*)''
read(5,*) eps
close(5)
write(6,1005) eps, jsz
write(6,1006) cr, ci
return
1000 format(20A1)
1001 format(I5,5X,3F10.7)
1002 format(2F9.5)
1003 format(2F15.10)
1004 format(5F10.4/4F10.4)
1005 format(//'eps= ',F15.6,10X,'k-points: ',i5//)
1006 format(//,'Convergence criterion (real and imaginary)'/,2F15.9,//)
end subroutine readin
