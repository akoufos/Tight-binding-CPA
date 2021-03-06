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
character(len=500) :: modestr
if (verbose) print 1000
open(5,file='cpaper.dat',blank='zero')
read(5,1001) title
write(7,1001) title 
read(5,1013) mode
read (5,*) nd,ne,cr,ci
select case (mode)
  case (1)
    write(modestr,'(A)')'Program is running in mode 1. This will run the &
      CPA program in its entirity and calculate the superconductivity &
      properties based on the Gaspari-Gyorffy and McMillan theories &
      using an approximation of the Hopfield parameter using the &
      calculated density of states.'
  case (2)
    write(modestr,'(A)')'Program is running in mode 2. This will only &
      run the cpaGG subroutine to calculate the superconductivity &
      properties based on the Gaspari-Gyorffy and McMillan theories &
      using an approximation of the Hopfield parameter using the &
      calculated density of states.'
  case (3)
    write(modestr,'(A)')'Program is running in mode 3. This will run a &
      virtual crystal approximation (VCA) of the system with DOS and &
      Gaspari-Gyorffy calculations.'
  case default
    write(modestr,'(A)')'No mode was selected. Please include a valid &
      mode in your input file and run the code again.'
    print 1012, mode, modestr
    stop
end select
write(7,1012) mode, modestr
write(7,1002) nd,ne,cr,ci
do i = 1, sec
  read(5,*)ons(1,i), ons(2,i)
end do
write(7,1008)
write(7,1009)ons(1,:)
write(7,1010)
write(7,1009)ons(2,:)
read(5,*) es,ep
write(7,1004) es,ep
read(5,*) con, emin, emax, del, epiv
read(5,*) (sagr1(i), sagi1(i), i=1,nse)
write(7,1005) con, emin, emax, del, epiv, (sagr1(i), sagi1(i), i=1,nse)
read(5,*) eps
close(5)
write(7,1006) eps, jsz
write(7,1007) cr, ci
if (verbose.and.vlvl.ge.2) then 
  print 1001, title
  print 1012, mode, modestr
  print 1002, nd,ne,cr,ci
  print 1008
  print 1009, ons(1,:)
  print 1010
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
1001 format(A75)
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
1012 format(//,'Mode: ',I1,//,A,//)
1013 format(5X,I5)
2000 format('End readin',/)
end subroutine readin
