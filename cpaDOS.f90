subroutine cpaDOS(dos,w,tot,e,eps)
!--------------------------------------------------------------------------
! Solves for the DOS of the system.
!--------------------------------------------------------------------------
! Variables:
! w(jsz) - Weights at each k-point
! tot - Total weight
! dos(sec+1) - Total and decomposed DOS of the system
! H(jsz,sec,sec) - Working Hamiltonian matrix
! ons_bar(sec) - Averaged onsite energies
! f1000 - Formatting string for printing matrices verbosely
!--------------------------------------------------------------------------
use global
implicit none
integer(4) :: i, l
real(8) :: gr(sec,sec), gi(sec,sec)
real(8), intent(in) :: w(jsz), tot, e, eps
real(8), intent(out) :: dos(sec+1)
complex(8) :: H(jsz,sec,sec)
character(len=100) :: f1000
if (verbose) print 1000
dos(:) = 0.0d0
H(:,:,:) = (0.0d0,0.0d0)
call setHam(H,e,eps)
if (verbose.and.vlvl.ge.3) then
  write(f1000,'(A,I1,A)') "(",sec,"(2(F10.6,1x)))"
  print*, 'H:'
  print f1000, transpose(H(1,:,:))
end if
if (verbose.and.vlvl.ge.3) then
  print*, 'H:'
  print f1000, transpose(H(1,:,:))
end if
do i = 1, jsz
  call cmplxInv(H(i,:,:),sec,verbose,vlvl)
  gr(:,:) = dble(H(i,:,:))
  gi(:,:) = aimag(H(i,:,:))
  do l = 1, sec
    dos(l) = dos(l) - pp*w(i)/tot*gi(l,l)
  end do
end do
dos(sec+1) = sum(dos(1:sec))
verbose = .true.; vlvl = 1
if (verbose.and.vlvl.ge.1) then
  print 1001, e, dos(:)
  print 1001, e, dos(sec+1)
end if
if (verbose) print 2000
1000 format(/,'Start subroutine cpaDOS')
1001 format('DOS for energy ',F10.6,':',/,37F10.6)
2000 format('End cpaDOS',/)
end subroutine cpaDOS
