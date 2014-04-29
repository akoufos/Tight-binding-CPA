subroutine cpaDOS(dos,dos2,w,tot,e,eps)
!--------------------------------------------------------------------------
! Solves for the DOS of the system.
!--------------------------------------------------------------------------
! Variables:
! w(jsz) - Weights at each k-point
! tot - Total weight
! dos(sec+1) - Total and decomposed DOS of the system
! H(jsz,sec,sec) - Working Hamiltonian matrix
! f1000 - Formatting string for printing matrices verbosely
!--------------------------------------------------------------------------
use global
use hamiltonians
use onsites
use sigma
implicit none
integer(kind=4) :: i, l, l2
real(kind=8) :: gr(sec,sec), gi(sec,sec)
real(kind=8), intent(in) :: w(jsz), tot, e, eps
real(kind=8), intent(out) :: dos(sec+1), dos2(ntype,sec)
complex(kind=8) :: H(jsz,sec,sec)
character(len=100) :: f1000
if (verbose) print 1000
dos(:) = 0.0d0
dos2(:,:) = 0.0d0
H(:,:,:) = -ham(:,:,:)
call setHam(H,e,eps)
write(f1000,'(A,I1,A)') "(",sec,"(2(F10.6,1x)))"
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
    if (l.ge.19.and.l.le.22.or.l.ge.28.and.l.le.31) then
      if (l.ge.19.and.l.le.22) l2 = l - 18
      if (l.ge.28.and.l.le.31) l2 = l - 23
      dos2(1,l) = dos2(1,l) - pp*w(i)/tot*(gr(l,l)*aimag(sig(l2))+ &
        gi(l,l)*dble(sig(l2))-gi(l,l)*ons(2,l))/(ons(1,l)-ons(2,l))
      dos2(2,l) = dos2(2,l) - pp*w(i)/tot*(gr(l,l)*aimag(sig(l2))+ &
        gi(l,l)*dble(sig(l2))-gi(l,l)*ons(1,l))/(ons(2,l)-ons(1,l))
    end if
  end do
end do
dos(sec+1) = sum(dos(1:sec))
if (verbose.and.vlvl.ge.1) then
  print 1001, e, dos(:)
  print 1003, e, dos2(1,:)
  print 1004, e, dos2(2,:)
end if
print 1002, e, dos(sec+1)
if (verbose) print 2000
1000 format(/,'Start subroutine cpaDOS')
1001 format('DOS for energy ',F10.6,':',/,37F10.6)
1002 format('DOS for energy ',F10.6,':',5X,F10.6)
1003 format("Atom type 1's DOS for energy ",F10.6,':',/,37F10.6)
1004 format("Atom type 2's DOS for energy ",F10.6,':',/,37F10.6)
2000 format('End cpaDOS',/)
end subroutine cpaDOS
