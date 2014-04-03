subroutine greens(H,G,wt,tot)
!--------------------------------------------------------------------------
! Calculates Green's function for the coherent potential approximation.
!--------------------------------------------------------------------------
! Variables:
! H(jsz,sec,sec) - Hamiltonian of the system
! wt(jsz) - Weight at the k-points
! tot - Total weight of all k-points
! G(sec,sec) - Green's function
! sig(nse) - Self-energies for disordered states
! f1000 - Formatting string for printing matrices verbosely
!--------------------------------------------------------------------------
use global
implicit none
integer(kind=4) :: i, l, l1
real(kind=8), intent(in) :: wt(jsz), tot
complex(kind=8), intent(out) :: G(sec,sec)
complex(kind=8), intent(inout) :: H(jsz,sec,sec)
character(len=100) :: f1000
if (verbose) print 1000
G(:,:) = (0.0d0,0.0d0)
write(f1000,'(A,I1,A)') "(",sec,"(2(F10.6,1x)))"
if (verbose.and.vlvl.ge.2) then
  print 1001
  print f1000, H(1,:,:)
end if
do i = 1, jsz
  call cmplxInv(H(i,:,:),sec,verbose,vlvl)
! this loop is only for Se/Te s & p disorder (currently)
  do l = 19, 22
    l1 = l + 9
    G(l,l) = G(l,l) + H(i,l,l)*wt(i)
    G(l1,l1) = G(l1,l1) + H(i,l1,l1)*wt(i)
  end do
end do
G(:,:) = G(:,:)/tot
if (verbose.and.vlvl.ge.2) print f1000, G
if (verbose) print 2000
return
1000 format(/,'Begin subroutine greens')
1001 format(/,'Hamiltonian at first k-point (usually Gamma)')
2000 format('End greens',/)
end subroutine greens
