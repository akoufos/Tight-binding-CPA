subroutine greens(ham,wt,tot,e,eps)
!--------------------------------------------------------------------------
! Calculates Green's function for the coherent potential approximation.
!--------------------------------------------------------------------------
! Variables:
! ham(jsz,sec,sec) - Hamiltonian of the system
! ons_bar(sec) - Averaged onsite energies
! wt(jsz) - Weight at the k-points
! tot - Total weight of all k-points
! grn(sec,sec) - Green's function
! sig(nse) - Self-energies for disordered states
! f1000 - Formatting string for printing matrices verbosely
!--------------------------------------------------------------------------
use omp_lib
use global
implicit none
common /d3/ sig, grn
integer(4) :: i, l, l1
real(8) :: ons_bar(sec)
real(8), intent(in) :: wt(jsz), tot, e, eps
complex(8) :: sig(nse), grn(sec,sec)
complex(8), intent(out) :: ham(jsz,sec,sec)
character(len=100) :: f1000
if (verbose) print 1001
ham(:,:,:) = (0.0d0,0.0d0)
ons_bar(:) = 0.0d0
call setHam(ham,ons_bar,e,eps)
write(f1000,'(A,I1,A)') "(",sec,"(2(F10.6,1x)))"
if (verbose.and.vlvl.ge.2) print f1000, ham(1,:,:)
do i = 1, jsz
! this loop is only for Se/Te s & p disorder (currently)
  do l = 19, 22
    l1 = l + 9
    ham(i,l,l) = ham(i,l,l) - sig(l-18) + ons_bar(l)
    ham(i,l1,l1) = ham(i,l1,l1) - sig(l1-23) + ons_bar(l1)
  end do
  call cmplxInv(ham(i,:,:),sec,verbose,vlvl)
! this loop is only for Se/Te s & p disorder (currently)
  do l = 19, 22
    l1 = l + 9
    grn(l,l) = grn(l,l) + ham(i,l,l)*wt(i)
    grn(l1,l1) = grn(l1,l1) + ham(i,l1,l1)*wt(i)
  end do
end do
grn(:,:) = grn(:,:)/tot
if (verbose.and.vlvl.ge.2) print f1000, grn
if (verbose) print 1002
return
1001 format(/,'Begin subroutine greens')
1002 format('End greens',/)
end subroutine greens
