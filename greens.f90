subroutine greens(ham,wt,tot,e,eps)
!--------------------------------------------------------------------------
! Calculates Green's function for the coherent potential approximation.
!--------------------------------------------------------------------------
! Variables:
! ham(jsz,sec,sec) - Hamiltonian of the system
! hmz(jsz,sec,sec) - Real part of initial Hamiltonian of the system
! vst(jsz,sec,sec) - Imaginary part of initial Hamiltonian of the system
! wt(jsz) - Weight at the k-points
! tot - Total weight of all k-points
! grn(sec,sec) - Green's function
! sig(nse) - Self-energies for disordered states
!--------------------------------------------------------------------------
use omp_lib
use global
implicit none
common /d3/ sig, grn
common /d4/ hmz, vst
common /d6/ ons
integer(4) :: i, l, l1
real(8) :: hmz(jsz,sec,sec), vst(jsz,sec,sec), ons(natom(2),sec)
real(8), intent(in) :: wt(jsz), tot, e, eps
complex(8) :: sig(nse), grn(sec,sec)
complex(8), intent(out) :: ham(jsz,sec,sec)
if (verbose) print 1001
ham(:,:,:) = cmplx(-hmz(:,:,:),-vst(:,:,:),8)
do l = 1, sec
  ham(:,l,l) = ham(:,l,l) + cmplx(e,eps)
end do
if (verbose.and.vlvl.ge.2) print 1000, ham(1,:,:)
do i = 1, jsz
! this loop is only for Se/Te s & p disorder (currently)
  do l = 19, 22
    l1 = l + 9
    ham(i,l,l) = ham(i,l,l) - sig(l-18) + ons(1,l)
    ham(i,l1,l1) = ham(i,l1,l1) - sig(l1-23) + ons(1,l1)
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
if (verbose.and.vlvl.ge.1) print 1000, grn
if (verbose) print 1002
return
1000 format(36(2(F10.6,1X)))
1001 format(/,'Begin subroutine greens')
1002 format('End greens',/)
end subroutine greens
