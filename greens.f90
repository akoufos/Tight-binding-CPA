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
integer(4) :: i, i1, i2, l
real(8) :: hmz(jsz,sec,sec), vst(jsz,sec,sec)
real(8), intent(in) :: wt(jsz), tot, e, eps
complex(8) :: sig(nse), grn(sec,sec)
complex(8), intent(out) :: ham(jsz,sec,sec)
ham(:,:,:) = cmplx(-hmz(:,:,:),-vst(:,:,:),8)
do l = 1, sec
  ham(:,l,l) = ham(:,l,l) + cmplx(e,eps)
end do
if (verbose) print 1000, ham(1,:,:)
do i = 1, jsz
! this loop is only for Se/Te s & p disorder (currently)
  do i1 = 19, 22
    i2 = i1 + 9
    ham(i,i1,i1) = ham(i,i1,i1) - sig(i1-18)
    ham(i,i2,i2) = ham(i,i2,i2) - sig(i2-23)
  end do
  call cmplxInv(ham(i,:,:),sec,verbose)
! this loop is only for Se/Te s & p disorder (currently)
  do i1 = 19, 22
    i2 = i1 + 9
    grn(i1,i1) = grn(i1,i1) + ham(i,i1,i1)*wt(i)
    grn(i2,i2) = grn(i2,i2) + ham(i,i2,i2)*wt(i)
  end do
end do
grn(:,:) = grn(:,:)/tot
if (verbose) print 1000, grn
return
1000 format(36(2(F10.6,1X)))
end subroutine greens
