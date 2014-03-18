subroutine cpaNR(ham,hmz,vst,wt,tot,e,eps,mchk,irep,numit)
!--------------------------------------------------------------------------
! Solves for the self-energies, and thus the Green's function, using the
! Newton-Raphson method.
!--------------------------------------------------------------------------
! Variables:
! grn(sec,sec) - Green's function
! sig(nse) - Self-energies of the disordered states
! ham(jsz,sec,sec) - Hamiltonian of the system
! hmz(jsz,sec,sec) - Real part of initial Hamiltonian of the system
! vst(jsz,sec,sec) - Imaginary part of initial Hamiltonian of the system
! wt(jsz) - Weights at each k-point in the Brillouin zone
! tot - Total weight from all k-points
! e - 
! eps - 
! mchk - 
! irep - 
! numit - Maximum number of iterations for the procedure
!--------------------------------------------------------------------------
use global
implicit none
common /d3/ sig, grn
integer(4) :: i, j, n
integer(4), intent(in) :: numit
integer(4), intent(out) :: mchk, irep
real(8), intent(in) :: hmz(jsz,sec,sec) vst(jsz,sec,sec), wt(jsz), tot, &
  e, eps
complex(8), intent(out) :: ham(jsz,sec,sec)
do n = 1, numit
  write(6,1000)n,sigs(1),(sigp(i),i=1,3)
  if (verbose) write(*,1000)n,sigs(1),(sigp(i),i=1,3)
  mchk = 0; irep = 0
  grn(:,:) = cmplx(0.0d0,0.0d0,8)
  sig(1) = sig(1) + cmplx(delsr1,delsi1,8)
  sig(5) = sig(5) + cmplx(delsr1,delsi1,8)
  sig(2) = sig(2) + cmplx(delpr1,delpi1,8)
  sig(3) = sig(3) + cmplx(delpr1,delpi1,8)
  sig(4) = sig(4) + cmplx(delpr1,delpi1,8)
  sig(6) = sig(6) + cmplx(delpr1,delpi1,8)
  sig(7) = sig(7) + cmplx(delpr1,delpi1,8)
  sig(8) = sig(8) + cmplx(delpr1,delpi1,8)
  call greens(ham,hmz,vst,weight,totvol,e,eps)
  call crete(dels,delp)
  if(verbose) print 1001, dels, delp
  
1000 format()
1001 format(4(2(F10.6,1X)))
end subroutine cpaNR
