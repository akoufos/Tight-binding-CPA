subroutine cpaNR(ham,wt,tot,e,eps,dels1,delp1,mchk,numit)
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
! dels(2) - Change in the self-energies of the disordered s states
! delp(6) - Change in the self-energies of the disordered p states
! mchk - 
! numit - Maximum number of iterations for the procedure
!--------------------------------------------------------------------------
use global
implicit none
common /d3/ sig, grn
common /d5/ dels, delp
common /d7/ cr, ci
integer(4) :: i, j, n, irep
integer(4), intent(in) :: numit
integer(4), intent(out) :: mchk
real(8) :: cr, ci
real(8), intent(in) :: wt(jsz), tot, e, eps
complex(8) :: z(nse), dels(2), delp(6), sig(nse), grn(sec,sec)
complex(8), intent(in) :: dels1(2), delp1(6)
complex(8), intent(out) :: ham(jsz,sec,sec)
do n = 1, numit
  if (verbose) print 1003
  write(6,1000)n,(sig(i),i=1,4)
  write(6,1001)(sig(i),i=5,8)
  if (verbose) then
    write(*,1000)n,(sig(i),i=1,4)
    write(*,1001)(sig(i),i=5,8)
  end if
  mchk = 0; irep = 0
  grn(:,:) = cmplx(0.0d0,0.0d0,8)
  sig(1) = sig(1) + dels1(1)
  sig(5) = sig(5) + dels1(2)
  sig(2) = sig(2) + delp1(1)
  sig(3) = sig(3) + delp1(2)
  sig(4) = sig(4) + delp1(3)
  sig(6) = sig(6) + delp1(4)
  sig(7) = sig(7) + delp1(5)
  sig(8) = sig(8) + delp1(6)
verbose = .false.
  call greens(ham,wt,tot,e,eps)
verbose = .true.
  call crete(dels,delp)
  if(verbose) print 1002, dels, delp
  !If any dels or delp parts are <= to convergence criterion set mchk=1
  if (abs(dble(dels(1))).le.cr.or.abs(aimag(dels(1))).le.ci.or. &
      abs(dble(dels(2))).le.cr.or.abs(aimag(dels(2))).le.ci.or. &
      abs(dble(delp(1))).le.cr.or.abs(aimag(delp(1))).le.ci) mchk = 1
! if (any(abs(dble(dels)).le.cr).or.any(abs(aimag(dels)).le.ci).or. &
!     any(abs(dble(delp)).le.cr).or.any(abs(aimag(delp)).le.ci)) mchk = 1
  if (mchk.eq.1) then
    if(aimag(sig(1)).gt.0.0d0.or.aimag(sig(5)).gt.0.0d0) then
      irep = 1
      sig(1) = cmplx(dble(sig(1)),-aimag(sig(1)),8)
      sig(5) = cmplx(dble(sig(5)),-aimag(sig(5)),8)
    end if
    if(aimag(sig(2)).gt.0.0d0.or.aimag(sig(6)).gt.0.0d0) then
      irep = 1
      sig(2) = cmplx(dble(sig(2)),-aimag(sig(2)),8)
      sig(3) = cmplx(dble(sig(3)),-aimag(sig(3)),8)
      sig(4) = cmplx(dble(sig(4)),-aimag(sig(4)),8)
      sig(6) = cmplx(dble(sig(6)),-aimag(sig(6)),8)
      sig(7) = cmplx(dble(sig(7)),-aimag(sig(7)),8)
      sig(8) = cmplx(dble(sig(8)),-aimag(sig(8)),8)
    end if
  end if
  if (irep.eq.1) cycle
  if (mchk.eq.1) exit
end do
verbose = .false.
return
1000 format(2X,I5,8F14.9)
1001 format(7X,8F14.9)
1002 format(4(2(F10.6,1X)))
1003 format('SUBROUTINE CPANR',/)
1003 format('END SUBROUTINE CPANR',/)
end subroutine cpaNR
