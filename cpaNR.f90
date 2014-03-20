subroutine cpaNR(ham,wt,tot,e,eps,dels1,delp1,mchk,numit)
!--------------------------------------------------------------------------
! Solves for the self-energies, and thus the Green's function, using the
! Newton-Raphson method.
!--------------------------------------------------------------------------
! Variables:
! grn(sec,sec) - Green's function
! sig(nse) - Self-energies of the disordered states
! ham(jsz,sec,sec) - Hamiltonian of the system
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
integer(4) :: i, n, irep
integer(4), intent(in) :: numit
integer(4), intent(out) :: mchk
real(8) :: cr, ci
real(8), intent(in) :: wt(jsz), tot, e, eps
complex(8) :: dels(2), delp(6), sig(nse), grn(sec,sec)
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
  call greens(ham,wt,tot,e,eps)
  call crete(dels,delp)
  if(verbose) print 1002, dels, delp
  if (abs(dble(dels(1))).le.cr.or.abs(aimag(dels(1))).le.ci.or. &
      abs(dble(dels(2))).le.cr.or.abs(aimag(dels(2))).le.ci.or. &
      abs(dble(delp(1))).le.cr.or.abs(aimag(delp(1))).le.ci.or. &
      abs(dble(delp(2))).le.cr.or.abs(aimag(delp(2))).le.ci.or. &
      abs(dble(delp(3))).le.cr.or.abs(aimag(delp(3))).le.ci.or. &
      abs(dble(delp(4))).le.cr.or.abs(aimag(delp(4))).le.ci.or. &
      abs(dble(delp(5))).le.cr.or.abs(aimag(delp(5))).le.ci.or. &
      abs(dble(delp(6))).le.cr.or.abs(aimag(delp(6))).le.ci) then
    mchk = 1
    if (verbose) print *,"Didn't you just converge???"
    do i = 1, nse
      if(aimag(sig(i)).gt.1.0d-20) then 
        irep = 1
        sig(i) = cmplx(dble(sig(i)),-aimag(sig(i)),8)
        if (verbose) print *,"Yeah, but ",aimag(sig(i))," of ",i," is &
          greater than zero"
      end if
    end do
  else
     sig(:) = sig(:) + dels1(1)
  end if
  if (irep.eq.1) cycle
  if (mchk.eq.1) exit
end do
if (verbose) print 1004
return
1000 format(2X,I5,8F14.9)
1001 format(7X,8F14.9)
1002 format(4(2(F12.8,1X)))
1003 format('SUBROUTINE CPANR',/)
1004 format('END SUBROUTINE CPANR',/)
end subroutine cpaNR
