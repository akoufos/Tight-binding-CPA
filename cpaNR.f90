subroutine cpaNR(ham,wt,tot,e,eps,del1,mchk,numit)
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
! e - Current energy level
! eps - Imaginary energy temperature broadening value
! dels(2) - Change in the self-energies of the disordered s states
! delp(6) - Change in the self-energies of the disordered p states
! mchk - Flag to specify whether routine converged or not
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
complex(8), intent(in) :: del1
complex(8), intent(out) :: ham(jsz,sec,sec)
do n = 1, numit
  if (verbose) print 1003
  write(6,1000)n,(sig(i),i=1,4)
  write(6,1001)(sig(i),i=5,8)
  if (verbose.and.vlvl.ge.1) then
    write(*,1000)n,(sig(i),i=1,4)
    write(*,1001)(sig(i),i=5,8)
  end if
  mchk = 0; irep = 0
  call greens(ham,wt,tot,e,eps)
  call crete(dels,delp)
!verbose = .true.; vlvl = 2;
  if(verbose.and.vlvl.ge.1) print 1002, dels, delp
  if (abs(dble(dels(1))).le.cr.and.abs(aimag(dels(1))).le.ci.and. &
      abs(dble(dels(2))).le.cr.and.abs(aimag(dels(2))).le.ci.and. &
      abs(dble(delp(1))).le.cr.and.abs(aimag(delp(1))).le.ci.and. &
      abs(dble(delp(2))).le.cr.and.abs(aimag(delp(2))).le.ci.and. &
      abs(dble(delp(3))).le.cr.and.abs(aimag(delp(3))).le.ci.and. &
      abs(dble(delp(4))).le.cr.and.abs(aimag(delp(4))).le.ci.and. &
      abs(dble(delp(5))).le.cr.and.abs(aimag(delp(5))).le.ci.and. &
      abs(dble(delp(6))).le.cr.and.abs(aimag(delp(6))).le.ci) then
    mchk = 1
    if (verbose.and.vlvl.ge.2) print 1005
    do i = 1, nse
      if(aimag(sig(i)).gt.1.0d-20) then 
        irep = 1
        sig(i) = cmplx(dble(sig(i)),-aimag(sig(i)),8)
        if (verbose.and.vlvl.ge.2) print 1006, i,abs(aimag(sig(i)))
      else
        if (verbose.and.vlvl.ge.2.and.i.eq.1) print 1007
      end if
    end do
  else
     sig(1) = sig(1) + del1*dels(1)
     sig(5) = sig(5) + del1*dels(2)
     sig(2) = sig(2) + del1*delp(1)
     sig(3) = sig(3) + del1*delp(2)
     sig(4) = sig(4) + del1*delp(3)
     sig(6) = sig(6) + del1*delp(4)
     sig(7) = sig(7) + del1*delp(5)
     sig(8) = sig(8) + del1*delp(6)
  end if
!verbose = .false.; vlvl = 2;
  if (irep.eq.1) cycle
  if (mchk.eq.1) exit
end do
if (verbose) print 1004
return
1000 format(2X,I5,8F14.9)
1001 format(7X,8F14.9)
1002 format(4(2(F12.8,1X)))
1003 format(/,'Begin subroutine cpaNR')
1004 format('End cpaNR',/)
1005 format("Didn't you just converge?")
1006 format("Yes, but |imaginary[sig(",I2,")]| is greater than zero. ", &
  E15.8)
1007 format("Yes, you are correct")
end subroutine cpaNR
