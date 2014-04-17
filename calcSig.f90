subroutine calcSig(wt,tot,e,eps,mchk,numit,method)
!--------------------------------------------------------------------------
! Solves for the self-energies, and thus the Green's function, using the
! Newton-Raphson method.
!--------------------------------------------------------------------------
! Variables:
! G(sec,sec) - Green's function
! sig(nse) - Self-energies of the disordered states
! H(jsz,sec,sec) - Hamiltonian of the system
! wt(jsz) - Weights at each k-point in the Brillouin zone
! tot - Total weight from all k-points
! dels(2) - Change in the self-energies of the disordered s states
! delp(6) - Change in the self-energies of the disordered p states
! mchk - Flag to specify whether routine converged or not
! numit - Maximum number of iterations for the procedure
! method - String used to decide which zero finding procedure to use
!--------------------------------------------------------------------------
use global
use converge
use hamiltonians, only : ham
use sigma
implicit none
integer(kind=4) :: i, n, irep
integer(kind=4), intent(in) :: numit
integer(kind=4), intent(out) :: mchk
real(kind=8), intent(in) :: wt(jsz), tot, e, eps
complex(kind=8) :: dels(2), delp(6), H(jsz,sec,sec), G(sec,sec), &
  sigs(2,nse)
character(len=100), intent(in) :: method
do n = 1, numit
  if (verbose) print 1000
  dels(:) = (0.0d0,0.0d0); delp(:) = (0.0d0,0.0d0)
  write(6,1001)n,(sig(i),i=1,4)
  write(6,1002)(sig(i),i=5,8)
  if (verbose.and.vlvl.ge.1) then
    write(*,1001)n,(sig(i),i=1,4)
    write(*,1002)(sig(i),i=5,8)
  end if
  mchk = 0; irep = 0
  sigs(:,:) = (0.0d0,0.0d0)
  H(:,:,:) = -ham(:,:,:)
  call setHam(H,e,eps)
  call greens(H,G,wt,tot)
!  if (n.eq.1) then
!    write(9,1004) e, G(19,19), G(20,20), G(21,21), G(22,22), G(28,28), &
!      G(29,29), G(30,30), G(31,31), (sig(i), i=1,4)
!  end if
  select case (method)
    case ('Newton')
      sigs(1,:) = sig(:)
      call newton(G,sigs(1,:),dels,delp)
    case ('Fixed')
      sigs(1,:) = sig(:)
      call fixpt(G,dels,delp)
    case ('False Position')
!      call falsi(G,dels,delp)
    case ('Bisect')
!      call bisect(G,dels,delp)
    case default
      sigs(1,:) = sig(:)
      call newton(G,sigs(1,:),dels,delp)
  end select
  sig(:) = sigs(1,:)
  if(verbose.and.vlvl.ge.1) print 1003, dels, delp
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
      if (aimag(sig(i)).gt.1.0d-15) then 
        irep = 1
        if (verbose.and.vlvl.ge.2) print 1006, i,abs(aimag(sig(i)))
        sig(i) = cmplx(dble(sig(i)),-aimag(sig(i)),8)
      elseif (aimag(sig(i)).gt.0.0d0.and.aimag(sig(i)).lt.1.0d0-15) then
        if (verbose.and.vlvl.ge.2) print 1006, i,abs(aimag(sig(i)))
        sig(i) = cmplx(dble(sig(i)),0.0d0,8)
      else
        if (verbose.and.vlvl.ge.2.and.i.eq.1) print 1007
      end if
    end do
  end if
  if (irep.eq.1) cycle
  if (mchk.eq.1.or.n.eq.numit) then
    write(9,1004) e, G(19,19), G(20,20), G(21,21), G(22,22), G(28,28), &
      G(29,29), G(30,30), G(31,31), (sig(i), i=1,4)
    if (verbose) print 2000
    exit
  end if
end do
if (verbose) print 2000
return
1000 format(/,'Begin subroutine calcSig')
1001 format(2X,I5,4(F14.9,E14.6))
1002 format(7X,4(F14.9,E14.6))
1003 format(4(2(F12.8,1X)))
1004 format(F8.5,1X,12(2(F12.8,1X)))
1005 format("Didn't you just converge?")
1006 format("Yes, but |imaginary[sig(",I2,")]| is greater than zero. ", &
  E15.8)
1007 format("Yes, you are correct")
2000 format('End calcSig',/)
end subroutine calcSig
