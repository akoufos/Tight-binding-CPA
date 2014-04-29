subroutine calcSig(wt,tot,e,eps,mchk,numit,method,sagr1,sagi1)
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
integer(kind=4) :: i, n, irep, reset(nse)
integer(kind=4), intent(in) :: numit
logical, intent(out) :: mchk(nse)
real(kind=8), intent(in) :: wt(jsz), tot, e, eps, sagr1(nse), sagi1(nse)
complex(kind=8) :: del(nse), H(jsz,sec,sec), G(sec,sec), &
  sigs(2,nse)
character(len=100), intent(in) :: method
reset = 0
do n = 1, numit
  if (verbose) print 1000
  del(:) = (0.0d0,0.0d0)
  write(7,1001)n,(sig(i),i=1,4)
  write(7,1002)(sig(i),i=5,8)
  if (verbose.and.vlvl.ge.1) then
    write(*,1001)n,(sig(i),i=1,4)
    write(*,1002)(sig(i),i=5,8)
  end if
  irep = 0
  mchk(:) = .false.
  sigs(:,:) = (0.0d0,0.0d0)
  H(:,:,:) = -ham(:,:,:)
  call setHam(H,e,eps)
  call greens(H,G,wt,tot)
!  if (n.eq.1) then
!    write(9,1004) e, G(19,19), G(20,20), G(21,21), G(22,22), G(28,28), &
!      G(29,29), G(30,30), G(31,31), (sig(i), i=1,4)
!  end if
  if (verbose.and.vlvl.ge.1) print 1008, trim(method)
  select case (method)
    case ('Newton')
      sigs(1,:) = sig(:)
      call newton(G,sigs(1,:),del)
    case ('Fixed')
      sigs(1,:) = sig(:)
      call fixpt(G,sigs,del)
    case ('False Position')
!      call falsi(G,del)
    case ('Bisect')
!      call bisect(G,del)
    case default
      sigs(1,:) = sig(:)
      call newton(G,sigs(1,:),del)
  end select
  sig(:) = sigs(1,:)
  if(verbose.and.vlvl.ge.1) print 1003, del
  do i = 1, nse
    if (abs(dble(del(i))).le.cr.and.abs(aimag(del(i))).le.ci) &
      mchk(i) = .true.
  end do
  if (all(mchk)) then
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
  else
    do i = 1, nse
      if (dble(sig(i)).lt.-4.0d0*abs(sagr1(i)).or. &
        dble(sig(i)).gt.4.0d0*abs(sagr1(i))) then
        reset(i) = reset(i) + 1
        print *, 'RESET: ',reset
        if (reset(i).ge.2) return
        sig(i) = cmplx(sagr1(i),sagi1(i),8)
      end if
    end do
  end if
  if (irep.eq.1) cycle
  if (all(mchk).or.n.eq.numit) then
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
1008 format(/,"Running ",A," method for finding the self-energies.",/)
2000 format('End calcSig',/)
end subroutine calcSig
