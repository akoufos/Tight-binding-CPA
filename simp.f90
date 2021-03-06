subroutine simp(xx, fx, ax, nx, verbose, vlvl)
!--------------------------------------------------------------------------
! Need to verify the purpose of this subroutine.
!--------------------------------------------------------------------------
! Variables:
! verbose - Logical for debugging flags (.true. = debug info on)
! vlvl - Level of debugging verboseness
!--------------------------------------------------------------------------
implicit none
integer(kind=4) :: ix
integer(kind=4), intent(in) :: nx, vlvl
logical :: verbose
real(kind=8) :: delx
real(kind=8), intent(in) :: fx(2000), xx(2000)
real(kind=8), intent(out) :: ax(2000)
if (verbose) print 1000
delx = xx(2) - xx(1)
ax(1) = 0.0d0
do ix = 2, nx, 2
  select case(nx-ix)
    case(0)
      ax(nx) = delx*(-fx(nx-2)+8.0d0*fx(nx-1)+5.0d0*fx(nx))/12.0d0 &
      + ax(nx-1)
      return
    case default
      ax(ix) = delx*(5.0d0*fx(ix-1)+8.0d0*fx(ix)-fx(ix+1))/12.0d0 &
      + ax(ix-1)
      ax(ix+1) = delx*(fx(ix-1)+4.0d0*fx(ix)+fx(ix+1))/3.0d0 + ax(ix-1)
   end select
end do
if (verbose.and.vlvl.ge.2) print 1002, ax
if (verbose) print 2000
return
1000 format(/,'Begin subroutine simp')
1002 format('Electrons',F10.7)
2000 format('End simp',/)
end subroutine simp
