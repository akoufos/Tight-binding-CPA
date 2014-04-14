!--------------------------------------------------------------------------
! Function for doing interpolation of the Fermi level
!--------------------------------------------------------------------------
real(kind=8) function interp(xx,x,f,j1,n)
implicit none
integer(kind=4) :: i, istart, j ,j2
integer(kind=4), intent(in) :: j1, n
real(kind=8) :: fx, p
real(kind=8), intent(in) :: f(1), x(1), xx
fx = 0.0d0
istart = j1 - n + n/2 + 1
j2 = istart + n - 1
do j = istart, j2
  p = f(j)
  do i = istart, j2
    if (i.eq.j) cycle
    p = p*(xx-x(i))/(x(j)-x(i))
  end do
  fx = fx + p
end do
interp = fx
return
end function interp

subroutine mainn()
!--------------------------------------------------------------------------
! Subroutine to run the main parts of the CPA program. Specifically this
! subroutine does the Newton-Raphson procedure, applies concentrations 
! and most of the initialization.
!--------------------------------------------------------------------------
! Variables:
! del - Temperature broadening ?
! epiv - Pivot energy (helps with N-R iterations [usually around E_F])
! eps - Imaginary part of energy shift ?
! emax/emin - Maximum and minimum of energy window, respectively
! sag**1 - Real and imaginary parts, respectively, of s & p initial onsite
  ! parameters (should be average between two substitution atoms)
! sig** - Real and imaginary parts, respectively, of s & p self-energies 
! sig - Complex self-energies
!--------------------------------------------------------------------------
use global
use hamiltonians
use sigma
implicit none
common /d2/ emin, emax, eps, del, epiv
integer(kind=4) :: i, l, ll
integer(kind=4) :: iss, iss1, itop, ixyz, l9, m, mchk, mc, mcm, mcount, &
  n, nchk, ndim, nmode, num99
integer(kind=4), parameter :: numit = 100
real(kind=8) :: del, dnorfl, e, efl, emax, emin, epiv, &
  eps, interp, nuelec, s, sagi1(nse), sagr1(nse), totvol, dos(sec+1)
real(kind=8) :: res(2000,9*ntype+2*nse+2), anumel(2000), dumm(2000), &
  edum(2000), dums1(2000), dump1(2000), dums2(2000), dump2(2000), & 
  dumxz(2000), dumxy(2000), dum3r(2000), dumx2(2000), densfl(10), &
  weight(jsz), qq(jsz,3)
complex(kind=8) :: sag(nse)
character(len=100) :: file1, file2
if (verbose) print 1000
open(6,file='cpaper.out',blank='zero')
nchk = 0
totvol = 0.0d0
write(file1,'(A)') 'cpamat1.dat'
write(file2,'(A)') 'cpamat2.dat'
verbose = .true.; vlvl = 3
call readin(ndim,nuelec,sagr1,sagi1)
verbose = .false.; vlvl = 2
call setOnsites
call kpts(jsz,qq,weight,totvol)
call readSec(hma,vsa,file1)
call readSec(hmb,vsb,file2)
call setInitHam
iss1 = 1
 1111 if(nchk.eq.2) goto 1234
sag(:) = cmplx(sagr1(:),sagi1(:),8)
if(nchk.eq.0) nmode = 1
if(nchk.eq.1) nmode = 2
if(nchk.eq.1) epiv = epiv - del
e = epiv
do ixyz = iss1, 2000
  sig(:) = sag(:)
  call cpaNR(weight,totvol,e,eps,mchk,numit)
  if (mchk.eq.1) then
    sag(:) = sig(:)
  else
    sag(:) = cmplx(sagr1(:),sagi1(:),8)
  end if
  if (mchk.eq.1) goto 9991
  write(6,5004)e
  write(*,5004)e
  goto 870 
  nchk = nchk + 1
  goto 1111
  9991 continue
  write(6,1015)n,(sig(i),i=1,4),e
  write(6,1016)  (sig(i),i=5,8),e
  verbose = .true.
  if (verbose) then
    write(*,1015)n,(sig(i),i=1,4),e
    write(*,1016)  (sig(i),i=5,8),e
  end if
  verbose = .false.
  call cpaDOS(dos,weight,totvol,e,eps)
870 continue
  res(ixyz,1) = e
  res(ixyz,2) = dble(sig(1))
  res(ixyz,3) = aimag(sig(1))
  res(ixyz,4) = dble(sig(2))
  res(ixyz,5) = aimag(sig(2))
  res(ixyz,6) = dble(sig(3))
  res(ixyz,7) = aimag(sig(3))
  res(ixyz,8) = dble(sig(4))
  res(ixyz,9) = aimag(sig(4))
  res(ixyz,10) = dble(sig(5))
  res(ixyz,11) = aimag(sig(5))
  res(ixyz,12) = dble(sig(6))
  res(ixyz,13) = aimag(sig(6))
  res(ixyz,14) = dble(sig(7))
  res(ixyz,15) = aimag(sig(7))
  res(ixyz,16) = dble(sig(8))
  res(ixyz,17) = aimag(sig(8))
  res(ixyz,18) = dos(1) + dos(10)
  res(ixyz,19) = sum(dos(2:4)) + sum(dos(11:13))
  res(ixyz,20) = sum(dos(5:6)) + sum(dos(14:15))
  res(ixyz,21) = dos(7) + dos(16)
  res(ixyz,22) = dos(8) + dos(17)
  res(ixyz,23) = dos(9) + dos(18)
  res(ixyz,24) = dos(19) + dos(28)
  res(ixyz,25) = sum(dos(20:22)) + sum(dos(29:31))
  res(ixyz,26) = dos(sec+1)
  write(6,1030)e,(sig(i),i=1,4),dos(sec+1)
  if(nmode.eq.2) goto 51
  e = e + del
  if(e.le.emax) cycle ! next iteration of ixyz loop
  iss = ixyz
  iss1 = iss + 1
  nchk = nchk + 1
  goto 1111
 51     e = e - del
  if(e.ge.emin) cycle ! next iteration of ixyz loop
  exit
end do ! End of ixyz loop
 1234 continue
itop = ixyz
write(6,3000)
write(6,3010)
l = itop + 1
l9 = 0
do ll = iss1, itop
  l = l - 1
  l9 = l9 + 1
  edum(l9) = res(l,1)
  dums1(l9) = res(l,18)*2.0d0
  dump1(l9) = res(l,19)*2.0d0
  dumxz(l9) = res(l,20)*2.0d0
  dumxy(l9) = res(l,21)*2.0d0
  dum3r(l9) = res(l,22)*2.0d0
  dumx2(l9) = res(l,23)*2.0d0
  dums2(l9) = res(l,24)*2.0d0
  dump2(l9) = res(l,25)*2.0d0
  dumm(l9) = res(l,26)*2.0d0
  write(6,3030)(res(l,m),m=1,9)
  write(6,3040)(res(l,m),m=10,17)
end do
do l = 1, iss
  l9 = l9 + 1
  edum(l9) = res(l,1)
  dums1(l9) = res(l,18)*2.0d0
  dump1(l9) = res(l,19)*2.0d0
  dumxz(l9) = res(l,20)*2.0d0
  dumxy(l9) = res(l,21)*2.0d0
  dum3r(l9) = res(l,22)*2.0d0
  dumx2(l9) = res(l,23)*2.0d0
  dums2(l9) = res(l,24)*2.0d0
  dump2(l9) = res(l,25)*2.0d0
  dumm(l9) = res(l,26)*2.0d0
  write(6,3030)(res(l,m),m=1,9)
  write(6,3040)(res(l,m),m=10,17)
end do
num99 = itop - 1
call simp(edum,dumm,anumel,num99,verbose,vlvl)
anumel(1) = 0.0d0
anumel(itop) = anumel(num99) + dumm(iss)*del
open(8,file='dosdat.cpa.plot')
write(6,3000)
write(6,3050)
write(8,3050)
l = itop + 1
l9 = 0
do ll = iss1, itop
  l = l - 1
  l9 = l9 + 1
  res(l,27) = anumel(l9)
  write(6,3070)res(l,1),(res(l,m),m=18,27),ll
  write(8,3070)res(l,1),(res(l,m),m=18,27),ll
end do
do l = 1, iss
  l9 = l9 + 1
  res(l,27) = anumel(l9)
  write(6,3070)res(l,1),(res(l,m),m=18,27),l
  write(8,3070)res(l,1),(res(l,m),m=18,27),l
end do
close(8)
mcount = 0
do l = 1, itop
  if(anumel(l).gt.nuelec) then
    if(mcount.ne.0) cycle
    mcount = l
  end if
end do
!     linear interpolation to find the fermi level
mc = mcount
mcm = mc - 1
s = nuelec
efl = interp(s,anumel,edum(1),mcm,2)
dnorfl = interp(s,anumel,dumm(1),mcm,2)/2.0d0
write(6,841)
densfl(1) = interp(s,anumel,dums1(1),mcm,2)/2.0d0
densfl(2) = interp(s,anumel,dump1(1),mcm,2)/2.0d0
densfl(3) = interp(s,anumel,dumxz(1),mcm,2)/2.0d0
densfl(4) = interp(s,anumel,dumxy(1),mcm,2)/2.0d0
densfl(5) = interp(s,anumel,dum3r(1),mcm,2)/2.0d0
densfl(6) = interp(s,anumel,dumx2(1),mcm,2)/2.0d0
densfl(7) = interp(s,anumel,dums2(1),mcm,2)/2.0d0
densfl(8) = interp(s,anumel,dump2(1),mcm,2)/2.0d0
write(6,839)
write(6,840)efl,nuelec,dnorfl,(densfl(i),i=1,8)
close(6)
if (verbose) print 2000
return
839  format(/,'Fermi energy   Electrons   Total DOS   Fe-s   Fe-p   Fe-&
  xy   Fe-xz   Fe-3r^2-z^2   Fe-x^2-y^2   Se-s   Se-p',//)
840  format(2f10.5   ,3x,7f10.5//)
841  format(//65x,11h (per spin)  )
1000 format(/,'Begin subroutine mainn')
2000 format('End mainn',/)
1015 format(5x,i5,8f15.8,f10.6)
1016 format(10X,8f15.8,f10.6)
1030 format(//1x,f9.5,8f12.8,10x, f10.5//)
3000 format(1h1,///)
3010 format(20x,'Energy Complex self-energies (s, px, py, pz)'///)
3030 format(F9.5,3X,4(2G12.5,2X))
3040 format(12X,4(2G12.5,2X))
3050 format(10x,'Energy, Fe-s, Fe-p(x,y,z), Fe-d(xz,yz), Fe-d(xz), Fe-&
  d(3r^2-z^2), Fe-d(x^2-y^2), Se-s, Se-p(x,y,z), Total DOS, electrons, &
  iteration',/)
3070 format(2x,12f10.4,1x,i5)
5004 format('   sigma didn t converge for e=',f10.4)
end subroutine mainn
