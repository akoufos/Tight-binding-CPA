!--------------------------------------------------------------------------
! Function for doing interpolation of the Fermi level
!--------------------------------------------------------------------------
real(8) function interp(xx,x,f,j1,n)
implicit none
integer(4) :: i, istart, j ,j2
integer(4), intent(in) :: j1, n
real(8) :: fx, p
real(8), intent(in) :: f(1), x(1), xx
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
! grn - Green's function matrix
! sag**1 - Real and imaginary parts, respectively, of s & p initial onsite
  ! parameters (should be average between two substitution atoms)
! sig** - Real and imaginary parts, respectively, of s & p self-energies 
! sig - Complex self-energies
!--------------------------------------------------------------------------
use omp_lib
use global
use hamiltonians
implicit none
common /d2/ emin, emax, eps, del, epiv, sagsr1, sagsi1, sagpr1, sagpi1
common /d3/ sig, grn
integer(4) :: i, l, ll
integer(4) :: iss, iss1, itop, ixyz, l9, m, mchk, mc, mcm, mcount, &
  n, nchk, ndim, nmode, num99
integer(4), parameter :: numit = 100
real(8) :: del, dnorfl, e, efl, emax, emin, epiv, &
  eps, interp, nuelec, s, sagpi1, sagpr1, sagsi1, sagsr1, totvol, &
  dos(sec+1)
real(8) :: res(2000,15), anumel(2000), dumm(2000), edum(2000), &
  dums1(2000), dump1(2000), dums2(2000), dump2(2000), dumxz(2000), &
  dumxy(2000), dum3r(2000), dumx2(2000), densfl(10), weight(jsz), &
  qq(jsz,3)
real(8), parameter :: dsig = 1.0d-4
complex(8) :: ham(jsz,sec,sec), grn(sec,sec), &
  sig(nse), sags(2), sagp(6), del1
character(len=100) :: file1, file2
if (verbose) print 1000
open(6,file='cpaper.out',blank='zero')
nchk = 0
totvol = 0.0d0
write(file1,'(A)') 'cpamat1.dat'
write(file2,'(A)') 'cpamat2.dat'
call readin(ndim, nuelec)
call kpts(jsz,qq,weight,totvol)
call readSec(hma,vsa,file1)
call readSec(hmb,vsb,file2)
iss1 = 1
 1111 if(nchk.eq.2) goto 1234
sags(:) = cmplx(sagsr1,sagsi1,8)
sagp(:) = cmplx(sagpr1,sagpi1,8)
if(nchk.eq.0) nmode = 1
if(nchk.eq.1) nmode = 2
if(nchk.eq.1) epiv = epiv - del
e = epiv
del1 = cmplx(dsig,dsig,8)
do ixyz = iss1, 2000
  do i = 1, 4
    if(i.eq.1) then
      sig(1) = sags(1)
      sig(5) = sags(2)
    else
      sig(i) = sagp(i-1)
      sig(i+4) = sagp(i+2)
    end if
  end do
  call cpaNR(ham,weight,totvol,e,eps,del1,mchk,numit)
  do i = 1, 4
    if(i.eq.1) then
      sags(1) = sig(1)
      sags(2) = sig(5)
    else
      sagp(i-1) = sig(i)
      sagp(i+2) = sig(i+4)
    end if
  end do
  if (mchk.eq.1) goto 9991
  write(6,5004)e
  goto 870 
  nchk = nchk + 1
  goto 1111
  9991 continue
  write(6,1015)n,(sig(i),i=1,4),e
  write(6,1016)  (sig(i),i=5,8),e
  call cpaDOS(dos,weight,totvol,e,eps)
verbose = .false.; vlvl = 1
870 continue
  res(ixyz,1) = e
  res(ixyz,2) = dble(sig(1))
  res(ixyz,3) = aimag(sig(1))
  res(ixyz,4) = dble(sig(2))
  res(ixyz,5) = aimag(sig(2))
  res(ixyz,6) = dos(1) + dos(10)
  res(ixyz,7) = sum(dos(2:4)) + sum(dos(11:13))
  res(ixyz,8) = sum(dos(5:6)) + sum(dos(14:15))
  res(ixyz,9) = dos(7) + dos(16)
  res(ixyz,10) = dos(8) + dos(17)
  res(ixyz,11) = dos(9) + dos(18)
  res(ixyz,12) = dos(19) + dos(28)
  res(ixyz,13) = sum(dos(20:22)) + sum(dos(29:31))
  res(ixyz,14) = dos(sec+1)
  write(6,1030)e,sig(1),sig(2),dos(sec+1)
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
  dums1(l9) = res(l,6)*2.0d0
  dump1(l9) = res(l,7)*2.0d0
  dumxz(l9) = res(l,8)*2.0d0
  dumxy(l9) = res(l,9)*2.0d0
  dum3r(l9) = res(l,10)*2.0d0
  dumx2(l9) = res(l,11)*2.0d0
  dums2(l9) = res(l,12)*2.0d0
  dump2(l9) = res(l,13)*2.0d0
  dumm(l9) = res(l,14)*2.0d0
  write(6,3030)(res(l,m),m=1,5)
end do
do l = 1, iss
  l9 = l9 + 1
  edum(l9) = res(l,1)
  dums1(l9) = res(l,6)*2.0d0
  dump1(l9) = res(l,7)*2.0d0
  dumxz(l9) = res(l,8)*2.0d0
  dumxy(l9) = res(l,9)*2.0d0
  dum3r(l9) = res(l,10)*2.0d0
  dumx2(l9) = res(l,11)*2.0d0
  dums2(l9) = res(l,12)*2.0d0
  dump2(l9) = res(l,13)*2.0d0
  dumm(l9) = res(l,14)*2.0d0
  write(6,3030)(res(l,m),m=1,5)
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
  res(l,15) = anumel(l9)
  write(6,3070)res(l,1),(res(l,m),m=6,15),ll
  write(8,3070)res(l,1),(res(l,m),m=6,15),ll
end do
do l = 1, iss
  l9 = l9 + 1
  res(l,15) = anumel(l9)
  write(6,3070)res(l,1),(res(l,m),m=6,15),l
  write(8,3070)res(l,1),(res(l,m),m=6,15),l
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
839  format(/,'Fermi energy   Electrons   Total DOS   Fe-s   Fe-p   Fe- &
  xy   Fe-xz   Fe-3r^2-z^2   Fe-x^2-y^2   Se-s   Se-p',//)
840  format(2f10.5   ,3x,7f10.5//)
841  format(//65x,11h (per spin)  )
1000 format(/,'Begin subroutine mainn')
2000 format('End mainn',/)
1015 format(5x,i5,8f15.8,f10.6)
1016 format(10X,8f15.8,f10.6)
1030 format(//1x,f9.5,4f12.8,10x, f10.5//)
3000 format(1h1,///)
3010 format(20x,' printing of --- energy --complex sigs  ---- complex &
  sigp'///)
3030 format(24x,f10.4,4f15.8)
3050 format(10x,'Energy, Fe-s, Fe-p(x,y,z), Fe-d(xz,yz), Fe-d(xz), Fe-&
  d(3r^2-z^2), Fe-d(x^2-y^2), Se-s, Se-p(x,y,z), Total DOS, electrons, &
  iteration',/)
3070 format(2x,11f10.4,1x,i5)
5004 format('   sigma didn t converge for e=',f10.4)
end subroutine mainn
