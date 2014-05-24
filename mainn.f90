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
! method - String used to decide which zero finding procedure to use
!--------------------------------------------------------------------------
use global
use hamiltonians
use sigma
implicit none
common /d2/ emin, emax, eps, del, epiv
integer(kind=4) :: i, l, ll, mode
integer(kind=4) :: iss, iss1, itop, ixyz, l9, m, mc, mcm, mcount, &
  n, nchk, ndim, nmode, num99
integer(kind=4), parameter :: numit = 100
logical :: mchk(nse)
real(kind=8) :: del, dnorfl, e, efl, emax, emin, epiv, eps, interp, &
  nuelec, s, sagi1(nse), sagr1(nse), totvol, dos(sec+1), &
  dos2(ntype,sec)
real(kind=8) :: res(8000,50), anumel(8000), dumm(8000), &
  edum(8000), dums1(8000), dump1(8000), dums2(8000), dump2(8000), & 
  dumxz(8000), dumxy(8000), dum3r(8000), dumx2(8000), dumsA(8000), &
  dumpA(8000), dumsB(8000), dumpB(8000), densfl(12), weight(jsz), &
  qq(jsz,3), spec(8000,jsz)
complex(kind=8) :: sag(nse), sagcon(nse)
character(len=100) :: file1, file2, method
if (verbose) print 1000
nchk = 0; iss1 = 1; n = 0
totvol = 0.0d0
open(7,file='cpaper.out',blank='zero')
write(file1,'(A)') 'cpamat1.dat'
write(file2,'(A)') 'cpamat2.dat'
verbose = .true.; vlvl = 3
call readin(ndim,nuelec,sagr1,sagi1,mode)
! If mode 2, program already run; Calculate superconductivity only
if (mode.eq.2) goto 9999
call setOnsites
verbose = .false.; vlvl = 2
call kpts(jsz,qq,weight,totvol)
call readSec(hma,vsa,file1)
call readSec(hmb,vsb,file2)
call setInitHam
open(9,file='green.dat',blank='zero')
 1111 if(nchk.eq.2) goto 1234
sag(:) = cmplx(sagr1(:),sagi1(:),8)
if(nchk.eq.0) nmode = 1
if(nchk.eq.1) nmode = 2
if(nchk.eq.1) epiv = epiv - del
e = epiv
sig(:) = sag(:)
do ixyz = iss1, 8000
!  sag(:) = cmplx(sagr1(:),sagi1(:),8)
  sig(:) = sag(:)
  write(method,'(A)') 'Newton'
  call calcSig(weight,totvol,e,eps,mchk,numit,method,sagr1,sagi1)
  if (all(mchk)) then
    sag(:) = sig(:)
    sagcon(:) = sagcon(:) + sig(:)
    n = n + 1
!    goto 9991
  else ! Didn't converge so estimate with some other self-energies
    do l = 1, nse
      if (mchk(l).eqv..false.) then
      ! Estimate with average of good self-energies
!        sig(l) = sagcon(l)/dble(n)
      ! Estimate next self-energy as the average of the good one
        sag(l) = sagcon(l)/dble(n)
        print *, 'Using average sig(l): ',l, sagcon(l)/dble(n)
      end if
!      if (mchk(l).eq..false.) sig(l) = cmplx(sagr1(l),sagi1(l),8)
    end do
    print 1002
    write(7,5004)e
    write(*,5004)e
    goto 870 
  end if
!  nchk = nchk + 1
!  goto 1111
!  9991 continue
  write(7,1015)n,(sig(i),i=1,4),e
  write(7,1016)  (sig(i),i=5,8),e
  verbose = .true.
  if (verbose) then
    write(*,1015)n,(sig(i),i=1,4),e
    write(*,1016)  (sig(i),i=5,8),e
  end if
  verbose = .false.
  call cpaDOS(dos,dos2,spec(n,:),weight,totvol,e,eps)
  res(n,1) = e
  res(n,2) = dble(sig(1))
  res(n,3) = aimag(sig(1))
  res(n,4) = dble(sig(2))
  res(n,5) = aimag(sig(2))
  res(n,6) = dble(sig(3))
  res(n,7) = aimag(sig(3))
  res(n,8) = dble(sig(4))
  res(n,9) = aimag(sig(4))
  res(n,10) = dble(sig(5))
  res(n,11) = aimag(sig(5))
  res(n,12) = dble(sig(6))
  res(n,13) = aimag(sig(6))
  res(n,14) = dble(sig(7))
  res(n,15) = aimag(sig(7))
  res(n,16) = dble(sig(8))
  res(n,17) = aimag(sig(8))
  res(n,18) = dos(1) + dos(10) !Fe-s
  res(n,19) = sum(dos(2:3)) + sum(dos(11:12)) !Fe-p(x,y)
  res(n,20) = dos(4) + dos(13) !Fe-p(z)
  res(n,21) = sum(dos(5:6)) + sum(dos(14:15)) !Fe-d(xz,yz)
  res(n,22) = dos(7) + dos(16) !Fe-d(xy)
  res(n,23) = dos(8) + dos(17) !Fe-d(3r^2-z^2)
  res(n,24) = dos(9) + dos(18) !Fe-d(x^2+y^2)
  res(n,25) = dos(19) + dos(28) !Se/Te-s
  res(n,26) = sum(dos(20:21)) + sum(dos(29:30)) !Se/Te-p(x,y)
  res(n,27) = dos(22) + dos(31) !Se/Te-p(z)
  res(n,28) = dos2(1,19) + dos2(1,28) !Se-s
  res(n,29) = sum(dos2(1,20:21)) + sum(dos2(1,29:30)) !Se-p(x,y)
  res(n,30) = dos2(1,22) + dos2(1,31) !Se-p(z)
  res(n,31) = dos2(2,19) + dos2(2,28) !Te-s
  res(n,32) = sum(dos2(2,20:21)) + sum(dos2(2,29:30)) !Te-p(x,y)
  res(n,33) = dos2(2,22) + dos2(2,31) !Te-p(z)
  res(n,34) = dos(sec+1) !Total DOS
  write(7,1030)e,(sig(i),i=1,4),dos(sec+1)
870 continue
  if(nmode.eq.2) goto 51
  e = e + del
  if(e.le.emax) cycle ! next iteration of n loop
  iss = n
  iss1 = iss + 1
  nchk = nchk + 1
  goto 1111
 51     e = e - del
  if(e.ge.emin) cycle ! next iteration of ixyz loop
  exit
end do ! End of ixyz loop
 1234 continue
itop = n
open(8,file='sigma.dat')
open(9,file='spectral.dat')
write(7,3000)
write(7,3010)
l = itop + 1
l9 = 0
do ll = iss1, itop
  l = l - 1
  l9 = l9 + 1
  edum(l9) = res(l,1)
  dums1(l9) = res(l,18)*2.0d0
  dump1(l9) = (res(l,19) + res(l,20))*2.0d0
  dumxz(l9) = res(l,21)*2.0d0
  dumxy(l9) = res(l,22)*2.0d0
  dum3r(l9) = res(l,23)*2.0d0
  dumx2(l9) = res(l,24)*2.0d0
  dums2(l9) = res(l,25)*2.0d0
  dump2(l9) = (res(l,26) + res(l,27))*2.0d0
  dumsA(l9) = res(l,28)*2.0d0
  dumpA(l9) = (res(l,29) + res(l,30))*2.0d0
  dumsB(l9) = res(l,31)*2.0d0
  dumpB(l9) = (res(l,32) + res(l,33))*2.0d0
  dumm(l9) = res(l,34)*2.0d0
  write(7,3030)(res(l,m),m=1,9)
  write(7,3040)(res(l,m),m=10,17)
  write(8,3030)(res(l,m),m=1,9)
  write(8,3040)(res(l,m),m=10,17)
  write(9,3045)res(l,1),(spec(l,m),m=1,jsz)
end do
do l = 1, iss
  l9 = l9 + 1
  edum(l9) = res(l,1)
  dums1(l9) = res(l,18)*2.0d0
  dump1(l9) = (res(l,19) + res(l,20))*2.0d0
  dumxz(l9) = res(l,21)*2.0d0
  dumxy(l9) = res(l,22)*2.0d0
  dum3r(l9) = res(l,23)*2.0d0
  dumx2(l9) = res(l,24)*2.0d0
  dums2(l9) = res(l,25)*2.0d0
  dump2(l9) = (res(l,26) + res(l,27))*2.0d0
  dumsA(l9) = res(l,28)*2.0d0
  dumpA(l9) = (res(l,29) + res(l,30))*2.0d0
  dumsB(l9) = res(l,31)*2.0d0
  dumpB(l9) = (res(l,32) + res(l,33))*2.0d0
  dumm(l9) = res(l,34)*2.0d0
  write(7,3030)(res(l,m),m=1,9)
  write(7,3040)(res(l,m),m=10,17)
  write(8,3030)(res(l,m),m=1,9)
  write(8,3040)(res(l,m),m=10,17)
  write(9,3045)res(l,1),(spec(l,m),m=1,jsz)
end do
close(8)
close(9)
num99 = itop - 1
call simp(edum,dumm,anumel,num99,verbose,vlvl)
anumel(1) = 0.0d0
anumel(itop) = anumel(num99) + dumm(iss)*del
open(8,file='dosdat.cpa.plot')
write(7,3000)
write(7,3050)
write(8,3060)
l = itop + 1
l9 = 0
do ll = iss1, itop
  l = l - 1
  l9 = l9 + 1
  res(l,35) = anumel(l9)
  write(7,3070)res(l,1),res(l,18),res(l,19)+res(l,20),res(l,21),&
    res(l,22),res(l,23),res(l,24),res(l,25),res(l,26)+res(l,27),&
    res(l,34),res(l,35),ll
  write(8,3080)res(l,1),(res(l,m),m=18,24),(res(l,m),m=28,35),ll
end do
do l = 1, iss
  l9 = l9 + 1
  res(l,35) = anumel(l9)
  write(7,3070)res(l,1),res(l,18),res(l,19)+res(l,20),res(l,21),&
    res(l,22),res(l,23),res(l,24),res(l,25),res(l,26)+res(l,27),&
    res(l,34),res(l,35),l
  write(8,3080)res(l,1),(res(l,m),m=18,24),(res(l,m),m=28,35),l
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
open(10,file='dosapw.itp')
write(10,1003)title
mc = mcount
mcm = mc - 1
s = nuelec
efl = interp(s,anumel,edum(1),mcm,2)
dnorfl = interp(s,anumel,dumm(1),mcm,2)/2.0d0
write(7,841)
densfl(1) = interp(s,anumel,dums1(1),mcm,2)/2.0d0
densfl(2) = interp(s,anumel,dump1(1),mcm,2)/2.0d0
densfl(3) = interp(s,anumel,dumxz(1),mcm,2)/2.0d0
densfl(4) = interp(s,anumel,dumxy(1),mcm,2)/2.0d0
densfl(5) = interp(s,anumel,dum3r(1),mcm,2)/2.0d0
densfl(6) = interp(s,anumel,dumx2(1),mcm,2)/2.0d0
densfl(7) = interp(s,anumel,dums2(1),mcm,2)/2.0d0
densfl(8) = interp(s,anumel,dump2(1),mcm,2)/2.0d0
densfl(9) = interp(s,anumel,dumsA(1),mcm,2)/2.0d0
densfl(10) = interp(s,anumel,dumpA(1),mcm,2)/2.0d0
densfl(11) = interp(s,anumel,dumsB(1),mcm,2)/2.0d0
densfl(12) = interp(s,anumel,dumpB(1),mcm,2)/2.0d0
write(7,839)
write(7,840)efl,nuelec,dnorfl,(densfl(i),i=1,8)
write(10,842)efl,nuelec,dnorfl,(densfl(i),i=1,12)
close(7)
close(9)
close(10)
9999 continue
call cpaGG(mode,efl,nuelec,dnorfl,densfl)
if (verbose) print 2000
return
839  format(/,'Fermi energy   Electrons   Total DOS   Fe-s   Fe-p   Fe-&
  xy   Fe-xz   Fe-3r^2-z^2   Fe-x^2-y^2   Se/Te-s   Se/Te-p',//)
840  format(2f10.5   ,3x,7f10.5//)
841  format(//65x,11h (per spin)  )
842  format(15(F10.6,1x))
1000 format(/,'Begin subroutine mainn')
1001 format(/,"Didn't converge for energy ",F8.5,/,"Trying ",A, &
  " method instead",/)
1002 format(/,"Still unable to converge. Something is wrong",/)
1003 format(A75)
1015 format(5X,I5,8F15.8,F10.6)
1016 format(10X,8F15.8,F10.6)
1030 format(//1X,F9.5,8F12.8,10X,F10.5//)
2000 format('End mainn',/)
3000 format(1H1,///)
3010 format(20X,'Energy Complex self-energies (s, px, py, pz)'///)
3030 format(F9.5,3X,4(2(G12.5,1X),2X))
3040 format(12X,4(2(G12.5,1X),2X))
3045 format(1F10.6,10X,1000(F10.6,1X))
3050 format(10x,'Energy, Fe-s, Fe-p(x,y,z), Fe-d(xz,yz), Fe-d(xz), Fe-&
  d(3r^2-z^2), Fe-d(x^2-y^2), Se/Te-s, Se/Te-p(x,y,z), &
  Total DOS, electrons, iteration',/)
3060 format(10x,'Energy, Fe-s, Fe-p(x,y), Fe-p(z), Fe-d(xz,yz), &
  Fe-d(xz), Fe-d(3r^2-z^2), Fe-d(x^2-y^2), Se-s, Se-p(x,y), Se-p(z), &
  Te-s, Te-p(x,y), Te-p(z), Total DOS, electrons, iteration',/)
3070 format(2x,11f10.4,1x,i5)
3080 format(2x,16f10.4,1x,i5)
5004 format('   sigma didn t converge for e=',f10.4,/)
end subroutine mainn
