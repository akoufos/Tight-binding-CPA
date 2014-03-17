subroutine crete(delsr,delsi,delpr,delpi)
use global
implicit none
common /d1/ con
common /d3/ sen, grn
integer(4) :: l
real(8), intent(out) :: delsr(2), delsi(2), delpr(2), delpi(2)
real(8) :: con
complex(8) :: ans(nse,nse), ge(nse,nse), su(nse,nse), term(nse,nse), &
  usi(nse,nse), grn(sec,sec), sen(nse)
usi(:,:) = (0.0d0,0.0d0)
ge(:,:) = (0.0d0,0.0d0)
su(:,:) = (0.0d0,0.0d0)
delsr(:) = 0.0d0
delsi(:) = 0.0d0
delpr(:) = 0.0d0
delpi(:) = 0.0d0
do l = 1, nse
  su(l,l) = cmplx(1.0d0,0.0d0,8)
end do
ge(1,1) = grn(19,19)
ge(2,2) = grn(20,20)
ge(3,3) = grn(21,21)
ge(4,4) = grn(22,22)
ge(5,5) = grn(28,28)
ge(6,6) = grn(29,29)
ge(7,7) = grn(30,30)
ge(8,8) = grn(31,31)
usi(1,1) = sen(1)
usi(2,2) = sen(2)
usi(3,3) = sen(3)
usi(4,4) = sen(4)
usi(5,5) = sen(5)
usi(6,6) = sen(6)
usi(7,7) = sen(7)
usi(8,8) = sen(8)
call herakl(ge,usi,su,ans)
term(:,:) = ans(:,:)*con
delsr(1) = dble(term(1,1))
delsi(1) = aimag(term(1,1))
delpr(1) = dble(term(2,2))
delpi(1) = aimag(term(2,2))
delsr(2) = dble(term(5,5))
delsi(2) = aimag(term(5,5))
delpr(2) = dble(term(6,6))
delpi(2) = aimag(term(6,6))
return
end subroutine crete
