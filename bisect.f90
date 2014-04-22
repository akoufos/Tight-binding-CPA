subroutine bisect(G,sig,dels,delp)
use global
use onsites
implicit none
integer(kind=4) :: l
complex(kind=8) :: ge(nse), term(nse), F(nse), dF(nse), alpha(nse), &
  beta(nse)
complex(kind=8), intent(in) :: G(sec,sec)
complex(kind=8), intent(inout) :: sig(2,nse)
complex(kind=8), intent(out) :: dels(2), delp(6)
if (verbose) print 1000
dels(:) = (0.0d0,0.0d0)
delp(:) = (0.0d0,0.0d0)
if(verbose.and.vlvl.ge.1) print 1001, sig
ge(1) = G(19,19)
ge(2) = G(20,20)
ge(3) = G(21,21)
ge(4) = G(22,22)
ge(5) = G(28,28)
ge(6) = G(29,29)
ge(7) = G(30,30)
ge(8) = G(31,31)
if (verbose) print 2000
1000 format(/,'Begin subroutine bisect')
1001 format('Sigma: ',/,4(2(F10.6,1X)))
2000 format('End bisect',/)
end subroutine bisect
