subroutine herakl(ge,u,su,term)
use global
implicit none
complex(8) :: ans(nse,nse), up(nse,nse)
complex(8), intent(in) :: ge(nse,nse), su(nse,nse), u(nse,nse)
complex(8), intent(out) :: term(nse,nse)
up = matmul(ge,u)
ans(:,:) = su(:,:) - up(:,:)
call cmplxInv(ans,nse,verbose)
term = matmul(u,ans)
return
end subroutine herakl
