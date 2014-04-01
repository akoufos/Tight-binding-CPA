subroutine herakl(ge,u,su,term)
!--------------------------------------------------------------------------
! Subroutine to calculate the inverted Green's matrix.
!--------------------------------------------------------------------------
! up(nse,nse) - Product matrix of Green's and self-energy matrices
! ans(nse,nse) - Matrix to be inverted.
!--------------------------------------------------------------------------
use global
implicit none
complex(kind=8) :: ans(nse,nse), up(nse,nse)
complex(kind=8), intent(in) :: ge(nse,nse), su(nse,nse), u(nse,nse)
complex(kind=8), intent(out) :: term(nse,nse)
up = matmul(ge,u)
ans(:,:) = su(:,:) - up(:,:)
call cmplxInv(ans,nse,verbose)
term = matmul(u,ans)
return
end subroutine herakl
