subroutine cpaGG(mode, fermi, nelec, nef, N)
!--------------------------------------------------------------------------
! This subroutine will read the CPA DOS at the Fermi level to approximate
! the eta, Hopfield, parameter for use with the Gaspari-Gyorffy theory.
! Using the eta parameter, the subroutine will then calculate the
! electron-phonon coupling constant to use in the McMillian equation.
!--------------------------------------------------------------------------
! Variables:
! effl - Total electron-phonon coupling constant
! effw2 - Effective average squared phonon frequency
! eta(ntype) - Hopfield parameter for each atom type
! fermi - Fermi level
! lambda(ntype) - Electron-phonon coupling constant for each atom type
! M - Atomic mass used in electron-phonon calculation
! mode - Subroutine run mode to calculate superconductivity properties
!   1 - Calculate during CPA program
!   2 - Calculate at a later date of CPA completion
! N(12) - Angular momentum components of the DOS
! nef - Total DOS at the Fermi level
! nelec - Number of electrons
! Tc - Critical superconducting temperature of the system
! w2(ntype) - Squared average phonon frequency of the atom types
!--------------------------------------------------------------------------
! Notes:
!
! 05-23-2014 - Subroutine only considers the major contributions to the 
! DOS at the Fermi level when calculating the Hopfield parameter. For 
! the FeSeTe system, those components are the Fe-d DOS and the Se/Te-p
! DOS.
!
!--------------------------------------------------------------------------
use global
implicit none
integer(kind=4) :: i, j
integer(kind=4), intent(in) :: mode
real(kind=8) :: eta(ntype), lambda(ntype), Tc, M, w2(ntype), effw2, effl, &
  ratio, sfratio, theta, w, lsf, must
real(kind=8), intent(inout) :: fermi, nelec, nef, N(12)
if (verbose) print 1000
open(11,file='CPA-GG.out')
select case (mode)
  case (1)
    print *, 'Continuing superconductivity calculations with results'
  case (2)
    open(10,file='dosapw.itp')
    read(10,*)
    read(10,1001)fermi,nelec,NEF,(N(i),i=1,12)
    close(10)
  case default
    print *, 'Problem with mode. Nothing was done'
    return ! Return to mainn subroutine
end select
write(11,1100)title
write(11,1101)
write(11,1102)fermi,nelec,nef,(N(i),i=1,12)
eta(1) = nef*(sum(N(3:6))/nef)
eta(2) = nef*(N(8)/nef)
write(11,1103)'Hopfield parameter of Fe:',eta(1)
write(11,1103)'Hopfield parameter of Se/Te:',eta(2)
write(11,1103)'Total eta:',sum(eta)
open(10,file='superCPA.in')
read(10,*)M
read(10,*)w2(1), w2(2)
read(10,*)must
read(10,*)lsf
close(10)
write(11,*)
write(11,1103)'Effective Atomic Mass:',M
write(11,*)
write(11,1103)'Average squared phonon frequency of FeSe:',w2(1)
write(11,1103)'Average squared phonon frequency of FeTe:',w2(2)
! Calculate the effective average sqaured phonon frequency
effw2 = sum(w2)
w = sqrt(effw2)
write(11,1103)'Effective squared average phonon frequency:',effw2
write(11,1103)'Effective phonon frequency:',w
write(11,*)
! Calculate the electron-phonon coupling constant
lambda(1) = eta(1)/(M*w2(1))
lambda(2) = eta(2)/(M*w2(2))
effl = sum(lambda)
write(11,1103)'Electron-phonon coupling constant of Fe:',lambda(1)
write(11,1103)'Electron-phonon coupling constant of Se/Te:',lambda(2)
write(11,1103)'Total electron-phonon coupling constant:',effl
! Calculate the critical superconducting temperature
ratio = -1.04d0*(1.0d0+effl) / (effl-must*(1.0d0+0.62d0*effl))
sfratio = -1.04d0*(1.0d0+effl+lsf) / &
  (effl-must*(1.0d0+0.62d0*(effl+lsf))-lsf)
Tc = w*exp(ratio)/1.2d0
write(11,*)
write(11,1103)'Critical superconductivity temperature w/o spin fluctuations:',Tc
Tc = w*exp(sfratio)/1.2d0
write(11,*)
write(11,1103)'Critical superconductivity temperature w spin fluctuations:',Tc
if (verbose) print 2000
return
1000 format(/,'Begin subroutine cpaGG')
1001 format(15(F10.6,1X))
1100 format('This file contains the superconductivity results of the &
  TB-CPA program using the Gaspari-Gyorffy and McMillan theories of &
  superconductivity.',//,A,/)
1101 format(/,'Fermi      # of elec  Total DOS  Fe-s       Fe-p       Fe-&
  xy/xz   Fe-xy      Fe-3z^2    Fe-x^2-y^2 Se/Te-s    Se/Te-p    &
  Se-s       Se-p       Te-s       Te-p')
1102 format(15(F10.6,1X),/)
1103 format(A50,1X,F10.6)
2000 format('End cpaGG',/)
end subroutine cpaGG
