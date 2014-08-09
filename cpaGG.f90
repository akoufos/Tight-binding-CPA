subroutine cpaGG(fermi, nelec, nef, N)
!--------------------------------------------------------------------------
! This subroutine will read the CPA DOS at the Fermi level to approximate
! the eta, Hopfield, parameter for use with the Gaspari-Gyorffy theory.
! Using the eta parameter, the subroutine will then calculate the
! electron-phonon coupling constant to use in the McMillian equation.
! Based on Papaconstantopoulos and Pickett, PRB 31, 7093 (1985)
!--------------------------------------------------------------------------
! Variables:
! C(ntype) - Constant for calculation of etas
! effl - Total electron-phonon coupling constant
! effM - Total effective atomic mass 
! effM2 - Effective atomic mass of Se/Te contribution
! effw2 - Effective average squared phonon frequency
! eta(ntype) - Hopfield parameter for each atom type (eV/ang^2)
! fermi - Fermi level
! lambda(ntype) - Electron-phonon coupling constant for each atom type
! M(3) - Atomic masses, in atomic units
! mode - Subroutine run mode to calculate superconductivity properties
!   1/3 - Calculate during CPA program
!   2 - Calculate at a later date of CPA completion
! mt1/2 - Muffin-tin ratio used to correct for lack of muffin-tins in the
!   TB-CPA formulism used (Version 3+: Not used for GG calculations.)
! N(12) - Angular momentum components of the DOS (Ry)
! nef - Total DOS at the Fermi level (Ry)
! nelec - Number of electrons
! Tc - Critical superconducting temperature of the system (K)
! Theta(ntype) - Debye temperature of the systems
! w2(ntype) - Squared average phonon frequency of the atom types
!--------------------------------------------------------------------------
! Units:
! [eta] = eV/angstroms^2
! [lambda] = unitless
! [M] = u (atomic mass)
! [w2] = K^2
! [Tc] = K
! [Theta] = K
!--------------------------------------------------------------------------
! Notes:
!
! V1 - 05-23-2014 - Subroutine only considers the major contributions to 
! the DOS at the Fermi level when calculating the Hopfield parameter. For 
! the FeSeTe system, those components are the Fe-d DOS and the Se/Te-p
! DOS.
!
! V2 - 06-09-2014 - Added more information for the input file, allowing for
! better approximations of the effective quanities used as final inputs for
! the Gaspari-Gyorffy method.
!
! V3 - 06-17-2014 - Changed the equations for calculating etas. Now have
! etas calculated for sp and pd components. The constant used, now come
! from the sp, pd, df components of eta calculated with LAPW results. This
! new constant is calculated by the equation 
!          eta = C*sum[N_l*N_{l+1}/N_t]
! where N_l is the decomposed DOS and N_t is the total DOS, all at
! the Fermi level. C is determined for l=1,3, where the CPA equation for
! eta only uses l=1,2 (assumes f-DOS is zero). CPA equation is given as
!          eta = C(N_s*N_p/N_t + N_p*N_d/N_t)
!--------------------------------------------------------------------------
use global
use unitconvert
use concentration
implicit none
integer(kind=4) :: i
real(kind=8) :: eta(3), lambda(3), Tc, M(3), w2(ntype), effw2, &
  effM, effM2, effl, ratio, sfratio, w, lsf, must, C(2*ntype), effTheta, &
  Theta(ntype), Theta2, mt, dos(12), eta1(2), eta2(2), eta3(2), C1, C2, C3
real(kind=8), intent(inout) :: fermi, nelec, nef, N(12)
real(kind=8), parameter :: mt1 = 1.d0/26.920d0
real(kind=8), parameter :: mt2 = 1.d0/36.4143d0
character(len=2) :: atom(ntype+1)
if (verbose) print 1000
open(11,file='CPA-GG.out')
select case (mode)
  case (1,3)
    write(6,*)'Continuing superconductivity calculations with results'
  case (2)
    open(10,file='dosapw.itp')
    read(10,*)
    read(10,1001)fermi,nelec,NEF,(N(i),i=1,12)
    close(10)
  case default
    write(6,*)'Problem with mode. Nothing was done'
    return ! Return to mainn subroutine
end select
write(11,1100)title
write(11,1101)
write(11,1102)fermi,nelec,2.0d0*nef,(2.0d0*N(i),i=1,12)
mt = mt2 + con*(mt1-mt2)
! Correct decomposed DOS to match LAPW results
dos(1) = N(1)*mt
dos(2) = N(2)*mt
dos(3:6) = N(3:6)
dos(7:) = N(7:)*mt
write(11,'(A)')'LAPW like DOS'
write(11,1101)
write(11,1102)fermi,nelec,2.0d0*nef,(2.0d0*dos(i),i=1,12)
open(10,file='superCPA.in')
read(10,*) (atom(i),i=1,3)
read(10,*) (C(i),i=1,4)
read(10,*) (M(i),i=1,3)
read(10,*) Theta(1), Theta(2)
read(10,*) must
read(10,*) lsf
close(10)
effTheta = Theta(2) + con*(Theta(1)-Theta(2))
Theta2 = effTheta**2.0d0/2.0d0
! Calculate the etas using LAPW constants and CPA DOS. 
!!!!!! No need for muffin-
!!!!!! tin corrections.
C1 = C(3) + con*(C(1)-C(3)) ! atom type 1 eta constant
!C2 = C(4) + con*(C(2)-C(4)) ! atom type 2 eta constant
C2 = C(2) ! atom type 2 eta constant
C3 = C(4) ! atom type 3 eta constant
eta1(1) = C1*(dos(1)*dos(2)/nef)
eta1(2) = C1*(dos(2)*sum(dos(3:6))/nef)
eta2(1) = C2*(dos(9)*dos(10)/nef)
eta2(2) = C2*(dos(10)*0.1d0/nef)
eta3(1) = C3*(dos(11)*dos(12)/nef)
eta3(2) = C3*(dos(12)*0.1d0/nef)
eta(1) = sum(eta1)
eta(2) = sum(eta2)
eta(3) = sum(eta3)
write(6,1104)'Eta '//atom(1)//' sp and pd: ',eta1
write(6,1104)'Eta '//atom(2)//' sp and pd: ',eta2
write(6,1104)'Eta '//atom(3)//' sp and pd: ',eta3
write(11,1104)'Eta '//atom(1)//' sp and pd: ',eta1
write(11,1104)'Eta '//atom(2)//' sp and pd: ',eta2
write(11,1104)'Eta '//atom(3)//' sp and pd: ',eta3
write(11,*)
write(11,1103)'Hopfield parameter of '//atom(1)//' (eV/angstrom^2):',eta(1)
write(11,1103)'Hopfield parameter of '//atom(2)//' (eV/angstrom^2):',eta(2)
write(11,1103)'Hopfield parameter of '//atom(3)//' (eV/angstrom^2):',eta(3)
!write(11,1103)'Total eta (eV/angstrom^2):',sum(eta)
! Calculate the effective atomic mass of the system
effM2 = M(3) + con*(M(2)-M(3))
effM = M(1) + effM2
write(11,*)
write(11,1103)'Effective Atomic Mass of '//atom(2)//'/'//atom(3)//' (u):', &
  effM2
write(11,1103)'Total Effective Atomic Mass (u):',effM
write(11,*)
write(11,1103)'Debye temperature of '//atom(1)//atom(2)//' (K^2):',Theta(1)
write(11,1103)'Debye temperature of '//atom(1)//atom(3)//' (K^2):',Theta(2)
write(11,1103)'Effective Debye temperature of system (K^2):',effTheta
! Calculate the average sqaured phonon frequencies and effective values
w2(:) = Theta(:)**2.0d0/2.0d0*K2meV*K2meV
effw2 = w2(2) + con*(w2(1)-w2(2))
w = sqrt(effw2)
write(11,*)
write(11,1103)'Average squared phonon frequency of '//atom(1)//atom(2)// &
  ' (meV^2):',w2(1)
write(11,1103)'Average squared phonon frequency of '//atom(1)//atom(3)// &
  ' (meV^2):',w2(2)
write(11,1103)'Effective squared average phonon frequency (meV^2):',effw2
write(11,1103)'Effective phonon frequency (meV):',w
write(11,*)
! Calculate the electron-phonon coupling constant
!lambda(1) = eta(1)/(M(1)*effw2*u2kg*meV2Hz*meV2Hz*kg_ss2eV_AA)
!lambda(2) = eta(2)/(effM2*effw2*u2kg*meV2Hz*meV2Hz*kg_ss2eV_AA)
lambda(1) = eta(1)/(M(1)*Theta2*uKK2eV_AA)
lambda(2) = eta(2)/(M(2)*Theta2*uKK2eV_AA)
lambda(3) = eta(3)/(M(3)*Theta2*uKK2eV_AA)
effl = sum(lambda)
print*,(u2kg*K2eV*K2eV*eV2Hz*eV2Hz*kg_ss2eV_AA)
print*,(M(1)*effw2*u2kg*meV2Hz*meV2Hz*kg_ss2eV_AA)
print*,(effM2*effw2*u2kg*meV2Hz*meV2Hz*kg_ss2eV_AA)
print*,(M(1)*Theta(1)**2.0d0/2.0d0*uKK2eV_AA)
print*,(effM2*Theta(1)**2.0d0/2.0d0*uKK2eV_AA)
print*,(M(1)*Theta(1)**2.0d0/2.0d0*uKK2eV_AA2)
print*,(effM2*Theta(1)**2.0d0/2.0d0*uKK2eV_AA2)
write(11,1103)'Electron-phonon coupling constant of '//atom(1)//':', &
  lambda(1)
write(11,1103)'Electron-phonon coupling constant of '//atom(2)//':', &
  lambda(2)
write(11,1103)'Electron-phonon coupling constant of '//atom(3)//':', &
  lambda(3)
write(11,1103)'Total electron-phonon coupling constant:',effl
! Calculate the critical superconducting temperature
ratio = -1.04d0*(1.0d0+effl) / (effl-must*(1.0d0+0.62d0*effl))
sfratio = -1.04d0*(1.0d0+effl+lsf) / &
  (effl-must*(1.0d0+0.62d0*(effl+lsf))-lsf)
Tc = effTheta*exp(ratio)/1.45d0
write(11,*)
write(11,1103)'Critical superconductivity temperature:'
write(11,1103)'without spin fluctuations:',Tc
Tc = effTheta*exp(sfratio)/1.45d0
write(11,*)
write(11,1103)'with spin fluctuations:',Tc
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
1103 format(A55,1X,F10.6)
1104 format(A20,(2(1X,F10.6)))
2000 format('End cpaGG',/)
end subroutine cpaGG
