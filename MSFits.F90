!*******************************************************************************
!
!  Routine:     MSFits
!
!  Description: Contains the fitting function to Tout et al. (1996)
!               MNRAS, 281, 257-262, for the luminosity and radius of a
!               main sequence star of solar metallicity. Formulas below can
!               be generalized to metallicities other than solar. See paper
!               for details.
!
!*******************************************************************************

module MSFits

contains

     subroutine calculate_L_main_sequence(M,Lms)
     !  Gives the luminosity of a main sequence star of mass M.
     !  Mass and luminosity are expressed in solar units.
     !
       implicit none
       double precision, parameter   :: alpha = 0.39704170
       double precision, parameter   :: beta = 8.52762600
       double precision, parameter   :: gamm = 0.00025546
       double precision, parameter   :: delta = 5.43288900
       double precision, parameter   :: epsil = 5.56357900
       double precision, parameter   :: zeta = 0.78866060
       double precision, parameter   :: eta = 0.00586685
       double precision, intent(in)  :: M
       double precision, intent(out) :: Lms
       double precision              :: mm
       !double precision, parameter   :: Msun = 1.9889225D33
       !double precision, parameter   :: Lsun = 3.8395D33
     
       mm = M   ! solar masses
       Lms = (alpha*mm**(5.5)+beta*mm**(11))
       Lms = Lms / (gamm + mm**(3) + delta*mm**(5) + epsil*mm**(7) + zeta*mm**(8) + eta*mm**(9.5))
     
     end subroutine calculate_L_main_sequence

     subroutine calculate_R_main_sequence(M,Rms)
     !  Gives the radius of a main sequence star of mass M.
     !  Mass and radius are expressed in solar units.
     !
       implicit none
       double precision, parameter   :: theta = 1.71535900
       double precision, parameter   :: iota = 6.59778800
       double precision, parameter   :: kappa = 10.08855000
       double precision, parameter   :: lambda = 1.01249500
       double precision, parameter   :: mu = 0.07490166
       double precision, parameter   :: nu = 0.01077422
       double precision, parameter   :: xi = 3.08223400
       double precision, parameter   :: omicron = 17.84778000
       double precision, parameter   :: ppi = 0.00022582
       double precision, intent(in)  :: M
       double precision, intent(out) :: Rms
       double precision              :: mm
       !double precision, parameter   :: Msun = 1.9889225E33

       mm = M   ! solar masses
       Rms = (theta*mm**(2.5) + iota*mm**(6.5) + kappa*mm**(11) + lambda*mm**(19) + mu*mm**(19.5))
       Rms = Rms / (nu + xi*mm**(2) + omicron*mm**(8.5) + mm**(18.5) + ppi*mm**(19.5)) 
     end subroutine calculate_R_main_sequence

end module
!===============================================================================
