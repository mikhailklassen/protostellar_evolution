!*******************************************************************************
!
!  Routine:     StarProperties
!
!  Description: Module defining star properties
!
!*******************************************************************************

module StarProperties

      type protostar
          double precision    :: mass     ! mass
          double precision    :: mdot     ! accretion rate
          double precision    :: radius   ! stellar radius
          double precision    :: polyn    ! polytropic index
          double precision    :: mdeut    ! deuterium mass
          double precision    :: lint     ! intrinsic luminosity
          double precision    :: lum      ! total luminosity (incl. accretion lum)
          integer :: stage    ! evolutionary stage (0, 5)
      end type protostar

      type(protostar), save :: star
   
      double precision, parameter :: Msun = 1.98892D33         ! in grams
      double precision, parameter :: Rsun = 6.9550D10          ! in cm
      double precision, parameter :: secyr = 31556926.D0       ! seconds in a year

end module StarProperties
