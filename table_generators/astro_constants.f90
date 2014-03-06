module astro_constants
    implicit none
    save
	! in CGS units
    double precision,parameter :: Msun = 1.98892D33                !g
    double precision,parameter :: Rsun = 6.955D10                  !cm
    double precision,parameter :: Lsun = 3.839D33                  !erg/s
    double precision,parameter :: G = 6.67428D-8                   !cm^3/g/s^2
    double precision,parameter :: planck = 6.62606896D-27          !erg
    double precision,parameter :: speed_of_light = 29979245800.0D0 !cm/s
    double precision,parameter :: stefan_boltzmann = 5.670400D-5   !erg/cm^2/s/K^4
    double precision,parameter :: kb = 1.3806504D-16               !erg/K
    double precision,parameter :: pi = 3.14159265358979D0
    double precision,parameter :: secyr = 31556926.0D0             !seconds in a year
    double precision,parameter :: mu = 0.6130000D0                 !mean molecular weight
    double precision,parameter :: H = 1.673534D-24                 !g
    double precision,parameter :: a_rad = 7.56591D-15              !erg/cm^3/K^4
    double precision,parameter :: f_k = 0.50000000D0


end module
