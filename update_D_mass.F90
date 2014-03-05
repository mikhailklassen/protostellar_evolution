subroutine update_D_mass(mdot, dt, Ld, &                        ! Input
                         md )                                   ! In/Out
                         
    implicit none
    double precision, intent(in)    :: mdot, dt, Ld
    double precision, intent(inout) :: md
    double precision, parameter     :: Msun = 1.98892D33         ! in grams
    double precision, parameter     :: Rsun = 6.9550D10          ! in cm
    double precision, parameter     :: Lsun = 3.839D33           ! in erg/s
    double precision, parameter     :: secyr = 31556926.D0       ! seconds in a year
    double precision                :: dmd

    dmd = mdot*dt - 1.D-5 * Msun * (Ld/(15.*Lsun)) * (dt/secyr)

    md = md + dmd

    return

end subroutine update_D_mass
