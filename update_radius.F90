subroutine update_radius(mass, mdot, ag, beta, dlogbeta, dt, &  ! Input
                         Ld, Lint, Li, stage, &                 ! Input
                         r )                                    ! In/Out
   
    use MSFits

    implicit none
    double precision, intent(in)    :: mass, mdot, ag, beta, dlogbeta, dt
    double precision, intent(in)    :: Ld, Lint, Li
    integer, intent(in) :: stage
    double precision, intent(inout) :: r
    double precision                :: dr,dm
    double precision                :: Rms
    double precision, parameter     :: G=6.674E-8
    double precision, parameter     :: fk=0.5
    double precision, parameter     :: Msun = 1.98892E33
    double precision, parameter     :: Rsun = 6.955E10
    
   
    if ( stage .EQ. 5 ) then
        call calculate_R_main_sequence(mass/Msun,Rms)
        r = Rms * Rsun
    else
        dm = mdot*dt
        dr = 2. * (dm/mass) * (1. - (1. - fk)/(ag*beta) + 0.5*dlogbeta) * r &
           - 2. * (dt/(ag*beta)) * (r/(G*mass*mass))*(Lint + Li - Ld) * r

        r = r + dr
    endif

    return 

end subroutine update_radius
