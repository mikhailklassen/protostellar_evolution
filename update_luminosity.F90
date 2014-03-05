subroutine update_luminosity(Lint, mass, mdot, r, &             ! Input
                             Ltot )                             ! Return
    implicit none
    double precision, intent(in)   :: Lint, mass, mdot, r
    double precision, intent(out)  :: Ltot
    double precision               :: Lacc, Ldisk
    double precision, parameter    :: fk=0.5
    double precision, parameter    :: facc=0.5
    double precision, parameter    :: G=6.674E-8

    Lacc = facc* fk * G * mass * mdot / r

    Ldisk = (1. - fk) * G * mass * mdot / r
   
    Ltot = Lint + Lacc + Ldisk

    return

end subroutine update_luminosity
