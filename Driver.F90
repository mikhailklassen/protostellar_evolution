!*******************************************************************************
!
!  Routine:     EvolveProtostellar
!
!  Description: Evolution subroutine for the protostellar_evolution sub-module.
!               Evolves the radius and luminosity of the sink particles so that
!               radiation codes can accurately determine feedback.
!
!  Written:     Mikhail Klassen 2010 (McMaster University)
!  Updated:     Mikhail Klassen 2011 (McMaster University)
!  ...    :     Mikhail Klassen 2014 (McMaster University)
!
!*******************************************************************************

program Driver

    use protostellar_interface, only: star, Msun, Rsun, secyr

    implicit none

    double precision :: dt
    double precision :: mdot
    double precision :: time, maxtime

    ! Initialize our star
    star%mass   = 0.0
    star%mdot   = 0.0
    star%radius = 0.0
    star%polyn  = 0.0
    star%mdeut  = 0.0
    star%lint   = 0.0
    star%lum    = 0.0
    star%stage  = 0

    ! Set the accretion rate (grams per second)
    mdot = 1.0D-4 * Msun / secyr
    star%mdot = mdot

    ! Set the timestep (seconds)
    dt = 100.0 * secyr

    ! Initialize the simulation time and set the maxtime
    time = 0.0
    maxtime = 1.D5 * secyr

    ! Evolve the star
    do while(time < maxtime)
        call EvolveProtostellar(dt)
        star%mass = star%mass + star%mdot*dt
        time = time + dt
        write(*,*) star%mass, star%mdot, star%radius, &
                        & star%polyn, star%mdeut, star%lint, &
                        & star%lum, star%stage
    end do

    print *, 'Simulation reach max time.'

end program Driver
!'(7G12.5,I3)'
