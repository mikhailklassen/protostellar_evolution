subroutine EvolveProtostellar(dt)

   use protostellar_interface, only: star, Msun, Rsun, secyr

   implicit none

   double precision, intent(in) :: dt

   integer          :: i, pno, nrOfStarsT, nrOfStarsL
   double precision :: mass, md, mdot, n, r
   double precision :: mmdot, dm
   integer          :: stage
   double precision :: Ld, Lint, Li
   double precision :: Lms, Lacc, Ldisk, Lh
   double precision :: Ltot
   double precision :: ag, beta, dlogbeta

   ! Get the star properties
   mass  = star%mass   
   mdot  = star%mdot   
   r     = star%radius 
   n     = star%polyn  
   md    = star%mdeut  
   Lint  = star%lint   
   Ltot  = star%lum    
   stage = star%stage
   ! If the the sink particle mass is above threshold mass, begin evolution 
   if (mass .GE. 0.1*Msun) then

     ! Convert mdot, which has units of g/s, into Msun/yr
     mmdot = mdot * secyr / Msun

     ! If sink particle ready to begin evolving, advance to next stage
     ! and initialize characteristic sink particle quantities
     if ( stage .EQ. 0 ) then
         stage = 1
         r  = 2.5 * Rsun * (mmdot / 1.E-5)**(0.2)
         n  = 5 - 3 / (1.475 + 0.07 * log10(mmdot))
         if ( n .LT. 1.5) n = 1.5
         if ( n .GT. 3.0) n = 3.0
         md = mass
     endif

     ! Calculate intermediate values
     call get_derived_values(stage, mass, mdot, r, n, &             ! Input
                             Ld, Lint, Li, Lms, beta, dlogbeta, ag) ! Return
    
     ! (1) Update the radius and deuterium mass
     call update_radius(mass, mdot, ag, beta, dlogbeta, dt, &       ! Input
                        Ld, Lint, Li, stage, &                      ! Input
                        r )                                         ! In/Out

     call update_D_mass(mdot, dt, Ld, &                             ! Input
                        md )                                        ! In/Out

     ! (2) Compute the new luminosity
     call update_luminosity(Lint, mass, mdot, r, &                  ! Input
                            Ltot )                                  ! Return
     
     ! (3) Advance to next evolutionary stage
     call update_stage(mass, md, Ld, Lms, &                         ! Input
                       n, r, stage)                                 ! In/Out

     
     ! Store dummy values into particle list
     star%mass   = mass
     star%mdot   = mdot
     star%radius = r
     star%polyn  = n
     star%mdeut  = md
     star%lint   = Lint
     star%lum    = Ltot
     star%stage  = stage

   endif ! mass>0.1Msun

end subroutine EvolveProtostellar
