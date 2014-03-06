subroutine get_derived_values(stage, mass, mdot, r, n, &             ! Input
                              Ld, Lint, Li, Lms, beta, dlogbeta, ag) ! Return
      
   use protostellar_interface, only: calculate_L_main_sequence, &
                        getNumModels, &
                        get_model_table_cols, &
                        get_beta_table_mvec, &
                        Msun, Rsun, Lsun, secyr, &
                        sigma, G, pi

   implicit none
   integer, intent(in)           :: stage
   double precision, intent(in)  :: mass,mdot,r,n
   double precision, intent(out) :: Ld,Lint,Li
   double precision              :: Lh,Lms
   double precision, intent(out) :: beta,dlogbeta,ag
   double precision              :: dlogbbc
   double precision, parameter   :: fk = 0.5
   double precision              :: mmdot
   ! Necessary variables for retrieving beta and dlogbeta
   integer                                     :: cols
   integer                                     :: models,nMs,mRows
   double precision,allocatable,dimension(:)   :: mvec
   double precision,allocatable,dimension(:)   :: ns
   double precision,allocatable,dimension(:,:) :: beta_table
   double precision,allocatable,dimension(:,:) :: dlogbeta_table
   double precision,allocatable,dimension(:,:) :: dlogbbc_table

   ! Convert accretion rate from g/s to Msun/yr
   mmdot = mdot * secyr / Msun

   ! Calculate coefficient of gravitational binding energy for a polytrop
   ag    = 3. / (5. - n)

   ! Lookup the polytropic beta and dlog(beta)/dlog(m)
   call getNumModels(models)
   nMs = 50
   mRows = nMs - 2
   cols = models * mRows
   allocate(beta_table(models,mRows))
   allocate(dlogbeta_table(models,mRows))
   allocate(dlogbbc_table(models,mRows))
   allocate(ns(models))
   
   call get_model_table_cols(models,n=ns)
   call get_beta_table_mvec(mvec)
   call get_beta_table(cols,models,mRows,beta_table,dlogbeta_table,dlogbbc_table)

   ! Interpolate within beta tables to get beta and dlogbeta given our mass and
   ! polytropic index n
   call interp2dlin(beta_table,shape(beta_table),ns,mvec,n,mass/Msun,beta)
   call interp2dlin(dlogbeta_table,shape(dlogbeta_table),ns,mvec,n,mass/Msun,dlogbeta)
   call interp2dlin(dlogbbc_table,shape(dlogbbc_table),ns,mvec,n,mass/Msun,dlogbbc)
   
   ! Evaluate intrinsic stellar luminosity
   call calculate_L_main_sequence(mass/Msun,Lms)
   Lms = Lms * Lsun
   if ( stage .EQ. 5 ) then
       Lint = Lms
   else
       Lh  = 4. * pi * r * r * sigma * (3000.)**4
       Lint = max(Lms,Lh)
   endif

   ! Evaluate ionization luminosity
   Li = 2.5 * Lsun * mmdot / 1.E-5

   ! Evaluate deuterium burning luminosity
   ! N.B.: dlogbbc represents the term d\log(beta/beta_c)/d\log(m)
   if (stage .LT. 2) then
      Ld = 0.
   else if (stage .EQ. 2) then
      Ld = Lint + Li + (G * mass)/r * mdot * &
           ( 1. - fk - ag*beta/2. * (1. + dlogbbc))
   else if (stage .GT. 2) then
      Ld = 15. * Lsun * mmdot / 1.E-5
   endif

   return

end subroutine get_derived_values
