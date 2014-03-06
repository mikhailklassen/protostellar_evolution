program make_beta_table
!=======================================================================================
! PROGRAM: make_beta_table.f90
!
! DESCRIPTION: Generates a lookup table and writes it to the file 'BetaTable.dat'. The
!     table consists of three columns (mass, beta, d\log(beta)/d\log(m)) for polytropic
!     stars. The volume-averaged beta is computed for a range of stellar masses and then
!     the derivate d\log(beta)/d\log(m) is computed numerically using the central
!     difference method.
!
!     This whole process is repeated once for each polytropic index contained in another
!     precomputed table, 'modeldata_table.dat'. If modeldata_table contains P number of 
!     models, and the number of star masses for which beta is calculated is Q, then the 
!     number of rows written to the 'BetaTable.dat' is P*(Q-2). The reduction in P by 
!     two is due to the trimming of the first and last rows in each section, which was 
!     necessary because of our use of the central difference method for computing the 
!     derivatives earlier.
!
! WRITTEN: Mikhail Klassen                                         (McMaster University)
!          15 / July / 2010
!
! UPDATED: Mikhail Klassen
!          08 / Jan / 2011
!          
!          - Added support for calculating d\log(beta/beta_c)/d\log(m)
!
! EMAIL: mikhail.klassen@gmail.com
!=======================================================================================
! Modules to include ===================================================================
    use astro_constants
    use polytools
! Variable declarations ================================================================
    implicit none
    double precision, parameter :: hmax = 1.D-5            ! Maximum step size in solver
    double precision, parameter :: tol  = 1.D-8            ! Maximum error tolerance of solver
    double precision :: n, M, Bn, beta, beta_c
    double precision, allocatable, dimension(:) :: Ms      ! Mass array
    integer, parameter                          :: nMs=50  ! Mass array granularity
    double precision, allocatable, dimension(:) :: betas   ! Temp array for storing betas
    double precision, allocatable, dimension(:) :: beta_cs ! Temp array for storing central betas
    double precision, allocatable, dimension(:) :: dlogBdlogM,dlogBBcdlogM
    double precision, allocatable, dimension(:) :: ns      ! Polytropic index array
    double precision, allocatable, dimension(:) :: Bns     ! Pressure coeff. array
    integer                                     :: models  ! Array granularity
    integer                                     :: i,j     ! Loop counters
    integer                                     :: ios     ! I/O status flag

! Allocate memory ======================================================================
    allocate(Ms(nMs))                                   ! Allocate memory for mass array
    allocate(betas(nMs))
    allocate(beta_cs(nMs))
    allocate(dlogBdlogM(nMs))
    allocate(dlogBBcdlogM(nMs))
    dlogBdlogM(1) = 0.D0
    dlogBdlogM(nMs) = 0.D0
    dlogBBcdlogM(1) = 0.D0
    dlogBBcdlogM(nMs) = 0.D0
! Look up other array granularities and allocate memory ================================
    call getNumModels(models)
    allocate(ns(models))                                ! Allocate memory for n array
    allocate(Bns(models))                               ! Allocate memory for Bn array
! Extract columns from modeldata_table.dat =============================================
    call get_model_table_cols(models,n=ns,Bn=Bns)
! Set mass array =======================================================================
    Ms = logspace(-2.3D0,2.6D0,nMs)   ! Create an array of logarithmically-spaced points
    Ms = Ms * Msun                    ! between 0.01 solar masses and ~250 solar masses

! Open file for writing ================================================================
    open(UNIT=100,FILE="beta_table.dat",STATUS="REPLACE",ACTION="WRITE",&
         & POSITION="REWIND",IOSTAT=ios)
    if (ios /=0) then
        print *,'Error occurred trying to open file BetaTable.dat for writing'
        stop
    end if
! Write the file header ================================================================
    write (UNIT=100,FMT=201,iostat=ios),models,nMs
    if (ios /= 0) then
        print *,'Error occurred trying to write header to file BetaTable.dat'
        stop
    end if
! Write to file ========================================================================
    do j = 1, models
        !write (UNIT=100,FMT=*,iostat=ios),'n = ',ns(j)
        !if (ios /= 0) then; print *,'Error writing to file.'; stop; endif
        
        do i = 1, nMs
            call beta_func_solver(ns(j),Ms(i),Bns(j),1.D0,beta_c)
            call solve_beta(ns(j),hmax,tol,Ms(i),Bns(j),beta)
            beta_cs(i) = beta_c
            betas(i) = beta
            if (i > 2 ) then
                dlogBdlogM(i-1) = (log10(betas(i))-log10(betas(i-2)))/(log10(Ms(i))-log10(Ms(i-2)))
                dlogBBcdlogM(i-1) = (log10(betas(i)/beta_cs(i))-log10(betas(i-2)/beta_cs(i-2)))/(log10(Ms(i))-log10(Ms(i-2)))
                write (UNIT=100,FMT=200,iostat=ios),Ms(i-1)/Msun,betas(i-1),dlogBdlogM(i-1),dlogBBcdlogM(i-1)
                if (ios /= 0) then; print *,'Error writing to file.'; stop; endif
            end if
        end do
    end do

! Close file ===========================================================================
    close (UNIT=100)

! Format Statements ====================================================================
    200 format (2X,f20.16,5X,f18.16,5X,f19.16,5X,f19.16)
    201 format (2X,I4,2X,I4)

end program make_beta_table

