program make_modeldata
!============================================================================
! PROGRAM: make_modeldata.f90
!
! DESCRIPTION: Generates a lookup table of polytrope parameters for
!		polytropes of different polytropic index and writes the table to a
!		file.
!
! WRITTEN: Mikhail Klassen                             (McMaster University)
!          9 / June / 2010
!
! EMAIL: mikhail.klassen@gmail.com
!============================================================================
    implicit none
! VARIABLE DECLARATIONS =====================================================
    double precision			 				:: n,hmax,tol,xi,Bn
    double precision							:: Dn,PcPbar,mx2dthetadx,zfin
    double precision,allocatable,dimension(:)	:: nn
    integer                                     :: i,ios
                                                ! SolveLE won't print any lines
    integer,parameter                           :: iout=99
    integer,parameter                           :: models=100
	double precision                            :: nmin,nmax
	double precision,dimension(models,6)        :: MODELDATA
! INTERACES TO EXTERNAL FUNCTIONS ===========================================
	interface
		function LINSPACE(A,B,Ns)
		double precision               :: A,B
		integer                        :: Ns
		double precision,dimension(Ns) :: linspace
		end function
	end interface
! ===========================================================================
    i = 0
	hmax = 1.D-1			! Sets the maximum step size for RFK45 method
    tol = 1.D-8				! Sets the maximum error tolerance

	allocate(nn(models))	! Allocate memory for array of polytropic indices
	nmin = 1.500D0			! Smallest index
	nmax = 3.000D0			! Largest index
	! Fill an array with equally spaced polytropic index values between the
	!  bounds set above
	nn = linspace( nmin, nmax, models )

!	Solve the polytrope defined by each of the polytropic indices in the
!   array nn
	do i = 1,models
		n = nn(i)
    	call solveLE(n,iout,hmax,tol,xi,Dn,mx2dthetadx,zfin)

    	Bn = (3.D0**((6.D0-2.D0*n)/(3.D0*n)))/(n+1.D0) * xi**(-(2.D0*n+6.D0)/(3.D0*n))
    	Bn = Bn * abs(zfin)**((6.D0-4.D0*n)/(3.D0*n))

		MODELDATA(i,1) = n			 !;PRINT *,'Index n: ',n
    	MODELDATA(i,2) = xi          !;PRINT *,'Radius: ',xi
    	MODELDATA(i,3) = zfin        !;PRINT *,'d(theta)/dx: ',zfin
    	MODELDATA(i,4) = Dn          !;PRINT *,'rho_c/rho_bar: ',Dn
    	MODELDATA(i,5) = Bn          !;PRINT *,'Bn: ',Bn
    	MODELDATA(i,6) = mx2dthetadx !;PRINT *,'-x^2 d(theta)/dx: ',mx2dthetadx
    end do

! OPEN INPUT FILE FOR WRITING ===============================================
	open(UNIT=100,FILE="modeldata_table.dat",STATUS="REPLACE",ACTION="WRITE",&
	     & POSITION="REWIND",IOSTAT=ios)
	if (ios /=0) THEN
		print *,'Error occurred trying to open file MODELDATA for writing'
		stop
	end if
! WRITE TO FILE =============================================================
	write (UNIT=100,FMT=*),models
	write (UNIT=100,FMT=200) &
		& 'Index n','Radius','d(theta)/dx',&
		& 'rho_c/rho_bar','Bn','-x^2 d(theta)/dx'
	do i = 1,models
		write (UNIT=100,FMT=201),MODELDATA(i,:)
	end do
	endfile 100
! CLOSE FILE ================================================================
	close (UNIT=100)

! FORMAT STATEMENTS =========================================================
	199 format (1X,A,2X,I10)
	200 format (2X,A,8X,A,8X,A,4X,A,6X,A,8X,A)
	201 format (1X,F12.9,3X,F12.9,3X,F12.9,3X,F12.9,3X,F13.9,3X,F12.9)
! ===========================================================================

	deallocate(nn)			! Deallocate the memory held by nn

end program make_modeldata


! EXTERNAL FUNCTIONS ========================================================
function linspace(A,B,Ns)
    ! Generate an array of Ns equally spaced grid points between A and B
    implicit none
    double precision               :: A,B
    integer                        :: Ns
    integer                        :: i
    double precision,dimension(Ns) :: LINSPACE

    linspace = (/ (i,i=0,Ns-1) /)
    linspace = linspace / (Ns-1.) * (B-A)+A
    return
end function linspace
! END OF LIST OF EXTERNAL FUNCTIONS =========================================
