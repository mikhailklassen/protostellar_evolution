program test
    use polytools
    implicit none
    integer :: cols,i,j
    integer :: models,nMs,mRows
    double precision, allocatable, dimension(:)   :: mvec
    double precision, allocatable, dimension(:)   :: ns
    double precision, allocatable, dimension(:,:) :: beta_table
    double precision, allocatable, dimension(:,:) :: dlogbeta_table
    double precision :: beta,dlogbeta,mass,n

    call getNumModels(models)
    nMs    = 50
    mRows  = nMs - 2

    cols = models * mRows

    allocate(beta_table(models,mRows))
    allocate(dlogbeta_table(models,mRows))
    allocate(ns(models))

    call get_model_table_cols(models,n=ns)

    call get_beta_table_mvec(mvec)

    call get_beta_table(cols,models,mRows,beta_table,dlogbeta_table)

    print *,beta_table

    mass = 10.D0
    n = 3.0D0

    call interp2dlin(beta_table,shape(beta_table),ns,mvec,n,mass,beta)
    call interp2dlin(dlogbeta_table,shape(dlogbeta_table),ns,mvec,n,mass,dlogbeta)

    ! Testing: Print a small table with interpolated output
    !print *,'------------------------------------------------'
    !print *,'|   Beta               |   dlog(beta)/dlog(m)  |'
    !print *,'|----------------------|-----------------------|'
    !print '(1X,A1,f20.16,2X,A1,1X,f20.16,2X,A1)','|',beta,'|',dlogbeta,'|'
    !print *,'------------------------------------------------'

end program test

subroutine get_beta_table(cols,models,mRows,beta_table,dlogbeta_table)
!============================================================================
! PROGRAM: get_beta_table.f90
!
! DESCRIPTION:
!
! WRITTEN: Mikhail Klassen                             (McMaster University)
!          15 / July / 2010
!
! EMAIL: mikhail.klassen@gmail.com
!============================================================================
    implicit none
    integer,intent(in) :: cols,models,mRows
    integer            :: nMs,tmp_models
    integer            :: i,ii,j,ios
    double precision,dimension(cols,3) :: betadat
    double precision,dimension(models,mRows),intent(out) :: beta_table
    double precision,dimension(models,mRows),intent(out) :: dlogbeta_table

! OPEN INPUT FILE FOR READING ===============================================
    open(UNIT=100,FILE="beta_table.dat",ACTION="read",&
        & POSITION="REWIND",IOSTAT=ios)
    if (ios /=0) then
        print *,'Error occurred trying to open file beta_table.dat for writing'
        stop
    end if
! READ FILE =================================================================
    read(UNIT=100,FMT=201) tmp_models,nMs
    if (cols /= models*(nMs-2)) then
        print *,'Error: mismatched dimensions in expected size of beta_table'
        print *,'Expected number of cols: ',cols
        print *,'Actualy number of cols: ',models*(nMs-2)
        stop
    endif

    do i=1,cols
        read(UNIT=100,FMT=200,IOSTAT=ios) betadat(i,:)
        if (ios<0) EXIT ! Check to see if end-of-file reached
    end do

! RECONFIGURE DATA ==========================================================
    ii = 1
    do i=1,cols,mRows
        do j=1,mRows
            beta_table(ii,j)=betadat(i+j-1,2)
            dlogbeta_table(ii,j)=betadat(i+j-1,3)
        end do
        ii = ii + 1
    end do

! CLOSE FILE ================================================================
    close (UNIT=100)


! FORMAT STATEMENTS =========================================================
    200 format (2X,f20.16,5X,f18.16,5X,f19.16)
    201 format (2X,I4,2X,I4)

    return

end subroutine get_beta_table
