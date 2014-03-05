subroutine get_beta_table(cols,models,mRows,beta_table,dlogbeta_table,dlogbbc_table)
!============================================================================
! PROGRAM: get_beta_table.f90
!
! DESCRIPTION:
!
! WRITTEN: Mikhail Klassen                             (McMaster University)
!          15 / July / 2010
!
! UPDATED: Mikhail Klassen
!          10 / Jan / 2011
!
!          - Added support for creating a 2D array of values for
!            d\log(beta/beta_c)/d\log(m) as a function of polytropic index n
!            and stellar mass m. This array is returned as dlogbbc_table
!
! EMAIL: mikhail.klassen@gmail.com
!============================================================================
    implicit none
    integer,intent(in)                                   :: cols,models,mRows
    integer                                              :: nMs,tmp_models
    integer                                              :: i,ii,j,ios
    double precision,dimension(cols,4)                   :: betadat
    double precision,dimension(models,mRows),intent(out) :: beta_table
    double precision,dimension(models,mRows),intent(out) :: dlogbeta_table
    double precision,dimension(models,mRows),intent(out) :: dlogbbc_table

! OPEN INPUT FILE FOR READING ===============================================
    open(UNIT=100,FILE="beta_table.dat",ACTION="read",&
        & POSITION="REWIND",IOSTAT=ios)
    if (ios /=0) then
        print *,'Error occurred trying to open file beta_table.dat for reading'
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
            dlogbbc_table(ii,j)=betadat(i+j-1,4)
        end do
        ii = ii + 1
    end do

! CLOSE FILE ================================================================
    close (UNIT=100)


! FORMAT STATEMENTS =========================================================
    200 format (2X,f20.16,5X,f18.16,5X,f19.16,5X,f19.16)
    201 format (2X,I4,2X,I4)

    return

end subroutine get_beta_table
