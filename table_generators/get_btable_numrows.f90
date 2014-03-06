subroutine get_btable_numrows(nrows)
!============================================================================
! PROGRAM: get_btable_numrows(nrows)
!
! DESCRIPTION: Reads the first line from the beta table to calculate the
!     number of rows in the data table.
!
! CALLING SEQUENCE:
!
!       CALL get_btable_numrows(nrows)
!
! WRITTEN: Mikhail Klassen                            (McMaster University)
!          15 / July / 2010
!
! EMAIL: mikhail.klassen@gmail.com
!============================================================================
    implicit none
    integer, intent(out) :: nrows
    integer              :: models,nMs,ios

! Open input file for reading ===============================================
    open(UNIT=100,FILE="beta_table.dat",ACTION="read",&
        & POSITION="REWIND",IOSTAT=ios)
    if (ios /=0) then
        print *,'Error occurred trying to open file beta_table.dat for writing'
        stop
    end if
! Read first line of file ===================================================
    read(UNIT=100,FMT=201),models,nMs

    nrows = models*(nMs-2)

    201 format (2X,I4,2X,I4)

end subroutine get_btable_numrows
