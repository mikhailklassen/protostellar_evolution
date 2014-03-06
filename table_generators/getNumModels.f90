subroutine getNumModels(models)
    implicit none
    integer :: models,ios

    open(UNIT=100,FILE="MODELDATA",ACTION="read",&
         & POSITION="REWIND",IOSTAT=ios)
    if (ios /=0) then
        print *,'Error occurred trying to open file MODELDATA for writing'
        stop
    end if
    read(UNIT=100,FMT=*),models
    close (UNIT=100)
    return
end subroutine getNumModels