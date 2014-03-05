subroutine interp2dlin(arr,arrsize,xvec,yvec,x0,y0,ans)
!============================================================================
! PROGRAM: interp2dlinear
!
! DESCRIPTION: Performs a bilinear interpolation on the data in the 2D array
!       'arr', which is considered a function of the independent variables
!       x and y. The x and y vectors must also be passed to the array, as
!       well as the point (x0,y0) where the function is to be evaluated.
!       Value is returned as variable 'ans'.
!
! CALLING SEQUENCE:
!
!       CALL interp2dlinear(arr,arrsize,xvec,yvec,x0,y0,ans)
!
! WRITTEN: Mikhail Klassen                            (McMaster University)
!          5 / July / 2010
!
! MODIFIED: Mikhail Klassen                           (McMaster University)
!          21 / July / 2010
!
! EMAIL: mikhail.klassen@gmail.com
!============================================================================
    implicit none
    integer,intent(in)           :: arrsize(2)
    double precision,intent(in)  :: arr(arrsize(1),arrsize(2))
    double precision,intent(in)  :: xvec(arrsize(1))
    double precision,intent(in)  :: yvec(arrsize(2))
    double precision,intent(in)  :: x0,y0
    double precision,intent(out) :: ans
    double precision             :: tmpshape(2)
    double precision             :: tmpm(arrsize(1)),tmpn(arrsize(2))
    double precision             :: y1,y2,y3,y4
    double precision             :: t,u
    logical                      :: xmask(arrsize(1)),ymask(arrsize(2))
    integer                      :: i,j,jj,k,kk
    integer                      :: tmp(1)

! ERROR CHECKING ============================================================
    do i=1,arrsize(1)-1
        if (xvec(i) > xvec(i+1)) then
            print *,'ERROR: x-vector not monotonically increasing.'
            stop
        end if
    end do
    do j=1,arrsize(2)-1
        if (yvec(j) > yvec(j+1)) then
            print *,'ERROR: y-vector not monotonically increasing.'
            stop
        end if
    end do
    if (x0 < xvec(1) .OR. x0 > xvec(arrsize(1))) then
        print *,'ERROR: x-variable out of bounds of input x-vector'
        print *,x0, xvec(1), xvec(arrsize(1))
        stop
    end if
    if (y0 < yvec(1) .OR. y0 > yvec(arrsize(2))) then
        print *,'ERROR: y-variable out of bounds of input y-vector'
        stop
    end if
    tmpshape = shape(arr)
    if (tmpshape(1) /= arrsize(1) .OR. tmpshape(2) /= arrsize(2)) then
        print *,'ERROR: Array shape does not equal declared array size'
        stop
    end if
    tmp = shape(xvec)
    if (tmp(1) /= arrsize(1)) then
        print *,'ERROR: X-vector length does not match declared number of matrix columns'
        stop
    end if
    tmp = shape(yvec)
    if (tmp(1) /= arrsize(2)) then
        print *,'ERROR: Y-vector length does not match declared number of matrix rows'
        stop
    end if
!============================================================================

    xmask = (/ (.TRUE.,i=1,arrsize(1)) /)
    ymask = (/ (.TRUE.,i=1,arrsize(2)) /)

    tmpm = xvec-x0
    tmpn = yvec-y0

    tmp = minloc(abs(tmpm),xmask)
    j = tmp(1)
    xmask(j) = .FALSE.
    tmp = minloc(abs(tmpm),xmask)
    jj = tmp(1)
    j = min(j,jj)

    tmp = minloc(abs(tmpn),ymask)
    k = tmp(1)
    ymask(k) = .FALSE.
    tmp = minloc(abs(tmpn),ymask)
    kk = tmp(1)
    k = min(k,kk)

    y1 = arr(j,k)
    y2 = arr(j+1,k)
    y3 = arr(j+1,k+1)
    y4 = arr(j,k+1)

    t = (x0 - xvec(j))/(xvec(j+1)-xvec(j))
    u = (y0 - yvec(k))/(yvec(k+1)-yvec(k))

    ans = (1.D0-t)*(1.D0-u)*y1 + t*(1.D0-u)*y2 + t*u*y3 + (1.D0-t)*u*y4

end subroutine interp2dlin
