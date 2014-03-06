!============================================================================
! PROGRAM: solve_beta.f90 (Modified version of solveLE.f90)
!
! DESCRIPTION: Solves the Lane-Emden Equation using the Runge-Kutta-Fehlberg
!       method with a tolerance of 1.0E-8, then solves for the average value
!       beta (the ratio of gas pressure to total pressure).
!
!       See manual page for detailed notes. Uses Brent's method to solve a
!       nonlinear equation for beta, given the mass.
!
! WRITTEN:  John Lattanzio                                (Monash University)
!
! MODIFIED: Mikhail Klassen                             (McMaster University)
!       7 / July / 2010
!
! EMAIL: mikhail.klassen@gmail.com
!============================================================================
subroutine solve_beta(nn,hmax,tol,M,Bn,beta)
    use astro_constants
    implicit none
    double precision :: x,y,z,hh,hmax,hmin
    double precision :: n,ff,gg,ynew,nn
    external :: ff,gg
    double precision :: dh(6),xarg,yarg,zarg
    double precision :: k1,k2,k3,k4,k5,k6
    double precision :: l1,l2,l3,l4,l5,l6
    double precision :: tol,Ry,Rz,RR
    double precision :: x2,x3,x4,x5,x6,delta
    double precision :: radius,rhocrhobar,mx2dthetadx,zfin
    double precision :: M,Bn,beta,vol
    common n
    integer :: icount,iout

!============================================================================
    vol = 0.D0
    dh(1)=0.D0
    dh(2)=1.D0/4.D0
    dh(3)=3.D0/8.D0
    dh(4)=12.D0/13.D0
    dh(5)=1.D0
    dh(6)=1./2.D0
    ! Set the frequency of output lines or 99 for none
    iout=1

    n = nn

    ! Use line below to define custom tolerance
    !tol=1.d-8

    hmin=tol

    ! Start at x=tolerance
    x=tol

    ! Use series solution to begin near origin
    x2=x*x
    x3=x2*x
    x4=x2*x2
    x5=x4*x
    x6=x4*x2
    y=1.D0-x2/6.D0 + n*x4/120.D0-n*(8.D0*n-5.D0)*x6/(42.D0*360.D0)
    z=-x/3. +n*x3/30.-n*(8*n-5)*x5/(7*360)

    ! Start with hh=hmin=tol
    hh=x

    icount=0

! Main Loop =================================================================
100 continue
    ! We should never use the next line !
    if(y <= 0) then; go to 200; endif

    k1=hh*ff(x,y,z)
    l1=hh*gg(x,y,z)

    xarg=x+dh(2)*hh
    yarg=y+k1*0.25D0
    if(yarg <= 0) then; go to 200; endif

    ! Then non-integer n will cause problems.
    ! Reduce hh and restart this step.
    zarg=z+l1*0.25D0
    k2=hh*ff(xarg,yarg,zarg)
    l2=hh*gg(xarg,yarg,zarg)

    xarg=x+dh(3)*hh
    yarg=y+3.D0*(k1+3.D0*k2)/32.D0
    if(yarg <= 0) then; go to 200; endif
    zarg=z+3.D0*(l1+3.D0*l2)/32.D0
    k3=hh*ff(xarg,yarg,zarg)
    l3=hh*gg(xarg,yarg,zarg)

    xarg=x+dh(4)*hh
    yarg=y+(1932.D0*k1-7200.D0*k2+7296.D0*k3)/2197.D0
    if(yarg <= 0) then; go to 200; endif
    zarg=z+(1932.D0*l1-7200.D0*l2+7296.D0*l3)/2197.D0
    k4=hh*ff(xarg,yarg,zarg)
    l4=hh*gg(xarg,yarg,zarg)

    xarg=x+dh(5)*hh
    yarg=y+439.D0*k1/216.D0-8.D0*k2+3680.D0*k3/513.D0-845.D0*k4/4104.D0
    if(yarg <= 0) then; go to 200; endif
    zarg=z+439.D0*l1/216.D0-8.D0*l2+3680.D0*l3/513.D0-845.D0*l4/4104.D0
    k5=hh*ff(xarg,yarg,zarg)
    l5=hh*gg(xarg,yarg,zarg)

    xarg=x+0.5D0*hh
    yarg=y-8.D0*k1/27.D0+2.D0*k2-3544.D0*k3/2565.D0+1859.D0*k4/4104.D0-11.D0*k5/40.D0
    if(yarg <= 0) then; go to 200; endif
    zarg=z-8.D0*l1/27.D0+2.D0*l2-3544.D0*l3/2565.D0+1859.D0*l4/4104.D0-11.D0*l5/40.D0
    k6=hh*ff(xarg,yarg,zarg)
    l6=hh*gg(xarg,yarg,zarg)

    ! Error-controlling R values
    Ry=abs(k1/360.D0-128.D0*k3/4275.D0-2197.D0*k4/75240.D0+k5/50.D0+2.D0*k6/55.D0)/hh
    Rz=abs(l1/360.D0-128.D0*l3/4275.D0-2197.D0*l4/75240.D0+l5/50.D0+2.D0*l6/55.D0)/hh
    RR=min(Ry,Rz)
    if(RR /= 0) delta=0.84D0*(tol/RR)**0.25D0
    ! if R < tol then this step is accepted.....
        if(RR <= tol) then
            ynew=y+25.D0*k1/216.D0+1408.D0*k3/2565.D0+2197.D0*k4/4104.D0-0.2D0*k5
            ! ... provided y is still positive
            ! if this y < 0 then reduce hh and re-compute this step...
            if(ynew < 0) then
                if(hh > tol) then
                hh=0.5D0*hh
                go to 100
            else
            ! ... unless hh < tol = hmin. In this case we found the surface
            ! so PRINT and return.
            beta = ( 3.D0 / (x*x*x)) * vol
            ! Set beta = 1 if volume integration results in beta slightly above 1
            if(beta>1.D0) beta=1.D0
            !print *,'Volume Averaged Beta = ',beta
            return
            endif
        endif

        ! Step accepted, y > 0
        ! Update x and z also
        y=ynew
        x=x+hh
        z=z+25.D0*l1/216.D0+1408.D0*l3/2565.D0+2197.D0*l4/4104.D0-0.2D0*l5

        ! If required, PRINT out the x,y,z and hh for this step
        icount=icount+1

        ! Calculate beta at this radius
        if(iout /= 99) then
            if(mod(icount,iout) == 0) then
            call beta_func_solver(n,M,Bn,y,beta)
            vol = vol + beta * x*x*hh
            !print 10,x,y,z,hh,beta
            endif
        endif
        else
        ! Here RR > tol, so we must reduce hh.
        ! This step is NOT accepted.
        ! Uncomment following statement to print adjustments to hh
        ! if(iout /= 99) then; PRINT 20,Ry,Rz,tol,hh; endif
    endif

    ! Now adjust hh
    if(delta <= 0.5D0) then
        hh=0.5D0*hh
    else if (delta > 2) then
        hh=2.D0*hh
    else
        hh=delta*hh
    endif

    if(hh > hmax) hh = hmax
    if(y < 0) return
    if(hh < hmin) hh=hmin
    go to 100

200 continue
    ! Must reduce hh as y went negative on this step
    hh=0.5D0*hh
    ! But if this makes hh < tol, then we have finished !
    if(hh < tol) then
        beta = (3.D0 / (x*x*x)) * vol
        ! Set beta = 1 if volume integration results in beta slightly above 1
        if(beta>1.D0) beta=1.D0
        !print *,'Volume Averaged Beta = ',beta
        return
    endif
    ! Reduced hh....go back to calculate the next step.
    go to 100

! End of Loop ===============================================================

! format statements
10  format(2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,F12.8)
11  format(2x,1pe14.6,2x,0pf8.6,2x,f8.5,2x,f8.5)
20  format('   *Ry=',e11.4,' Rz=',e11.4,' tol=',e10.3,'   hh=',f6.4)
30  format(//,' Solution: radius            = ',f10.5,/,&
              '           rho(c)/rho(bar)   = ',f10.5,/,&
              '           -x^2 d(theta)/dx  = ',f10.5,//)

end subroutine solve_beta

! FUNCTIONS =================================================================

function ff(xx,yy,zz)
!   ff is z=ff(x,y,z)=dy/dz
    implicit none
    double precision :: xx,yy,zz
    double precision :: ff
    ff=zz
    return
end function ff

function gg(xx,yy,zz)
!   gg is z'=gg(x,y,z)=dz/dz
    implicit none
    double precision :: xx,yy,zz,n
    double precision :: gg
    common n

    gg=-2.D0*zz/xx -yy**n

    return
end function gg

function beta_func(beta,n,Bn,M,theta)
    ! Nonlinear function of beta depending only on stellar mass and polytropic
    ! index
    use astro_constants
    implicit none
    double precision :: beta_func
    double precision :: n,Bn,M,theta,beta

    beta_func = (a_rad / 3.D0) * (mu * H / kb)**(4.D0) * &
              (4.D0 * pi * Bn * Bn * Bn * G * G * G * &
              M * M * theta**(3.D0 - n)) * beta**(4.D0) + &
              beta - 1.D0

end function beta_func

! SOLVE FOR BETA USING BRENT'S METHOD =======================================

subroutine beta_func_solver(n,M,Bn,theta,beta)
    ! Calling program, which sets the bounds and the tolerance
    use astro_constants
    implicit none
    double precision, intent(in)  :: n,M,Bn,theta
    double precision, intent(out) :: beta
    double precision              :: a,b,tol,beta_func
    external beta_func
    a   = 0.D0
    b   = 1.D0
    tol = 1.D-8
    call zbrent_beta(a,b,tol,beta_func,n,M,Bn,theta,beta)

end subroutine beta_func_solver

subroutine zbrent_beta(x1,x2,tol,func,n,M,Bn,theta,beta)
    ! Subroutine to solve for beta, given the stellar mass M, dimensionless
    ! pressure theta, parameter Bn, and polytropic index n. Based on Brent's
    ! method. Returns beta.
    use astro_constants
    implicit none
    double precision, intent(in)   :: x1,x2,tol,n,M,Bn,theta
    double precision, intent(out)  :: beta
    double precision               :: func
    integer, parameter             :: itmax=100
    double precision, parameter    :: eps=epsilon(x1)
    integer                        :: iter
    double precision               :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

    a=x1                 ! Lower limit
    b=x2                 ! Upper limit

    fa=func(a,n,Bn,M,theta)  ! Lower function value limit
    fb=func(b,n,Bn,M,theta)  ! Upper function value limit


    if ((fa > 0.D0 .and. fb > 0.D0) .or. (fa < 0.D0 .and. fb < 0.D0)) then
    print *,'Error: Root must be bracketed for zbrent'
    stop
    end if
    c=b
    fc=fb
    do iter=1,itmax
    if ((fb > 0.D0 .and. fc > 0.D0) .or. (fb < 0.D0 .and. fc < 0.D0)) then
      c=a
      fc=fa
      d=b-a
      e=d
    end if
    if (abs(fc) < abs(fb)) then
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    end if
    tol1=2.D0*eps*abs(b)+0.5D0*tol
    xm=0.5D0*(c-b)
    if (abs(xm) <= tol1 .or. fb == 0.D0) then
      beta=b
      return
    end if
    if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
      s=fb/fa
      if (a == c) then
        p=2.D0*xm*s
        q=1.D0-s
      else
        q=fa/fc
        r=fb/fc
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.D0))
        q=(q-1.D0)*(r-1.D0)*(s-1.D0)
      end if
      if (p > 0.D0) q=-q
      p=abs(p)
      if (2.D0*p  <  min(3.D0*xm*q-abs(tol1*q),abs(e*q))) then
        e=d
        d=p/q
      else
        d=xm
        e=d
      end if
    else
      d=xm
      e=d
    end if
    a=b
    fa=fb
    b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
    fb=func(b,n,Bn,M,theta)
    end do
    print *,'zbrent: Exceeded maximum iterations'
    beta=b
end subroutine zbrent_beta
