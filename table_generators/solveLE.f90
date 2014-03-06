subroutine solveLE(nn,iout,hmax,tol,radius,rhocrhobar,mx2dthetadx,zfin)
!============================================================================
! SUBROUTINE: SolveLE.F90
!
! DESCRIPTIONS: Solves the Lane-Emden Equation using the Runge-Kutta-Fehlberg
!       method with a tolerance of 1.0E-8 or user specified.
!
! WRITTEN:  John Lattanzio                                (Monash University)
!
! MODIFIED: Mikhail Klassen                             (McMaster University)
!       2 / July / 2010
! EMAIL: mikhail.klassen@gmail.com
!============================================================================
      implicit none
! VARIABLE DECLARATIONS =====================================================
      double precision x,y,z,h,hmax,hmin
      double precision n,f,g,ynew,nn
      external f,g
      double precision dh(6),xarg,yarg,zarg
      double precision k1,k2,k3,k4,k5,k6
      double precision l1,l2,l3,l4,l5,l6
      double precision tol,Ry,Rz,R
      double precision x2,x3,x4,x5,x6,delta
      double precision radius,rhocrhobar,mx2dthetadx,zfin
      common n
      integer icount,iout
!============================================================================
      dh(1)=0.D0
      dh(2)=1.D0/4.D0
      dh(3)=3.D0/8.D0
      dh(4)=12.D0/13.D0
      dh(5)=1.D0
      dh(6)=1./2.D0

! USER INPUT ================================================================
!      WRITE(6,*) ' Enter polytropic index n'
!      READ(5,*) n
       n = nn
!      WRITE(6,*) '      '
!      WRITE(6,*) ' Output every iout steps: Enter iout (99 for none)'
!      READ(5,*) iout
!      tol=1.d-8
!!     Place restrictions on hmax and hmin
!!     These are also used to find the edge accurately
!      WRITE(6,*) '      '
!      WRITE(6,*) ' Enter maximum step-length h (larger for larger n)'
!      READ(5,*) hmax
!      WRITE(6,*) ' Enter tolerance = maximum error'
!      READ(5,*) tol
      hmin=tol
!============================================================================

!     Start at x=tolerance
      x=tol

!     Use series solution to begin near origin
      x2=x*x
      x3=x2*x
      x4=x2*x2
      x5=x4*x
      x6=x4*x2
      y=1.D0-x2/6.D0 + n*x4/120.D0-n*(8.D0*n-5.D0)*x6/(42.D0*360.D0)
      z=-x/3. +n*x3/30.-n*(8*n-5)*x5/(7*360)

!     Start with h=hmin=tol
      h=x

!     Print first lines
      if(iout /= 99) PRINT 11
 11   format('          x            y         z        h')
      if(iout /= 99) PRINT 10,x,y,z,h

      icount=0

! Main Loop =================================================================
 100  continue
!     We should never use the next line !
      if(y <= 0) then; go to 200; endif

      k1=h*f(x,y,z)
      l1=h*g(x,y,z)

      xarg=x+dh(2)*h
      yarg=y+k1*0.25D0
      if(yarg <= 0) then; go to 200; endif
!     Then non-integer n will cause problems.
!     Reduce h and restart this step.
      zarg=z+l1*0.25D0
      k2=h*f(xarg,yarg,zarg)
      l2=h*g(xarg,yarg,zarg)

      xarg=x+dh(3)*h
      yarg=y+3.D0*(k1+3.D0*k2)/32.D0
      if(yarg <= 0) then; go to 200; endif
      zarg=z+3.D0*(l1+3.D0*l2)/32.D0
      k3=h*f(xarg,yarg,zarg)
      l3=h*g(xarg,yarg,zarg)

      xarg=x+dh(4)*h
      yarg=y+(1932.D0*k1-7200.D0*k2+7296.D0*k3)/2197.D0
      if(yarg <= 0) then; go to 200; endif
      zarg=z+(1932.D0*l1-7200.D0*l2+7296.D0*l3)/2197.D0
      k4=h*f(xarg,yarg,zarg)
      l4=h*g(xarg,yarg,zarg)

      xarg=x+dh(5)*h
      yarg=y+439.D0*k1/216.D0-8.D0*k2+3680.D0*k3/513.D0-845.D0*k4/4104.D0
      if(yarg <= 0) then; go to 200; endif
      zarg=z+439.D0*l1/216.D0-8.D0*l2+3680.D0*l3/513.D0-845.D0*l4/4104.D0
      k5=h*f(xarg,yarg,zarg)
      l5=h*g(xarg,yarg,zarg)

      xarg=x+0.5D0*h
      yarg=y-8.D0*k1/27.D0+2.D0*k2-3544.D0*k3/2565.D0+1859.D0*k4/4104.D0-11.D0*k5/40.D0
      if(yarg <= 0) then; go to 200; endif
      zarg=z-8.D0*l1/27.D0+2.D0*l2-3544.D0*l3/2565.D0+1859.D0*l4/4104.D0-11.D0*l5/40.D0
      k6=h*f(xarg,yarg,zarg)
      l6=h*g(xarg,yarg,zarg)


!     Error-controlling R values
      Ry=abs(k1/360.D0-128.D0*k3/4275.D0-2197.D0*k4/75240.D0+k5/50.D0+2.D0*k6/55.D0)/h
      Rz=abs(l1/360.D0-128.D0*l3/4275.D0-2197.D0*l4/75240.D0+l5/50.D0+2.D0*l6/55.D0)/h
      R=min(Ry,Rz)
      if(R /= 0) delta=0.84D0*(tol/R)**0.25D0
!       if R < tol then this step is accepted.....
        if(R <= tol) then
           ynew=y+25.D0*k1/216.D0+1408.D0*k3/2565.D0+2197.D0*k4/4104.D0-0.2D0*k5
!          ... provided y is still positive
!          if this y < 0 then reduce h and re-compute this step...
           if(ynew < 0) then
              if(h > tol) then
                 h=0.5D0*h
                 go to 100
              else
!               ... unless h < tol = hmin. In this case we found the surface
!               so PRINT and return.
300             radius = x
      	        rhocrhobar = -x/(3.D0*z)
      	        mx2dthetadx = -x*x*z
      	        zfin = z
      	        !PRINT 30,x,-x/(3*z),-x*x*z
                return
              endif
           endif

!          Step accepted, y > 0
!          Update x and z also
           y=ynew
           x=x+h
           z=z+25.D0*l1/216.D0+1408.D0*l3/2565.D0+2197.D0*l4/4104.D0-0.2D0*l5

!          If required, PRINT out the x,y,z and h for this step
        icount=icount+1
           if(iout /= 99) then
              if(mod(icount,iout) == 0) PRINT 10,x,y,z,h
           endif
        else
!       Here R > tol, so we must reduce h.
!       This step is NOT accepted.
!       Uncomment following statement to print adjustments to h
!       if(iout /= 99) then; PRINT 20,Ry,Rz,tol,h; endif
      endif

!     Now adjust h
      if(delta <= 0.5D0) then
         h=0.5D0*h
      else if (delta > 2) then
         h=2.D0*h
      else
         h=delta*h
      endif

      if(h > hmax) h = hmax
      if(y < 0) return
      if(h < hmin) h=hmin
      go to 100

200   continue
!     Must reduce h as y went negative on this step
      h=0.5D0*h
!     But if this makes h < tol, then we have finished !
      if(h < tol) then
      	 radius = x
      	 rhocrhobar = -x/(3.D0*z)
      	 mx2dthetadx = -x*x*z
      	 zfin = z
         !PRINT 30,x,-x/(3*z),-x*x*z
             return
      endif
!     Reduced h....go back to calculate the next step.
      go to 100
! End of Loop ===============================================================


! format statements
10    format(2x,1pe14.6,2x,0pf8.6,2x,f8.5,2x,f8.5)
20    format('   *Ry=',e11.4,' Rz=',e11.4,' tol=',e10.3,'   h=',f6.4)
30    format(//,' Solution: radius            = ',f10.5,/,&
               '           rho(c)/rho(bar)   = ',f10.5,/,&
               '           -x^2 d(theta)/dx  = ',f10.5,//)

end subroutine solveLE

! FUNCTIONS =================================================================
!     f is z=f(x,y,z)=dy/dz
double precision function f(xx,yy,zz)
      double precision xx,yy,zz
      f=zz
      return
end function f

!     g is z'=g(x,y,z)=dz/dz
double precision function g(xx,yy,zz)
      common n
      double precision xx,yy,zz,n
      g=-2.D0*zz/xx -yy**n
      return
end function g
