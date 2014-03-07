subroutine update_stage(mass, md, Ld, Lms, &                    ! Input
                         n, r, stage)                           ! In/Out
  
  use protostellar_interface

  implicit none
  double precision                              :: mass, md, Ld, Lms, n, r
  integer                                       :: stage
  integer                                       :: models
  double precision                              :: Rms
  double precision                              :: Dn,Bn
  double precision                              :: rho_avg, rho_c
  double precision                              :: Pc, Tc
  double precision                              :: t1,t2,tol,corepress
  double precision, dimension(:),   allocatable :: Dns, Bns, ns
  double precision, dimension(:,:), allocatable :: nDns, nBns
  double precision, parameter                   :: frad = 0.33
  external                                      :: corepress

  !/*********************************************************************
  ! As a first test, we must solve for the central temperature of our
  ! polytropic protostar. To do this, we must first consult a lookup
  ! table with precomputed values for such things as the density gradient
  ! (Dn) and another relevant parameter (Bn), each solved for a
  ! particular polytropic index n. We then interpolate between these
  ! values to find appropriate values for the polytropic index of the
  ! protostar being tested.
  !
  ! A variety of tools for managing the lookup tables are contained in
  ! the module 'polytools', e.g. getNumModels, which retrieves a number
  ! 'models' that describes the granularity of our lookup table. The
  ! table itself may be recomputed to finer precision, in which case,
  ! subroutines such as this one need to be able to handle finer or
  ! coarser precomputed tables. After 'models' is retrieved, we can 
  ! allocate memory for our arrays.
  ! 
  ! 'get_model_table_cols' is another tool within the 'polytools'
  ! module that is used to extract particular columns from the lookup
  ! table. Here we extract indices n density gradient Dn and parameter
  ! Bn. We then construct arrays of size 2x'models' and interpolate
  ! along the columns to get values needed for our current protostar.
  !*********************************************************************/
  if ( stage .EQ. 1 ) then
      call getNumModels(models)
      allocate(ns(models))
      allocate(Dns(models))
      allocate(Bns(models))
      allocate(nDns(2,models))
      allocate(nBns(2,models))
      ! Get columns
      call get_model_table_cols(models,n=ns,Dn=Dns,Bn=Bns)

      ! Prepare arrays
      nDns(1,:) = ns
      nDns(2,:) = Dns
      nBns(1,:) = ns
      nBns(2,:) = Bns

      ! Interpolate to get values for our current protostar
      call interp2dlin(nDns,shape(nDns),[1.D0,2.D0],ns,2.D0,n,Dn)
      call interp2dlin(nBns,shape(nBns),[1.D0,2.D0],ns,2.D0,n,Bn)

      ! Calculate average density, central density, and central pressure
      rho_avg = 3. * mass / (4. * pi * r * r * r)
      rho_c   = Dn * rho_avg
      Pc     = (4. * pi)**(1./3.) * Bn * G * mass**(2./3.) * rho_c**(4./3.)

      ! Find the central temperature of the protostar using Brent's method
      t1  = 0.0                ! Lower limit (Kelvin)
      t2  = 1.D8               ! Upper limit (Kelvin)
      tol = 1.D-1              ! Tolerance
      call zbrent_coretemp(t1,t2,tol,corepress,Pc,rho_c,Tc)
  endif

  ! If the core temperature exceeds 1.5*10^6 K, upgrade stage 1->2
  ! and change the polytropic index n = 1.5
  if (( stage .EQ. 1 ) .AND. ( Tc .GE. 1.5E6 )) then
      stage = 2
      n     = 1.5
  endif

  ! If the deuterium exhausts, we switch to stage 3:
  ! 'core burning at variable Tc'
  if (( stage .EQ. 2 ) .AND. ( md .LE. 0.0 )) then
      stage = 3
  endif

  ! Check if radiative zone has formed by comparing deuterium
  ! luminosity Ld to the luminosity of an equal-mass main-sequence star
  ! if Ld/Lms > frad = 0.33, switch to n = 3 and increase radius by
  ! factor of 2.1. This represents the swelling of the star at the
  ! formation of the radiative barrier.
  if (( stage .EQ. 3 ) .AND. ( Ld/Lms .GT. frad )) then
      stage = 4
      r = 2.1 * r
      n = 3.0
  endif

  ! Finally, check if the radius has approached that of an equal-mass
  ! main-sequence star
  call calculate_R_main_sequence(mass/Msun,Rms)
  Rms = Rms * Rsun                 ! Scaling to cm from solar units
  
  if (( stage .EQ. 4 ) .AND. ( r .LE. Rms )) then
      stage = 5
      r = Rms
  endif

  ! At stage 5, the protostar is now approximated as a star on the
  ! main sequence
  if ( stage .EQ. 5) then
      r = Rms
  endif

  if ( r .LE. Rms ) then
      r = Rms
  endif

  return

end subroutine update_stage

function corepress(Pc, rho_c, Tc)
    implicit none
    double precision            :: corepress
    double precision            :: Pc, rho_c, Tc
    double precision, parameter :: kb=1.380658E-16 ! erg/K
    double precision, parameter :: mu=0.613
    double precision, parameter :: mH=1.673534E-24 ! g
    double precision, parameter :: a = 7.56591E-15 ! erg/cm^3/K^4
  
    corepress = (rho_c * kb * Tc)/(mu * mH) + (1./3.)*a*Tc**4. - Pc

end function corepress

subroutine zbrent_coretemp(T1,T2,tol,func,Pc,rho_c,rt)
  implicit none
  double precision, intent(in)  :: T1,T2,tol,Pc,rho_c
  double precision, intent(out) :: rt
  double precision              :: func
  integer, parameter            :: itmax=100
  double precision, parameter   :: eps=epsilon(T1)
  integer                       :: iter
  double precision              :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  
  a=T1
  b=T2
  fa=func(Pc,rho_c,a)
  fb=func(Pc,rho_c,b)
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
    print *,'root must be bracketed for zbrent'
    stop
  end if
  c=b
  fc=fb
  do iter=1,itmax
    if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
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
    tol1=2.0*eps*abs(b)+0.5*tol
    xm=0.5*(c-b)
    if (abs(xm) <= tol1 .or. fb == 0.0) then
      rt=b
      return
    end if
    if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
      s=fb/fa
      if (a == c) then
        p=2.0*xm*s
        q=1.0-s
      else
        q=fa/fc
        r=fb/fc
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
        q=(q-1.0)*(r-1.0)*(s-1.0)
      end if
      if (p > 0.0) q=-q
      p=abs(p)
      if (2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
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
    fb=func(Pc,rho_c,b)
  end do
  print *,'zbrent: exceeded maximum iterations'
  rt=b
end subroutine zbrent_coretemp
