module analytic_ionization_solutions_mod
use myf90_mod
use hui_gnedin_atomic_rates_mod
use ion_table_class, only: mini_spec_type
use cen_atomic_rates_mod, only: Verner_HI_photo_cs
use ion_table_class, only: gammaHI_from_mini_spec_shield_plane
use config_mod, only: CV
implicit none
private


real(r8b), parameter :: zero = 0.0d0
real(r8b), parameter :: half = 0.5d0
real(r8b), parameter :: one = 1.0d0
real(r8b), parameter :: two = 2.0d0
real(r8b), parameter :: four = 4.0d0

real(r8b), parameter :: TAUHI_EFF_MAX = 20.0d0
real(r8b), parameter :: GHI_RATIO_MIN = exp(-TAUHI_EFF_MAX)

public :: analytic_xHIeq
public :: bracket_xHI_plane_parallel
public :: bracket_GHI_plane_parallel
public :: numerical_integration_of_slab

contains


!> Ionization State Diff EQ, Analytic expressions
!================================================================
!================================================================


! smoothing of transition from caseA to caseB
!
! we use tauHI_eff = -ln(GHI/GHI_thin) to determine 
! if we should use caseA, caseB, or something in between
! first we scale tauHI_th by a user defined value
! tf = tauHI_eff / CV%ShldTauThresh
! then we set the range of tf between which we transition
! between caseA and caseB using t_lo and t_hi and do a 
! linear interpolation to get the recomb. rate
! RC = RCa - frac * RC1 
! where frac goes from 0 to 1 as tf goes from t_lo to t_hi
!----------------------------------------------------
function tauHI_eff_to_RC( tauHI_eff, T ) result( RC )
  real(r8b) :: tauHI_eff
  real(r8b) :: T
  real(r8b) :: RC

  real(r8b) :: RCa, RCb, RC1
  real(r8b) :: frac, tf

  real(r8b), parameter :: t_lo = zero
  real(r8b), parameter :: t_hi = 3*one
  real(r8b), parameter :: dtau = t_hi - t_lo

  if ( tauHI_eff >= TAUHI_EFF_MAX ) then
     RCb = Hui_HII_recombB( T )
     RC = RCb
     return
  endif

  tf = tauHI_eff - CV%ShldTauThresh

  if ( tf > t_hi ) then
     RCb = Hui_HII_recombB( T )
     RC = RCb
  else if ( tf < t_lo ) then
     RCa = Hui_HII_recombA( T )
     RC = RCa
  else
     RCa = Hui_HII_recombA( T )
     RCb = Hui_HII_recombB( T )
     RC1 = RCa - RCb
     frac = ( tf - t_lo ) / dtau
     RC = RCa - frac * RC1

     if ( RC /= RC ) then
        write(*,*) 'RC is NaN'
        write(*,*) 'tauHI_eff: ', tauHI_eff
        write(*,*) 'tf: ', tf
        stop
     endif

  end if

end function tauHI_eff_to_RC


!> time derivative, neutral fraction
!! returns dxHI/dt
!-------------------------------------------------------------
function dxHIdt( T, nH, GHI, y, xHI, tauHI_eff) result(dxdt)
  real(r8b) :: T, nH, GHI, y, xHI, tauHI_eff
  real(r8b) dxdt
  real(r8b) :: RC,CI,xHII,ne
  
  RC = tauHI_eff_to_RC( tauHI_eff, T )
  CI = Hui_HI_col_ion( T )

  xHII = one - xHI
  ne = (xHII + y) * nH

  dxdt = -(GHI + CI * ne) * xHI + RC * ne * xHII 
  
end function dxHIdt


!> time derivative, ionized fraction
!! returns dxHII/dt
!-------------------------------------------------------------
function dxHIIdt( T, nH, GHI, y, xHII, tauHI_eff ) result(dxdt)
  real(r8b) :: T, nH, GHI, y, xHII, tauHI_eff
  real(r8b) dxdt
  real(r8b) :: RC,CI,xHI,ne

  RC = tauHI_eff_to_RC( tauHI_eff, T )
  CI = Hui_HI_col_ion( T )

  xHI = one - xHII
  ne = (xHII + y) * nH

  dxdt = (GHI + CI * ne) * xHI - RC * ne * xHII
  
end function dxHIIdt



!> implements analytic solution for xH's 
!-----------------------------------------------------------------------
function analytic_xHIeq( T, GHI, nH, y, tauHI_eff ) result(xHI)
  real(r8b), intent(in) :: T, GHI, nH, y
  real(r8b), intent(in) :: tauHI_eff
  real(r8b) :: xHI
  real(r8b) :: CI, RC, R, Q, P, d, d2, qq

  real(r8b), parameter :: xHI_min = 1.0d-12
  real(r8b), parameter :: T_MIN = 2.2d3
  real(r8b), parameter :: GHI_MIN = 1.0d-35

  ! if cool and shielded set to maximum neutral
  !-------------------------------------------------
  if ( T < T_MIN .and. GHI < GHI_MIN ) then
     xHI = one - xHI_min
     return
  endif

  ! get recombination rate and collisional ionization rate
  !--------------------------------------------------------
  RC = tauHI_eff_to_RC( tauHI_eff, T )
  CI = Hui_HI_col_ion( T )


  ! calculate coefficients
  !----------------------------------------------------
  R = ( CI + RC ) * nH
  Q = -( GHI + (CI + 2*RC) * nH + (CI + RC) * nH * y )
  P = RC * (1+y) * nH 

  d2 = Q * Q - four * R * P

  if (d2 <= zero) then
     d = zero
  else
     d = sqrt(d2)
  endif


  ! if d is very small that means neutral particle
  !----------------------------------------------------
  if ( d <= zero ) then
     write(*,*) 
     write(*,*) 'T: ', T
     write(*,*) 'GHI: ', GHI
     write(*,*) 'nH: ', nH
     write(*,*) 'y: ', y
     write(*,*) 
     write(*,*) 'RC: ', RC
     write(*,*) 'CI: ', CI
     write(*,*) 'R: ', R
     write(*,*) 'Q: ', Q
     write(*,*) 'P: ', P
     write(*,*) 'sign(Q): ', sign(1.0d0,Q)
     write(*,*) 'Q^2: ', Q*Q
     write(*,*) '4PR: ', 4*P*R
     write(*,*) 'd2: ', d2
     write(*,*) 'd: ', d
     write(*,*) 'd <= 0'
     write(*,*)
     xHI = one - xHI_min
     return
  endif


  ! do calculation if its floating point safe
  !----------------------------------------------------
!  qq = -0.5d0 * ( Q + sign(1.0d0,Q) * d )
  qq = -half * ( Q - d )  ! Q is always negative
  xHI = P / qq



  ! still better dummy check the results for NaN
  !----------------------------------------------------  
  if (.not. xHI >= one .and. .not. xHI < one)  then
     write(*,*) 
     write(*,*) 'T: ', T
     write(*,*) 'GHI: ', GHI
     write(*,*) 'nH: ', nH
     write(*,*) 'y: ', y
     write(*,*) 
     write(*,*) 'RC: ', RC
     write(*,*) 'CI: ', CI
     write(*,*) 'R: ', R
     write(*,*) 'Q: ', Q
     write(*,*) 'P: ', P
     write(*,*) 'sign(Q): ', sign(one,Q)
     write(*,*) 'Q^2: ', Q*Q
     write(*,*) '4PR: ', 4*P*R
     write(*,*) 'd2: ', d2
     write(*,*) 'd: ', d
     write(*,*) 'xHI is NaN'
     write(*,*) 'xHI = ', xHI
     write(*,*)
     xHI = one - xHI_min
     return
  endif


  ! still better dummy check the results for bounds
  !----------------------------------------------------  
  if (xHI >= one .or. xHI <= zero)  then
     write(*,*) 
     write(*,*) 'T: ', T
     write(*,*) 'GHI: ', GHI
     write(*,*) 'nH: ', nH
     write(*,*) 'y: ', y
     write(*,*) 
     write(*,*) 'RC: ', RC
     write(*,*) 'CI: ', CI
     write(*,*) 'R: ', R
     write(*,*) 'Q: ', Q
     write(*,*) 'P: ', P
     write(*,*) 'sign(Q): ', sign(one,Q)
     write(*,*) 'Q^2: ', Q*Q
     write(*,*) '4PR: ', 4*P*R
     write(*,*) 'd2: ', d2
     write(*,*) 'd: ', d

     if (xHI >= one) then
        write(*,*) 'xHI >= 1 in analytic'
        xHI = one - xHI_min
     else
        write(*,*) 'xHI <= 0 in analytic'
        xHI = xHI_min
     endif

     write(*,*) 'xHI = ', xHI
     write(*,*)
  endif



end function analytic_xHIeq








!> Plane parallel monochromatic radiation incident on a slab of 
!! of uniform density and temperature.  These analytic solutions 
!! only consider pure hydrogen chemistry.
!!
!! x0 is the solution at the surface of the slab
!! AA is collisional ionization equilibrium (i.e. GHI=0)
!==================================================================
!==================================================================


! zero of plane parallel solution, neutral fraction
! (y=0, monochromatic).  In this case, xHI is bound 
! between x0 < xHI < AA = RC/(CI+RC) = xHI^CIeq,
!--------------------------------------------------------------
function zero_plane_parallel_xHI( T, nH, GHI0, tauH, xHI, tauHI_eff ) result(eval)
  real(r8b) :: T, nH, GHI0, tauH, xHI, tauHI_eff
  real(r8b) eval
  real(r8b) :: RC, CI, xHI0, AA
  real(r8b) :: top, bot, topA, botA
  real(r8b) :: yy 

  logical :: err

  yy = zero

  ! xHI at surface of slab
  xHI0 = analytic_xHIeq( T, GHI0, nH, yy, tauHI_eff )

  RC = tauHI_eff_to_RC( tauHI_eff, T )
  CI = Hui_HI_col_ion( T )

  ! collisional ionization equilibrium
  AA = RC / (RC+CI)

  ! calculate zero
  top = xHI * (one-xHI0)
  bot = xHI0 * (one-xHI)

  topA = xHI * (AA-xHI0)
  botA = xHI0 * (AA-xHI)

  eval = (one / xHI0 - one / xHI)    &
       + log( top /bot  )            &
       + log( topA/botA ) / AA       &
       - tauH


  !write(*,*) 'xHI = ', xHI
  !write(*,*) 'xHI0 = ', xHI0
  !write(*,*) 'top/bot = ', top/bot
  !write(*,*) 'topA/botA = ', topA/botA
  !write(*,*) 'AA = ', AA
  !write(*,*) 'eval = ', eval

end function zero_plane_parallel_xHI






! given tauH (depth of slab) and GHI0 (photoionizarion rate at 
! surface of slab) find the plane parallel xHI at tauH
!-------------------------------------------------------------
function bracket_xHI_plane_parallel( T, nH, GHI0, tauH, GHI_thin, GHI ) result( xHI )
  real(r8b), intent(in) :: T, nH, GHI0, tauH, GHI_thin
  real(r8b) :: GHI
  real(r8b) :: xHI

  real(r8b) :: tauHI_eff
  real(r8b) :: CI, xHI0, AA
  real(r8b) :: RC, RCa, RCb
  real(r8b) :: x_lo, x_mid, x_hi
  real(r8b) :: f_lo, f_mid, f_hi, f1
  real(r8b) :: tauH_max, tauH_min
  real(r8b) :: yy

  real(r8b) :: tauHI_high
  real(r8b) :: xHI_max

  real(r8b), parameter :: tol = 1.0d-6
  !real(r8b), parameter :: x_max = 1.0d0 - tol

  ! calculate rates
  !-------------------------------------------------------
  RCa = Hui_HII_recombA( T )
  RCb = Hui_HII_recombB( T )
  CI = Hui_HI_col_ion( T )

  ! calculate minimum possible xHI (at surface of slab)
  ! use caseB here -> use very high optical depth
  !--------------------------------------------------------
  yy = zero
  tauHI_high = TAUHI_EFF_MAX
  xHI0 = analytic_xHIeq( T, GHI0, nH, yy, tauHI_high  )

  ! calculate maximum possible xHI (GHI=0)
  ! use caseA > caseB here
  !--------------------------------------------------------
  AA = RCa / (RCa + CI)
  xHI_max = min( AA, one ) - tol*1.0d-1

  ! if the difference between them is less than the tol
  ! then return their mean as the solution
  !------------------------------------------------------
  if ( abs(AA - xHI0) <= tol ) then
     xHI = half * (AA + xHI0)
     return
  endif


  ! if the supplied tauH is greater than the max 
  ! tauH allowed, then return the max neutral fraction
  !------------------------------------------------------
  tauH_max = &
       (one/xHI0 - one/xHI_max) + &
       log( xHI_max * (one-xHI0) / xHI0 / (one-xHI_max) ) +    &
       log( xHI_max * (AA-xHI0)  / xHI0 / (AA-xHI_max) ) / AA

  if (tauH_max /= tauH_max) then
     write(*,*) 'tauH_max is NaN - green'
     stop
  endif

  if (tauH .ge. tauH_max) then
     xHI = xHI_max
     return
  endif

  ! set bounds for bracketing
  !-------------------------------------------------------
  x_lo = xHI0  
  x_hi = xHI_max ! min( x_max, AA - tol )

  ! make sure x_hi > x_lo
  !-------------------------------------------------------
  if (xHI0 >= x_hi) then
     xHI = x_hi
     return
  endif

  ! we dont know the xHI or the GHI at tauH yet so we take
  ! a first guess at tauHI_eff,RC and update as we go
  !-------------------------------------------------------
  tauHI_eff = -log( GHI0/GHI_thin )
  if ( tauHI_eff /= tauHI_eff ) then
     write(*,*) 'a tauHI_eff is NaN'
     write(*,*) 'GHI0: ', GHI0
     write(*,*) 'GHI_thin: ', GHI_thin
     stop
  endif


  RC = tauHI_eff_to_RC( tauHI_eff, T )

  f_lo = zero_plane_parallel_xHI( T, nH, GHI0, tauH, x_lo, tauHI_high )
  f_hi = zero_plane_parallel_xHI( T, nH, GHI0, tauH, x_hi, tauHI_high )

  !f_lo = zero_plane_parallel_xHI( T, nH, GHI0, tauH, x_lo, tauHI_eff=zero )
  !f_hi = zero_plane_parallel_xHI( T, nH, GHI0, tauH, x_hi, tauHI_eff=zero )
 
  if (f_lo * f_hi > zero) then
     write(*,*) 'no zero bracketed in bracket_xHI_plane_parallel'
     write(*,*) 'T:     ', T
     write(*,*) 'nH:    ', nH
     write(*,*) 'GHI0:   ', GHI0
     write(*,*) 'tauH:  ', tauH
     write(*,*) 'xHI0:    ', xHI0
     write(*,*) 'AA: ', AA
     write(*,*) 'xHI_max: ', xHI_max
     write(*,*) 'x_lo, x_hi: ', x_lo, x_hi
     write(*,*) 'f_lo, f_hi: ', f_lo, f_hi
     write(*,*) 'tauH_max: ', tauH_max
     write(*,*) 'RC: ', RC
     write(*,*) 'RCa: ', RCa
     write(*,*) 'RCb: ', RCb
     write(*,*)
     stop
  endif


  if (x_hi - x_lo <= tol) then
     x_mid = half * (x_hi + x_lo)
     xHI = x_mid
     return
  endif
 
  do while (x_hi - x_lo > tol)
     
     x_mid = half * (x_hi + x_lo)

     GHI = ( RC * nH * (one-x_mid)**2 - CI * nH * (one-x_mid) * x_mid ) / x_mid

     if (GHI <= zero) then
        GHI = 1.0d-35
        xHI = x_mid
        return 
     endif

     if ( GHI/GHI_thin <= GHI_RATIO_MIN ) then
        RC = RCb
     else
        tauHI_eff = -log( GHI/GHI_thin )
        RC = tauHI_eff_to_RC( tauHI_eff, T )
     endif

     if ( tauHI_eff /= tauHI_eff ) then
        write(*,*) 'b tauHI_eff is NaN'
        write(*,*) 'GHI0: ', GHI0
        write(*,*) 'tauH: ', tauH
        write(*,*) 'xmid: ', x_mid
        write(*,*) 'GHI: ', GHI
        write(*,*) 'GHI_thin: ', GHI_thin
        write(*,*) 'ratio: ', GHI/GHI_thin 
        write(*,*) 'log: ', -log( GHI/GHI_thin )
        stop
     endif

     f_mid = zero_plane_parallel_xHI( T, nH, GHI0, tauH, x_mid, tauHI_eff )
     
     if (f_lo * f_mid > zero) then 
        x_lo = x_mid
        f_lo = f_mid
     else 
        x_hi = x_mid
        f_hi = f_mid
     endif
     
  end do

  x_mid = half * (x_hi + x_lo)
  xHI = x_mid

  return


end function bracket_xHI_plane_parallel



! given tau_0 (tau up to the particle, but not including the 
! contribution of its neighbors) and gamma_HI find the plane parallel xHI
!-------------------------------------------------------------
function bracket_GHI_plane_parallel( T, nH, GHI0, tauH, GHI_thin ) result( GHI )
  real(r8b), intent(in) :: T, nH, GHI0, tauH, GHI_thin
  real(r8b) :: GHI

  real(r8b) :: tauHI_eff
  real(r8b) :: xHI
  real(r8b) :: RC, CI

  real(r8b), parameter :: tol = 1.0d-6
  real(r8b), parameter :: x_max = 1.0d0 - tol


  ! get xHI and calculate GHI
  !------------------------------------------------------
  xHI = bracket_xHI_plane_parallel( T, nH, GHI0, tauH, GHI_thin, GHI )

  ! get recombination and collisional ionization rates
  !------------------------------------------------------
  tauHI_eff = -log( GHI/GHI_thin )
  if ( tauHI_eff /= tauHI_eff ) then
     write(*,*) 'c tauHI_eff is NaN'
     write(*,*) 'GHI0: ', GHI0
     write(*,*) 'tauH: ', tauH
     write(*,*) 'xHI: ', xHI
     write(*,*) 'GHI: ', GHI
     write(*,*) 'GHI_thin: ', GHI_thin
     write(*,*) 'ratio: ', GHI/GHI_thin 
     write(*,*) 'log: ', -log( GHI/GHI_thin )
     stop
  endif

  RC = tauHI_eff_to_RC( tauHI_eff, T )
  CI = Hui_HI_col_ion( T )

  GHI = return_GHI_for_xHI_plane_parallel( nH, RC, CI, xHI )

  return


end function bracket_GHI_plane_parallel


! use analytic relation to get GHI at given xHI for slab
!-------------------------------------------------------------
function return_GHI_for_xHI_plane_parallel( nH, RC, CI, xHI ) result( GHI )
  real(r8b) :: nH, RC, CI, xHI
  real(r8b) :: GHI

  GHI = ( RC * nH * (one-xHI)**2 - CI * nH * (one-xHI) * xHI ) / xHI

end function return_GHI_for_XHI_plane_parallel

! expensive numerical integration of slab to get 
! GHI at a given tauH into slab
!---------------------------------------------------
function numerical_integration_of_slab( mini_spec, tauH_slab, tauHI_pre, Tslab, nH ) &
     result(GHI_out)
  type(mini_spec_type), intent(in) :: mini_spec
  real(r8b), intent(in) :: tauH_slab   ! optical depth of slab in total H
  real(r8b), intent(in) :: tauHI_pre   ! optical depth before slab in neutral H
  real(r8b), intent(in) :: Tslab       ! temperature of slab
  real(r8b), intent(in) :: nH          ! number density of Hydrogen in slab
  real(r8b) :: GHI_out

  real(r8b), parameter :: tauH_per_layer = 1.0d1
  real(r8b), parameter :: y = zero
  logical, save :: first = .true.
  real(r8b), save :: sigmaHI_th


  integer(i8b) :: i   ! 
  integer(i8b) :: N   ! number of layers
  real(r8b) :: Lslab
  real(r8b) :: dz
  real(r8b) :: tauHI
  real(r8b) :: dtauH
  real(r8b) :: dtauHI
  real(r8b) :: tauH, tauH_prior
  real(r8b) :: xHI,  xHI_prior, xHI_tmp
  real(r8b) :: GHI,  GHI_prior
  real(r8b) :: frac, tauH_lo, tauH_hi
  real(r8b) :: tauHI_eff

  if (first) then
     sigmaHI_th = Verner_HI_photo_cs(one)    
  endif
  first = .false.


  ! set number of layers so that each must be optically thin
  !----------------------------------------------------------
  N     = tauH_slab / tauH_per_layer + 1
  Lslab = tauH_slab / ( nH * sigmaHI_th )   ! [cm]
  dz    = Lslab / N                         ! [cm]
  dtauH = tauH_per_layer


  ! initialize tracking variables
  !----------------------------------------------------------
  tauH   = zero       ! to surface of layer
  tauHI  = tauHI_pre   ! to surface of layer

  dtauHI  = zero
  xHI     = one
  xHI_tmp = half

  do while ( abs(xHI/xHI_tmp-one) > 1.0d-4 )
     xHI_tmp   = xHI
     GHI       = gammaHI_from_mini_spec_shield_plane( mini_spec, tauHI_th=tauHI+dtauHI*half ) 
     tauHI_eff = -log(GHI/mini_spec%gammaHI_int)
     xHI       = analytic_xHIeq( Tslab, GHI, nH, y, tauHI_eff )
     dtauHI    = dz * nH * xHI * sigmaHI_th
  end do


  do i = 1, N

     ! remember previous values if on last step
     !-------------------------------------------
     if (i==N) then
        tauH_prior = tauH
        xHI_prior  = xHI
        GHI_prior  = GHI
     endif

     ! increment values
     !-------------------------------------------
     tauH   = tauH + dtauH
     tauHI  = tauHI + dtauHI

     xHI_tmp = xHI+1.0d3
     do while ( abs(xHI/xHI_tmp-1.0d0) > 1.0d-4 )
        xHI_tmp   = xHI
        GHI       = gammaHI_from_mini_spec_shield_plane( mini_spec, tauHI_th=tauHI+dtauHI*0.5d0 ) 
        tauHI_eff = -log(GHI/mini_spec%gammaHI_int)
        xHI       = analytic_xHIeq( Tslab, GHI, nH, y, tauHI_eff )
        dtauHI    = dz * nH * xHI * sigmaHI_th
     end do

  end do


  if ( tauH_slab < tauH_prior .or. tauH_slab > tauH ) then
     write(*,*) 'problem tracking tauH **** '
     stop
  endif

      
  tauH_lo = tauH_prior
  tauH_hi = tauH
     
  frac = (tauH_slab - tauH_lo) / (tauH_hi - tauH_lo)        
     
  if (GHI > 1.0d-99) then
     GHI_out = log10(GHI_prior) + frac * ( log10(GHI) - log10(GHI_prior) )
     GHI_out = 10.0d0**GHI_out
  else
     GHI_out = 0.0d0
  endif
     

!  write(*,*) 'GHI_out: ', GHI_out


  


end function numerical_integration_of_slab






end module analytic_ionization_solutions_mod
