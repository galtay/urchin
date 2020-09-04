!> \file ionpar.F90

!> \brief the module that handles ionization particles
!<

module ionpar_mod
use myf90_mod
use atomic_rates_mod, only: atomic_rates_type
use particle_system_mod, only: particle_type
use raylist_mod, only: raylist_type
use raylist_mod, only: intersection_type

use particle_returns_mod, only: return_H_mf
use particle_returns_mod, only: return_H2_mf
use particle_returns_mod, only: return_eos
use particle_returns_mod, only: return_fhot

use cen_atomic_rates_mod, only: Verner_HI_photo_cs
use cen_atomic_rates_mod, only: Osterbrok_HeI_photo_cs
use cen_atomic_rates_mod, only: Osterbrok_HeII_photo_cs
use cen_atomic_rates_mod, only: Haiman_Bremss_cool
use cen_atomic_rates_mod, only: Haiman_Comp_Heol
use slab_lookup_table_class, only: interpolate_slab_table

use analytic_ionization_solutions_mod, only: analytic_xHIeq
use analytic_ionization_solutions_mod, only: bracket_xHI_plane_parallel

use global_mod, only: CV
use global_mod, only: GV
use global_mod, only: gconst
use hui_gnedin_atomic_rates_mod
use w4_gadget_spline_kernel_class
implicit none
private


public :: tau_particle_type
public :: bckgnd_particle_type

public :: par_to_tau_par
public :: tauHIxHI_from_erf
public :: tauHI_from_erf
public :: tauH_from_erf
public :: initialize_bckgnd_particle
public :: solve_bckgnd_particle_xH_eq






 
!> particle type for solving the analytic ionization equations
!---------------------------------------------------------------
type bckgnd_particle_type

   ! initialized from particle
   !------------------------------
   real(r8b) :: sigmaHIth       !< HI photoionization cross-section @ nu_th
   real(r8b) :: H_mf            !< H mass fraction (mass in H2 + HI + HII) / total mass
   real(r8b) :: fH2             !< H2 mass fraction (mass of H2 / total H mass)
   real(r8b) :: fhot            !< hot phase mass fraction (mass in hot phase / total mass )
   real(r8b) :: T               !< Temperature (warm phase for multi-phase)

   real(r8b) :: xHI             !< Hydrogen neutral fraction = nHI  / nH 
   real(r8b) :: xHII            !< Hydrogen ionized fraction = nHII / nH
   real(r8b) :: hsml            !< smoothing length of particle
   real(r8b) :: eos             !< eos variable (0= never, 1= now, <0= scale fac when off)

   real(r8b) :: mass_cgs        !< total mass in cgs
   real(r8b) :: H2_mass_cgs     !< molecular Hydrogen mass (H2) in cgs
   real(r8b) :: H_mass_cgs      !< atomic Hydrogen mass (HI + HII) in cgs
   real(r8b) :: Hcnt            !< total number of atomic H nuclei

   real(r8b) :: rho_cgs         !< total density in cgs
   real(r8b) :: H2_rho_cgs      !< molecular Hydrogen density (H2) in cgs
   real(r8b) :: H_rho_cgs       !< atomic Hydrogen density (HI + HII) in cgs
   real(r8b) :: nH              !< atomic hydrogen number density


   ! for solution of xHI
   !-------------------------------
   real(r8b) :: GHI             !< H photoionization rate
   real(r8b) :: y               !< n_e = (xHII + y) n_H
   real(r8b) :: tauHI_ray       !< tauHI along ray, exluding neighbors at nu_th
   real(r8b) :: tauH_nbr        !< tauH from neighbors only at nu_th
   real(r8b) :: tauHI_eff       !< tauHI_eff = -ln( Gamma_shld / Gamma_thin )
   logical   :: caseA           !< use case A rates?  if not case B

end type bckgnd_particle_type


!> particle type for calculating optical depth
!---------------------------------------------------------------
type tau_particle_type

   ! initialized from particle
   !------------------------------
   real(r8b) :: sigmaHIth       !< HI photoionization cross-section @ nu_th
   real(r8b) :: H_mf            !< H mass fraction (mass in H2 + HI + HII) / total mass
   real(r8b) :: fH2             !< H2 mass fraction (mass of H2 / total H mass)
   real(r8b) :: fhot            !< hot phase mass fraction (mass in hot phase / total mass )

   real(r8b) :: xHI             !< Hydrogen neutral fraction = nHI  / nH 
   real(r8b) :: xHII            !< Hydrogen ionized fraction = nHII / nH
   real(r8b) :: hsml            !< smoothing length of particle
   real(r8b) :: eos             !< eos variable (0= never, 1= now, <0= scale fac when off)

   real(r8b) :: mass_cgs        !< total mass in cgs
   real(r8b) :: H2_mass_cgs     !< molecular Hydrogen mass (H2) in cgs
   real(r8b) :: H_mass_cgs      !< atomic Hydrogen mass (HI + HII) in cgs
   real(r8b) :: Hcnt            !< total number of atomic H nuclei
   real(r8b) :: HIcnt           !< total number of HI nuclei

   real(r8b) :: rho_cgs         !< total density in cgs
   real(r8b) :: H2_rho_cgs      !< molecular Hydrogen density (H2) in cgs
   real(r8b) :: H_rho_cgs       !< atomic Hydrogen density (HI + HII) in cgs
   real(r8b) :: nH              !< atomic hydrogen number density
   real(r8b) :: nHI             !< HI number density


   ! contribution to column depth
   !-------------------------------
   real(r8b) :: d               !< distance from ray start to particle perp. drop 
   real(r8b) :: b               !< impact parameter
   real(r8b) :: l               !< = sqrt(h^2 - b^2)
   real(r8b) :: bnorm           !< impact parameter / hsml
   real(r8b) :: raylen          !< length of ray 
   real(r8b) :: left            !< distance from ray start to left edge of par path
   real(r8b) :: right           !< distance from ray start to right edge of par path

   real(r8b) :: tauHI_th        !< HI optical depth thru particle at nu_th
   real(r8b) :: tauH_th         !< H optical depth thru particle at nu_th

   real(r8b) :: tau_erf         !< most recently calculated optical depth

end type tau_particle_type




real(r8b), parameter, private :: zero = 0.0d0
real(r8b), parameter, private :: one = 1.0d0
real(r8b), parameter, private :: two = 2.0d0
real(r8b), parameter, private :: four = 4.0d0


contains






!> moves particle variables over to a tau particle
!-----------------------------------------------------------------------
function par_to_tau_par( par, intersection, raylen ) result (tpar)
  type(particle_type), intent(in) :: par
  type(intersection_type), intent(in) :: intersection
  real(r8b), intent(in) :: raylen
  type(tau_particle_type) :: tpar

  logical, save :: first = .true.
  real(r8b), save :: sigmaHI_th

  if (first) then
     sigmaHI_th = Verner_HI_photo_cs(one)    
  endif
  tpar%sigmaHIth = sigmaHI_th
  first = .false.


  tpar%H_mf = return_H_mf( par )
  tpar%fH2  = return_H2_mf( par )
  tpar%eos  = return_eos( par )
  tpar%fhot = return_fhot( par )

  tpar%xHI  = par%xHI
  tpar%xHII = par%xHII
  tpar%hsml = par%hsml
  tpar%nH   = par%nH

  tpar%nHI   = tpar%nH * tpar%xHI

  tpar%mass_cgs    = par%mass * GV%cgs_mass * (one-tpar%fhot)     ! cold & warm gas mass
  tpar%H2_mass_cgs = tpar%mass_cgs * tpar%H_mf * tpar%fH2         ! molecular
  tpar%H_mass_cgs  = tpar%mass_cgs * tpar%H_mf * (one-tpar%fH2)   ! atomic

  tpar%Hcnt  = tpar%H_mass_cgs / gconst%PROTONMASS    ! atomic
  tpar%HIcnt = tpar%Hcnt * tpar%xHI                   ! atomic
 
  tpar%rho_cgs    = par%rho * GV%cgs_rho * (one-tpar%fhot)
  tpar%H2_rho_cgs = tpar%rho_cgs * tpar%H_mf * tpar%fH2
  tpar%H_rho_cgs  = tpar%rho_cgs * tpar%H_mf * (1.0d0 - tpar%fH2)

  tpar%d = intersection%d   
  tpar%b = intersection%b   
  tpar%l = sqrt( tpar%hsml*tpar%hsml - tpar%b*tpar%b )

  tpar%raylen = raylen

  tpar%left  = tpar%d - tpar%l
  tpar%right = tpar%d + tpar%l



end function par_to_tau_par




!> calculates optical depth using error function integrations
!! this assumes that tpar%length is the ray length excluding buffers
!-----------------------------------------------------------------------
subroutine tauHI_from_erf( tpar, cgs_ilen2 )
  type(tau_particle_type) :: tpar
  real(r8b), intent(in) :: cgs_ilen2

  real(r8b) :: zi
  real(r8b) :: zf


  ! particle path completely on ray           
  !---------------------------------------------------
  if ( tpar%left .ge. zero .and. tpar%right .le. tpar%raylen ) then

     zi = -tpar%l
     zf = tpar%l
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%HIcnt * tpar%sigmaHIth * cgs_ilen2
                                  
  ! not on ray at all
  !---------------------------------------------------
  else if ( tpar%right .lt. zero .or. tpar%left .gt. tpar%raylen ) then

     tpar%tau_erf = zero
     return

  ! ray totally inside particle
  !---------------------------------------------------
  else if ( tpar%left .lt. zero .and. tpar%right .gt. tpar%raylen ) then

     zi = -tpar%d
     zf = zi + tpar%raylen
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%HIcnt * tpar%sigmaHIth * cgs_ilen2
     return

  ! partial in front (already cycled for total miss and ray inside par)
  !---------------------------------------------------
  else if ( tpar%left .lt. zero ) then
                 
     zf = tpar%l
     zi = zf - tpar%right 
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%HIcnt * tpar%sigmaHIth * cgs_ilen2   
     return
                       
  ! partial at end (already cycled for total miss and ray inside par)
  !---------------------------------------------------
  else if ( tpar%right .gt. tpar%raylen ) then 
                 
     zi = -tpar%l
     zf = zi + ( tpar%raylen - tpar%left )
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%HIcnt * tpar%sigmaHIth * cgs_ilen2
     return
      
  else
     
     stop 'no branch found for erf solution'
     
  endif

end subroutine tauHI_from_erf




!> calculates optical depth using error function integrations
!! this assumes that tpar%length is the ray length 
!-----------------------------------------------------------------------
subroutine tauH_from_erf( tpar, cgs_ilen2 )
  type(tau_particle_type) :: tpar
  real(r8b), intent(in) :: cgs_ilen2

  real(r8b) :: zi
  real(r8b) :: zf


!  write(*,*) 'par mass cgs  = ', tpar%mass_cgs
!  write(*,*) 'par H mass cgs = ', tpar%H_mass_cgs
!  write(*,*) 'par raylen =     ', tpar%raylen
!  write(*,*) 'par left/right = ', tpar%left, tpar%right

!  stop 'tauH from erf'


  ! particle path completely on ray           
  !---------------------------------------------------
  if ( tpar%left .ge. zero .and. tpar%right .le. tpar%raylen ) then

!     write(*,*) 'ray begin and end outside par'
     zi = -tpar%l
     zf = tpar%l
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%Hcnt * tpar%sigmaHIth * cgs_ilen2
                                  
  ! not on ray at all
  !---------------------------------------------------
  else if ( tpar%right .lt. zero .or. tpar%left .gt. tpar%raylen ) then

!     write(*,*) 'ray miss par'
     tpar%tau_erf = zero
     return

  ! ray totally inside particle
  !---------------------------------------------------
  else if ( tpar%left .lt. zero .and. tpar%right .gt. tpar%raylen ) then

!     write(*,*) 'ray inside par'
     zi = -tpar%d
     zf = zi + tpar%raylen
!     write(*,*) 'zi,zf: ', zi,zf
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%Hcnt * tpar%sigmaHIth * cgs_ilen2
     return

  ! partial in front (already cycled for total miss and ray inside par)
  !---------------------------------------------------
  else if ( tpar%left .lt. zero ) then
                 
!     write(*,*) 'ray starts in par but exits'
     zf = tpar%l
     zi = zf - tpar%right 
!     write(*,*) 'zi,zf: ', zi,zf
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%Hcnt * tpar%sigmaHIth * cgs_ilen2   
     return
                       
  ! partial at end (already cycled for total miss and ray inside par)
  !---------------------------------------------------
  else if ( tpar%right .gt. tpar%raylen ) then 

!     write(*,*) 'ray starts outside of par but ends in it'                 
     zi = -tpar%l
     zf = zi + ( tpar%raylen - tpar%left )
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%Hcnt * tpar%sigmaHIth * cgs_ilen2
     return
      
  else
     
     stop 'no branch found for erf solution'
     
  endif

end subroutine tauH_from_erf



!> calculates optical depth using error function integrations
!! this assumes that tpar%length is the ray length excluding buffers
!-----------------------------------------------------------------------
subroutine tauHIxHI_from_erf( tpar, cgs_ilen2 )
  type(tau_particle_type) :: tpar
  real(r8b), intent(in) :: cgs_ilen2

  real(r8b) :: zi
  real(r8b) :: zf


  ! particle path completely on ray           
  !---------------------------------------------------
  if ( tpar%left .ge. zero .and. tpar%right .le. tpar%raylen ) then

     zi = -tpar%l
     zf = tpar%l
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%HIcnt * tpar%xHI * tpar%sigmaHIth * cgs_ilen2
                                  
  ! not on ray at all
  !---------------------------------------------------
  else if ( tpar%right .lt. zero .or. tpar%left .gt. tpar%raylen ) then

     tpar%tau_erf = zero
     return

  ! ray totally inside particle
  !---------------------------------------------------
  else if ( tpar%left .lt. zero .and. tpar%right .gt. tpar%raylen ) then

     zi = -tpar%d
     zf = zi + tpar%raylen
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%HIcnt * tpar%xHI * tpar%sigmaHIth * cgs_ilen2
     return

  ! partial in front (already cycled for total miss and ray inside par)
  !---------------------------------------------------
  else if ( tpar%left .lt. zero ) then
                 
     zf = tpar%l
     zi = zf - tpar%right 
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%HIcnt * tpar%xHI * tpar%sigmaHIth * cgs_ilen2   
     return
                       
  ! partial at end (already cycled for total miss and ray inside par)
  !---------------------------------------------------
  else if ( tpar%right .gt. tpar%raylen ) then 
                 
     zi = -tpar%l
     zf = zi + ( tpar%raylen - tpar%left )
     tpar%tau_erf = integrate_G3_line( tpar%hsml, tpar%b, zi, zf )
     tpar%tau_erf = tpar%tau_erf * tpar%HIcnt * tpar%xHI * tpar%sigmaHIth * cgs_ilen2
     return
      
  else
     
     stop 'no branch found for erf solution'
     
  endif

end subroutine tauHIxHI_from_erf







!> initializes a particle used just for calculating the xH's
!-----------------------------------------------------------------------
function initialize_bckgnd_particle( par, GHI, tauHI_eff ) result( bpar )
  type(particle_type), intent(in) :: par
  real(r8b), intent(in) :: GHI
  real(r8b), intent(in) :: tauHI_eff
  type(bckgnd_particle_type) :: bpar

  logical, save :: first = .true.
  real(r8b), save :: sigmaHI_th

  if (first) then
     sigmaHI_th = Verner_HI_photo_cs(one)    
  endif
  first = .false.

  bpar%sigmaHIth = sigmaHI_th
  

!  bpar%H_mf = return_H_mf( par )
!  bpar%fH2  = return_H2_mf( par )
!  bpar%eos  = return_eos( par )
!  bpar%fhot = return_fhot( par )

  bpar%GHI       = GHI
  bpar%tauHI_eff = tauHI_eff

  bpar%hsml    = par%hsml
  bpar%T       = par%T
  bpar%xHI     = par%xHI
  bpar%xHII    = par%xHII
  bpar%nH      = par%nH

!  bpar%mass_cgs    = par%mass * GV%cgs_mass * (one-bpar%fhot)     ! cold & warm gas mass
!  bpar%H2_mass_cgs = bpar%mass_cgs * bpar%H_mf * bpar%fH2         ! molecular
!  bpar%H_mass_cgs  = bpar%mass_cgs * bpar%H_mf * (one-bpar%fH2)   ! atomic

!  bpar%Hcnt  = bpar%H_mass_cgs / gconst%PROTONMASS    ! atomic
 
!  bpar%rho_cgs    = par%rho * GV%cgs_rho * (one - bpar%fhot)
!  bpar%H2_rho_cgs = bpar%rho_cgs * bpar%H_mf * bpar%fH2
!  bpar%H_rho_cgs  = bpar%rho_cgs * bpar%H_mf * (one - bpar%fH2)

!  bpar%nH    = bpar%H_rho_cgs / gconst%PROTONMASS

  
end function initialize_bckgnd_particle


!> implements an analytic solution for xHI in particles.
!-----------------------------------------------------------------------
function solve_bckgnd_particle_xH_eq( bpar ) result(err)
 
  type(bckgnd_particle_type) :: bpar
  integer(i4b) :: err
  logical :: err_flag

  bpar%y = zero 
  bpar%xHI = analytic_xHIeq( bpar%T, bpar%GHI, bpar%nH, bpar%y, bpar%tauHI_eff )
  bpar%xHII = 1.0d0 - bpar%xHI

  if (bpar%xHII > one .or. bpar%xHII < zero) then
     err = -1
  else
     err = 0
  endif

end function solve_bckgnd_particle_xH_eq







end module ionpar_mod








module newton_solver_mod
use myf90_mod
use hui_gnedin_atomic_rates_mod
implicit none
private

public :: solve_xHI_low_tau0
public :: solve_xHI_high_tau0
public :: set_newton_constants_med_tau0
public :: newton_xHIeq

integer, parameter :: MAX_ITERATIONS = 10000
real(r8b), parameter :: zero = 0.0d0
real(r8b), parameter :: one = 1.0d0
real(r8b), parameter :: two = 2.0d0

real(r8b) :: T, GHI, nH, y, tau0
real(r8b) :: RC, CI
real(r8b) :: R, Q, P, S

contains


  subroutine solve_xHI_low_tau0( T_in, GHI_in, nH_in, y_in, tau0_in, xHI )
    real(r8b), intent(in) :: T_in, GHI_in, nH_in, y_in, tau0_in
    real(r8b), intent(out) :: xHI
    real(r8b) :: fac, qq, d

    T = T_in
    GHI = GHI_in
    nH = nH_in
    y = y_in
    tau0 = tau0_in

    RC = Hui_HII_recombA( T )
    CI = Hui_HI_col_ion( T )

    R = GHI * tau0 + ( CI + RC ) * nH 
    Q = -( GHI + (CI + two*RC) * nH + (CI + RC) * nH * y ) 
    P = RC * nH * (one + y)

    if (Q*Q < 4.0d0*R*P) then
       write(*,*) 'imbalance'
       write(*,*) 'Q*Q = ', Q*Q
       write(*,*) '4*R*P = ', 4.0d0*R*P
       stop
    endif
    d = sqrt( Q*Q - 4.0d0*R*P )

    qq = - 0.5d0 * ( Q + sign(1.0d0,Q) * d )
    xHI = P / qq

    if (xHI > one)  then
       write(*,*) 'xHI > 1 in low tau'
       write(*,*) 'xHI = ', xHI
       xHI = 1.0d0
       stop
    endif
    
    if (xHI < zero) then
       write(*,*) 'xHI < 0 in low tau'
       write(*,*) 'xHI = ', xHI
       
       write(*,*) 'f(xHI) = ', f(xHI)
       
       xHI = zero
       stop
    endif

  end subroutine solve_xHI_low_tau0

  subroutine solve_xHI_high_tau0( T_in, GHI_in, nH_in, y_in, tau0_in, xHI )
    real(r8b), intent(in) :: T_in, GHI_in, nH_in, y_in, tau0_in
    real(r8b), intent(out) :: xHI
    real(r8b) :: fac, qq, d

    T = T_in
    GHI = GHI_in
    nH = nH_in
    y = y_in
    tau0 = tau0_in

    RC = Hui_HII_recombA( T )
    CI = Hui_HI_col_ion( T )

    if (tau0 > 100.0d0) then
       fac = zero
    else
       fac = exp(-tau0)
    endif
    R = -( GHI * tau0 * fac + ( CI + RC ) * nH ) 
    Q = GHI * fac * (tau0-one) + CI * nH - (CI + RC) * nH * y 
    P = GHI * fac + CI * nH * y
    d = sqrt( Q*Q - 4.0d0*R*P )

    qq = - 0.5d0 * ( Q + sign(1.0d0,Q) * d )
    xHI = qq / R

    if (xHI > one)  then
       write(*,*) 'xHI > 1 in high tau'
       write(*,*) 'xHI = ', xHI
       xHI = 1.0d0
       stop
    endif
    
    if (xHI < zero) then
       write(*,*) 'xHI < 0 in high tau'
       write(*,*) 'xHI = ', xHI
       
       write(*,*) 'f(xHI) = ', f(xHI)
       
       xHI = zero
       stop
    endif

  end subroutine solve_xHI_high_tau0


  subroutine set_newton_constants_med_tau0( T_in, GHI_in, nH_in, y_in, tau0_in )
    real(r8b) :: T_in, GHI_in, nH_in, y_in, tau0_in
    real(r8b) :: fac

    T = T_in
    GHI = GHI_in
    nH = nH_in
    y = y_in
    tau0 = tau0_in

    RC = Hui_HII_recombA( T )
    CI = Hui_HI_col_ion( T )

    R = (CI + RC) * nH
    Q = -( (CI + two*RC) * nH + (CI+RC) * nH * y )
    P = RC * nH * (one+y)
    S = -GHI

  end subroutine set_newton_constants_med_tau0


  ! for tau0 * xHI ~ 1
  !-------------------------------------------
  function f( xHI ) result( y )
    real(r8b) :: xHI, y    
    y = R * xHI*xHI + Q * xHI + P + S * exp(-tau0 * xHI) * xHI
  end function f
  
  function fp( xHI ) result( y )
    real(r8b) :: xHI, y
    y = 2 * R * xHI + Q + S * exp(-tau0 * xHI) * (one - xHI * tau0)
  end function fp
  


  subroutine newton_xHIeq(xHI, tol, iter_flag, neg_flag, pos_flag)
    
    real(r8b), intent(inout) :: xHI
    real(r8b), intent(in) :: tol
    logical, intent(out) :: iter_flag
    logical, intent(out) :: neg_flag
    logical, intent(out) :: pos_flag

    integer :: iter
    
    iter = 1
    iter_flag = .false.
    neg_flag = .false.
    pos_flag = .false.


    if (xHI > one) then 
       write(*,*) 'initial guess high!'
       stop
    endif
    
    if (xHI < zero) then 
       write(*,*) 'initial guess low!'
       stop
    endif

    do 

       xHI = xHI - f(xHI)/fp(xHI)
       if (xHI > one)  then
          pos_flag = .true.
          write(*,*) 'xHI > 1 in newton'
          write(*,*) 'xHI = ', xHI
          write(*,*) 'f(xHI) = ', f(xHI)
          write(*,*) 
          return
       endif
       if (xHI < zero) then 
          neg_flag = .true.
          write(*,*) 'xHI < 0 in newton'
          write(*,*) 'xHI = ', xHI
          write(*,*) 'f(xHI) = ', f(xHI)
          write(*,*) 
          return
       endif

       if ( (abs(f(xHI)) < tol ) ) then                 
          return
       endif

       if (iter > MAX_ITERATIONS) then
          iter_flag = .true.
          return
       endif
       iter = iter + 1

    end do
    
    
  end subroutine newton_xHIeq



end module newton_solver_mod




  










