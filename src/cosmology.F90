!> \file cosmology.F90

!> \brief Module contains simple analytic cosmological functions
!<
module cosmology_mod
use myf90_mod
implicit none
private

public :: tsinceBB
public :: da2dt

real(r8b), parameter :: one_third = 1.0d0 / 3.0d0
real(r8b), parameter :: three_halves = 3.0d0 / 2.0d0
real(r8b), parameter :: one = 1.0d0
real(r8b), parameter :: two = 2.0d0
real(r8b), parameter :: three = 3.0d0


  contains


!> returns the time since the big bang in a flat LCDM Universe
!! in units of the Hubble time (t_H = 1/H0) 
!-------------------------------------------------------------------
  function tsinceBB(a, OmegaM) result (t)
    real(r8b), intent(in) :: a      !< scale factor
    real(r8b), intent(in) :: OmegaM !< Omega Matter
    real(r8b) :: t                  !< time since big bang (t * H0)

    real(r8b) :: OmegaL,aeq
    real(r8b) :: pre,arg

    OmegaL = one - OmegaM
    aeq    = (OmegaM / OmegaL)**(one_third)

    pre = two / ( three * sqrt(OmegaL) )
    arg = (a/aeq)**(three_halves) + sqrt(one + (a/aeq)**3)
    t = pre * log(arg)

  end function tsinceBB



!> returns the time between two scale factors in a flat LCDM Universe
!! in units of the Hubble time (t_H = 1/H0)
!-------------------------------------------------------------------
  function da2dt(a1, a2, OmegaM) result (dt) 
    real(r8b), intent(in) :: a1     !< initial scale factor
    real(r8b), intent(in) :: a2     !< final scale factor
    real(r8b), intent(in) :: OmegaM !< Omega Matter
    real(r8b) :: dt                 !< time interval (t * H0)
    
    real(r8b) :: t1,t2
    
    if (OmegaM <= 0.0) stop "cosmology.f90>> OmegaMatter must be >= 0.0"

    if (a1 <= 0.0) then
       write(*,*) "cosmology.f90>> da2dt - a1 out of range "
       write(*,*) "a1 = ", a1
       stop
    end if

    if (a2 <= 0.0) then
       write(*,*) "cosmology.f90>> da2dt - a2 out of range "
       write(*,*) "a2 = ", a2
       stop
    end if

    if (a2 < a1) then
       write(*,*) "cosmology.f90>> da2dt -  a2 must be > a1"
       write(*,*) "a1,a2 = ", a1, a2
       stop
    end if

    t1 = tsinceBB(a1,OmegaM)
    t2 = tsinceBB(a2,OmegaM)
    dt = t2-t1
 
  end function da2dt



end module cosmology_mod


