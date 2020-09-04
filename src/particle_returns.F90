!> \file particle_returns.F90

!> \brief Returns particle quantities that take on different
!! values depending on Makefile flags
!! 
!<

module particle_returns_mod
use myf90_mod
use particle_system_mod, only: particle_type
use config_mod, only: CV
use global_mod, only: GV
use global_mod, only: gconst
implicit none
private


public :: return_H_mf
public :: return_He_mf
public :: return_H2_mf
public :: return_eos
public :: return_fhot
public :: return_nH

real(r8b), parameter :: zero = 0.0d0
real(r8b), parameter :: one  = 1.0d0

contains


!> returns a Hydrogen mass fraction for a single particle
!=================================================================
function return_H_mf( par ) result( H_mf )
  type(particle_type), intent(in)  :: par       
  real(r8b) :: H_mf
  
#ifdef incHmf
  H_mf = par%Hmf
#else     
  H_mf = CV%H_mf
#endif
  
end function return_H_mf


!> returns a Helium mass fraction for a single particle
!=================================================================
function return_He_mf( par ) result( He_mf )
  type(particle_type), intent(in)  :: par       
  real(r8b) :: He_mf
  
#ifdef incHemf
  He_mf = par%Hemf
#else     
  He_mf = CV%He_mf
#endif
  
end function return_He_mf


!> returns a Molecular Hydrogen mass fraction for a single particle
!=================================================================
function return_H2_mf( par ) result( fH2 )
  type(particle_type), intent(in)  :: par       
  real(r8b) :: fH2

#ifdef incH2
  if (CV%IsmAddH2) then     
     fH2 = par%fH2
  else
     fH2 = zero
  endif
#else
  fH2 = zero
#endif  

end function return_H2_mf


!> returns Equation of State variable for a single particle
!=================================================================
function return_eos( par ) result( eos )
  type(particle_type), intent(in)  :: par       
  real(r8b) :: eos
  
#ifdef incEOS
  eos = par%eos
#else
  eos = zero
#endif
  
end function return_eos


!> returns Hot mass fraction for a single particle
!=================================================================
function return_fhot( par ) result( fhot )
  type(particle_type), intent(in)  :: par       
  real(r8b) :: fhot
  
#ifdef incEOS
  if (CV%AdjustIsmPars) then
     fhot = CV%IsmHotFrac
  else
     fhot = zero
  endif
#else
  fhot = zero
#endif
  
end function return_fhot



!> returns atomic Hydrogen number density
!=================================================================
function return_nH( par, h, a, adjust_ism ) result( nH )
  type(particle_type), intent(in)  :: par       
  real(r8b), intent(in) :: h      ! little Hubble 
  real(r8b), intent(in) :: a      ! scale factor
  logical, intent(in) :: adjust_ism
  real(r8b) :: nH
  real(r8b) :: H_mf
  real(r8b) :: fH2
  real(r8b) :: fhot

  H_mf = return_H_mf( par )        ! H mass / total mass
  nH = par%rho * H_mf * GV%cgs_rho ! [g/cm^3 h^2 / a^3] 

  nH = nH * (h*h) / (a*a*a)

  nH = nH / gconst%PROTONMASS      ! [atoms/cm^3]
  
  if (adjust_ism) then
     fhot = return_fhot( par )
     fH2 = return_H2_mf( par )
     nH = nH * (one - fhot) * (one - fH2)
  endif
     
end function return_nH



  







end module particle_returns_mod
