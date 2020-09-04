!> \file main_input.F90

!> \brief The module that calls the specific input routines
!<

module main_input_mod
use myf90_mod
use gadget_public_input_mod
use gadget_cosmoBH_input_mod
use gadget_owls_input_mod
use gadget_eagle_input_mod
use gadget_vbromm_input_mod
use gadget_public_input_hdf5_mod
use particle_system_mod, only: particle_type
use particle_system_mod, only: box_type
use particle_system_mod, only: scale_comoving_to_physical
use particle_system_mod, only: scale_physical_to_comoving
use particle_system_mod, only: set_ye
use particle_system_mod, only: enforce_x_and_T_minmax
use particle_system_mod, only: particle_info_to_screen
use particle_system_mod, only: constrain_psys_to_selection
use particle_returns_mod, only: return_nH
use particle_returns_mod, only: return_H_mf
use particle_returns_mod, only: return_He_mf
use atomic_rates_mod, only: get_atomic_rates
use config_mod, only: CV
use global_mod, only: GV
use global_mod, only: psys 
use global_mod, only: gconst
use global_mod, only: h1_itab
use ion_table_class, only: return_gammaHI_at_z
implicit none
private

public :: get_planning_data
public :: readin_snapshot


real(r8b), parameter :: zero = 0.0d0
real(r8b), parameter :: one = 1.0d0

contains


!> Read in planning data from the headers of all snapshots 
!===========================================================
subroutine get_planning_data()  
  character(clen), parameter :: myname="get_planning_data"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=1

  call mywrite("getting planning data:", verb)
  call mywrite("",verb) 

  ! branch on input type
  !-------------------------------------------------------------------
  select case (trim(CV%InputType))
  case('public')
     call get_planning_data_gadget_public()
  case('cosmobh')
     call get_planning_data_gadget_cosmoBH()
  case('owls')
     call get_planning_data_gadget_owls()
  case('vbromm')
     call get_planning_data_gadget_vbromm()
  case('publichdf5')
     call get_planning_data_gadget_public_hdf5()
  case('eagle')
     call get_planning_data_gadget_eagle()
  end select


end subroutine get_planning_data


!> read in particle and box data 
!============================================
subroutine readin_snapshot()

  character(clen), parameter :: myname="readin_snapshot"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt
  
  real(r8b) :: MB
  integer(i8b) :: ipar
  real(r8b) :: a     !< scale factor
  real(r8b) :: h     !< Hubble paraemter (little H)
  character(clen) :: snpbase

  real(r8b) :: H_mf
  real(r8b) :: He_mf
  real(r8b) :: nH

  
  call mywrite("reading in particle snapshots:", verb-1)
  call mywrite("",verb-1) 

  ! set local variables
  !======================
  a = GV%ScaleFac
  h = GV%LittleH

  ! read in the particle data
  !=============================
  call mywrite('   input type = ', verb, adv=.false.)

  select case (trim(CV%InputType))
  case('public')
     call mywrite(" Gadget Public (SnapFormat=1)", verb)
     call read_Gpublic_particles()

  case('cosmobh')
     call mywrite(" Gadget CosmoBH", verb)
     call read_GcosmoBH_particles()

  case('owls')
     call mywrite(" Gadget OWLS/GIMIC", verb)
     call read_Gowls_particles()

  case('vbromm')
     call mywrite(" Gadget V. Bromm", verb)
     call read_Gvbromm_particles()

  case('publichdf5')
     call mywrite(" Gadget Public HDF5", verb)
     call read_Gpubhdf5_particles()

  case('eagle')
     call mywrite(" Gadget EAGLE", verb)
     call read_Geagle_particles()

  case default
     write(str,*) "input type, ", CV%InputType, "not recognized" 
     call myerr(str,myname,crash)

  end select



  ! default output mask to TRUE
  !==================================================
  psys%par(:)%mask = .true. 


  ! copy over box properties
  !==================================================
  psys%box%top = GV%BoxUprsComoh
  psys%box%bot = GV%BoxLwrsComoh

  psys%box%tbound(1) = CV%XcoordBndryCond
  psys%box%bbound(1) = CV%XcoordBndryCond

  psys%box%tbound(2) = CV%YcoordBndryCond
  psys%box%bbound(2) = CV%YcoordBndryCond

  psys%box%tbound(3) = CV%ZcoordBndryCond
  psys%box%bbound(3) = CV%ZcoordBndryCond



  ! use selected region if we have one
  !==================================================
#ifndef EAGLE
  if (CV%DoSelection) then
     call constrain_psys_to_selection(psys)
  endif
#endif

  ! initialize to zero
  !===============================================================
#ifdef outGammaHI
  psys%par(:)%gammaHI = zero
#endif


  ! calculate atomic nH without adjustments for log output
  !=======================================================
  do ipar = 1, size(psys%par)
     psys%par(ipar)%nH = return_nH( psys%par(ipar), h, a, adjust_ism=.false. )
  end do


  ! set par file bases for output to logfiles
  !====================================================
  fmt = "(A,'/',A)"
  write(snpbase,fmt) trim(CV%SnapPath),   trim(CV%ParFileBase)


  ! write fresh reads to the particle_data log file
  !========================================================================  
  fmt = "(A,A)"
  write(str,fmt) "Fresh read from ", trim(snpbase)
  call particle_info_to_screen(psys,str,GV%pardatalun)
  write(GV%pardatalun,*)
  write(GV%pardatalun,*)
  flush(GV%pardatalun)


  ! re-calculate atomic nH with ISM adjustments
  !==================================================
  do ipar = 1, size(psys%par)
     psys%par(ipar)%nH = return_nH( psys%par(ipar), h, a, adjust_ism=CV%AdjustIsmPars )
  end do





 
  ! set constant temperature if we have one
  !=======================================================
  if (CV%IsoTemp > 0.0) psys%par(:)%T = CV%IsoTemp


  ! scale the mass / density if we want
  !=======================================================
  if (CV%RhoMultiplier > 0.0d0) then
     psys%par(:)%rho  = psys%par(:)%rho * CV%RhoMultiplier
     psys%par(:)%mass = psys%par(:)%mass * CV%RhoMultiplier
     psys%par(:)%nH  = psys%par(:)%nH * CV%RhoMultiplier
  endif



  ! write data after above conditionals to the particle log file
  !==========================================================================
  fmt = "(A,A)"
  write(str,fmt) "After test conditionals from ", trim(snpbase)
  call particle_info_to_screen(psys,str,GV%pardatalun)
  write(GV%pardatalun,*)
  write(GV%pardatalun,*)
  flush(GV%pardatalun)
 

  ! scale the data if we need to
  !=====================================================================
  if(CV%Comoving) then
     call scale_comoving_to_physical(psys, a, h )
  endif


  ! write data after rescaling to the particle_data log file
  !=====================================================================
  fmt = "(A,F5.3,A,F5.3,A,A)"
  write(str,fmt) "After rescaling (a=",a,",h=",h,") from ", trim(snpbase)
  call particle_info_to_screen(psys,str,GV%pardatalun)
  write(GV%pardatalun,*)
  write(GV%pardatalun,*)
  flush(GV%pardatalun)


 
  ! and the rest of the stuff
  !===============================================================  
  GV%UVB_gammaHI_cloudy = return_gammaHI_at_z( h1_itab, GV%RedShift )


  GV%total_mass = zero
  GV%total_H_atoms = zero
  GV%total_He_atoms = zero

  do ipar = 1,size(psys%par)

     GV%total_mass = GV%total_mass + psys%par(ipar)%mass

     H_mf = return_H_mf( psys%par(ipar) )
     GV%total_H_atoms = GV%total_H_atoms + &
          psys%par(ipar)%mass * GV%cgs_mass * H_mf  / gconst%protonmass

     He_mf = return_He_mf( psys%par(ipar) )    
     GV%total_He_atoms = GV%total_He_atoms + &
          psys%par(ipar)%mass * GV%cgs_mass * He_mf / (4.0d0*gconst%protonmass)

  end do


  call mywrite("",verb)


 
end subroutine readin_snapshot





end module main_input_mod
