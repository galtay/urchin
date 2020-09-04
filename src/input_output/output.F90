!> \file output.F90

!> \brief The module that handles output
!<

module output_mod
! libs
use myf90_mod
! types
use particle_system_mod, only: particle_system_type
! routines
use particle_system_mod, only: scale_comoving_to_physical
use particle_system_mod, only: scale_physical_to_comoving
use particle_system_mod, only: set_ye
use output_gadget_public_mod, only: output_snap_gadget_public
use output_gadget_hdf5_mod, only: output_snap_gadget_hdf5
! variables
use config_mod, only: CV
use global_mod, only: GV
use global_mod, only: saved_gheads, snap_file_names
use config_mod, only: write_config_hdf5_lun

implicit none
private



public :: output_total_snap


contains
 
!=======================================================================
!=======================================================================
!
! OUTPUT ROUTINES
!
!=======================================================================
!=======================================================================






!> outputs the whole shebang to a file
!=======================================

  subroutine output_total_snap(psys)

     character(clen), parameter :: myname="output_total_snap"
     logical, parameter :: crash = .true.
     integer, parameter :: verb = 2

     type(particle_system_type), intent(inout) :: psys

     integer(i4b) :: NumFiles

     real(r8b) :: scale
     real(r8b) :: hub


     ! set local variables
     !==============

     scale = GV%ScaleFac
     hub   = GV%LittleH


     ! output to single file if this is a selection
     ! this must be 1 or CV%ParFilesPerSnap
     !------------------------------------------------
     if (CV%DoSelection) then
        NumFiles=1
     else          
        NumFiles = CV%ParFilesPerSnap
     end if


     write(*,*) 'output number: ', GV%OutputIndx
     write(*,*) 'output type:   ', trim(CV%OutputType)
     write(*,*) 'output a selection?: ', CV%DoSelection
     write(*,*) 'output Nfiles:       ', NumFiles
     write(*,*) "writing total state of system"


     if (CV%Comoving) then
        call scale_physical_to_comoving(psys, scale, hub)
     endif
     

     call set_ye( psys )


     !================================================================
     ! GADGET public formatted output
     !================================================================
     if (trim(CV%OutputType) == 'binary') then

        call output_snap_gadget_public(psys, NumFiles)

     end if
 

     !================================================================     
     ! GADGET HDF5 formatted output 
     !================================================================
     if (trim(CV%OutputType) == 'hdf5') then

        call output_snap_gadget_hdf5(psys,NumFiles)
      
     end if
 

     !================================================================
     !================================================================
     if (CV%Comoving) then
        call scale_comoving_to_physical(psys, scale, hub)
     endif


  end subroutine output_total_snap



  





end module output_mod
