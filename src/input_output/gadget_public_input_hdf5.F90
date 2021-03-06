!> \file gadget_public_input_hdf5.F90

!> \brief Handles readin of GADGET Public HDF5 formatted files
!<

module gadget_public_input_hdf5_mod
use myf90_mod
use hdf5_wrapper
use gadget_general_class
use gadget_public_header_class
use gadget_public_header_hdf5_class
use gadget_sphray_header_class
use particle_system_mod
use gadget_public_input_mod, only: set_temp_from_u

use config_mod, only: CV
use global_mod, only: psys, GV
use global_mod, only: saved_gheads, snap_file_names





implicit none
private

public :: get_planning_data_gadget_public_hdf5
public :: read_Gpubhdf5_particles


contains



  !>   gets run planning data from Gadget Public HDF5 Headers
  !=============================================================
  subroutine get_planning_data_gadget_public_hdf5()
    
    character(clen), parameter :: myname = 'get_planning_data_gadget_public_hdf5'
    logical, parameter :: crash = .true.
    integer, parameter :: verb = 2
    
    type(gadget_public_header_type) :: ghead
    type(gadget_units_type) :: gunits
    type(gadget_constants_type) :: gconst
    

    integer(i4b) :: pfiles          ! files/snap for particles    
    integer(i4b) :: i,j             ! counters
    character(clen) :: snapfile     ! snapshot file name
    
    integer(i4b) :: loglun
    character(clen) :: logfile
    

    
    
    ! these global variables are read from the config file
    !======================================================
    pfiles = CV%ParFilesPerSnap
    
    if ( allocated(saved_gheads) ) deallocate(saved_gheads)
    allocate( saved_gheads(0:pfiles-1) )
    
    if ( allocated(snap_file_names) ) deallocate(snap_file_names)
    allocate( snap_file_names(0:pfiles-1) )

    ! set global units
    !===================================================
    GV%cgs_len  = gunits%cgs_length
    GV%cgs_mass = gunits%cgs_mass
    GV%cgs_vel  = gunits%cgs_velocity
    GV%cgs_time = gunits%cgs_time
    GV%cgs_rho  = gunits%cgs_density
    GV%cgs_prs  = gunits%cgs_pressure
    GV%cgs_enrg = gunits%cgs_energy
    
    ! read all particle headers and write to log file
    !===================================================
    do j = 0,pfiles-1
          
       call form_gadget_snapshot_file_name(CV%SnapPath,CV%ParFileBase,j,snapfile,hdf5bool=.true.)
       snap_file_names(j) = trim(snapfile)

       call gadget_public_header_hdf5_read_file(ghead, snapfile )
       call gadget_sphray_header_copy_public( saved_gheads(j), ghead )
          
       saved_gheads(j)%OmegaB = 0.045 ! this should be moved to the config file
                                      ! but its not used in the code now
          
       saved_gheads(j)%time_gyr = gadget_public_header_return_gyr(ghead)
          
       ! make sure there is gas in this snapshot
       !-------------------------------------------
       if (.not. ghead%npar_all(0) > 0) then
          write(*,*) "Gadget snapshot does not contain any gas particles"
          write(*,*) "Sphray cannot read dark matter particles directly, "
          write(*,*) "please calculate smoothing lengths for these particles"
          write(*,*) "and write them as gas particles. "
          stop
       end if
                    
    end do
    

    ! use one header to set global variables
    !===================================================
    GV%OmegaM   = saved_gheads(0)%OmegaM
    GV%OmegaL   = saved_gheads(0)%OmegaL
    GV%OmegaB   = saved_gheads(0)%OmegaB 
    GV%LittleH  = saved_gheads(0)%h
    GV%TimeGyr  = saved_gheads(0)%time_gyr
    GV%RedShift = saved_gheads(0)%z
    GV%ScaleFac = saved_gheads(0)%a
    GV%Tcmb_cur = gconst%t_cmb0 / saved_gheads(0)%a
    GV%RhoConv  = GV%LittleH * GV%LittleH / (GV%ScaleFac * GV%ScaleFac * GV%ScaleFac)
    
    GV%BoxLwrsComoh(:) = 0.0d0
    GV%BoxUprsComoh(:) = saved_gheads(0)%boxlen    
    
       
    ! write units to log file
    !===================================================
    
    logfile = trim(CV%OutputDir) // "/" // "code_units.log"
    call open_formatted_file_w(logfile,loglun)
    call gadget_units_print_lun(gunits, loglun, saved_gheads(0)%h)
    close(loglun)
    
    
  end subroutine get_planning_data_gadget_public_hdf5






  
  !> reads a Gadget Public HDF5 snapshot into a particle array  
  !========================================================================
  subroutine read_Gpubhdf5_particles()
    
    character(clen), parameter :: myname="read_Gpubhdf5_particles" 
    logical, parameter :: crash=.true.
    integer, parameter :: verb=2
    character(clen) :: str,fmt
    
    real(r4b), allocatable :: rblck(:)
    real(r4b), allocatable :: rblck3(:,:)
    integer(i4b), allocatable :: iblck(:)
    integer(i8b) :: ngasread
    
    type(gadget_public_header_type) :: ghead
    type(gadget_sphray_header_type) :: shead
    character(clen) :: snapfile, VarName, GroupName
    integer(i4b) :: fh
    integer(i8b) :: i
    integer(i4b) :: err
    
    integer(i8b) :: npar, ngas, nmass
    integer(i8b) :: npar1, ngas1, nmass1
    logical :: varmass(0:5)
    integer(i4b) :: fn
    
    real(r8b) :: meanweight
    logical :: caseA(2)
    real(r8b) :: MB 
    
    logical :: hdf5bool
    
    ! set hdf5 boolean
    !======================================================
    hdf5bool = .true.
    
    ! set local particle numbers
    !============================
    shead = saved_gheads( 0 )
    varmass = (shead%npar_all > 0 .and. shead%mass == 0)
    npar = sum(shead%npar_all)
    ngas = shead%npar_all(0)
    nmass = sum(shead%npar_all, mask=varmass)
    
    ! do Gadget dummy checks
    !============================
    if (ngas == 0) call myerr("snapshot has no gas particles",myname,crash)
    
    ! calculate bytes per particle and allocate particle array
    !===========================================================
    MB = GV%bytesperpar * real(ngas) / 2.0d0**20
    GV%MB = GV%MB + MB
    
    fmt="(A,F10.4,A,I10,A)"
    write(str,fmt) "   allocating ", MB, " MB for ", ngas, " particles"
    call mywrite(str,verb) 
    
    allocate (psys%par(ngas), stat=err)
    if (err /= 0) call myerr("failed to allocate par",myname,crash)
    
    
    ! now read all snapshot files
    !==============================          
    ngasread = 0
    GroupName = 'PartType0/'
    files: do fn = 0, shead%nfiles-1
       
       
       ! recall the header info
       !-----------------------------------------------------------!  
       shead   = saved_gheads( fn )
       varmass = (shead%npar_file > 0 .and. shead%mass == 0)
       npar1   = sum(shead%npar_file)
       ngas1   = shead%npar_file(0)
       nmass1  = sum(shead%npar_file, mask=varmass)
       if (ngas1 == 0) cycle
       
       
       ! begin read
       !-----------------------------------------------------------!  
       call form_gadget_snapshot_file_name(CV%SnapPath,CV%ParFileBase,fn,snapfile,hdf5bool)
       call mywrite("   reading particle snapshot file: "//trim(snapfile), verb)
       call hdf5_open_file(fh, snapfile, readonly=.true.)
       call gadget_public_header_hdf5_read_lun(ghead,fh)

       
       ! read positions 
       !-----------------------------------------------------------!  
       allocate(rblck3(3,ngas1), stat=err)
       if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
       VarName = 'Coordinates'
       call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)
       forall(i=1:ngas1) psys%par(ngasread+i)%pos(1) = rblck3(1,i)
       forall(i=1:ngas1) psys%par(ngasread+i)%pos(2) = rblck3(2,i)
       forall(i=1:ngas1) psys%par(ngasread+i)%pos(3) = rblck3(3,i)
       deallocate(rblck3)
       
       ! read velocities 
       !-----------------------------------------------------------!  
#ifdef incVel
       allocate(rblck3(3,ngas1), stat=err)
       if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
       VarName = 'Velocity'
       call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)
       forall(i=1:ngas1) psys%par(ngasread+i)%vel(1) = rblck3(1,i)
       forall(i=1:ngas1) psys%par(ngasread+i)%vel(2) = rblck3(2,i)
       forall(i=1:ngas1) psys%par(ngasread+i)%vel(3) = rblck3(3,i)
       deallocate(rblck3)
#endif
       
       ! read id's 
       !-----------------------------------------------------------!  
       allocate(iblck(ngas1), stat=err )
       if(err/=0) call myerr("allocating iblck for ID",myname,crash)
       VarName = 'ParticleIDs'
       call hdf5_read_data(fh,trim(GroupName)//trim(VarName),iblck)
       forall(i=1:ngas1) psys%par(ngasread+i)%id = iblck(i)
       deallocate(iblck)
       
       ! read masses 
       !-----------------------------------------------------------!  
       
       ! if gas particles are variable mass
       if (varmass(0)) then  
          allocate(rblck(ngas1), stat=err)
          if(err/=0) call myerr("allocating rblck for mass",myname,crash)
          VarName = 'Mass'
          call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
          forall(i=1:ngas1) psys%par(ngasread+i)%mass = rblck(i)
          deallocate(rblck)
          
          ! if gas particles are isomass
       else
          psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(0)       
       end if
       
       ! read temperature (internal energy / unit mass for now)
       !-----------------------------------------------------------!  
       if (CV%IsoTemp <= 0.0d0) then
          allocate(rblck(ngas1), stat=err)
          if(err/=0) call myerr("allocating rblck for u",myname,crash)
          VarName = 'InternalEnergy'
          call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
          forall(i=1:ngas1) psys%par(ngasread+i)%T = rblck(i)
          deallocate(rblck)
       endif
       
       ! read density 
       !-----------------------------------------------------------!  
       allocate(rblck(ngas1), stat=err)
       if(err/=0) call myerr("allocating rblck for rho",myname,crash)
       VarName = 'Density'
       call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
       forall(i=1:ngas1) psys%par(ngasread+i)%rho = rblck(i)
       deallocate(rblck)
       
       
       ! read smoothing lengths 
       !-----------------------------------------------------------!  
       allocate(rblck(ngas1), stat=err)
       if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
       VarName = 'SmoothingLength'
       call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
       forall(i=1:ngas1) psys%par(ngasread+i)%hsml = rblck(i)
       deallocate(rblck)
       
       
       ngasread = ngasread + ngas1
       call hdf5_close_file(fh)
       
       
    end do files
    
    ! There is no cooling in the public version of Gadget-2 so we
    ! will make the same assumptions as were made in the routine
    ! show_hotgas.pro included in the source distribution when 
    ! calculating temperature from internal energy / unit mass. 
    ! Specifically mu = 1.0 * PROTONMASS or fully neutral gas. 
    !---------------------------------------------------------------
    psys%par(:)%xHI = 1.0d0
    psys%par(:)%xHII = 0.0d0
    psys%par(:)%ye = psys%par(:)%xHII

    if (CV%IsoTemp <= 0.0d0) then
       call set_temp_from_u(psys, CV%H_mf, GV%cgs_enrg, GV%cgs_mass)
    else
       psys%par(:)%T = CV%IsoTemp
    endif

    ! now we have temperatures, we will set the ionization
    ! state to collisional equilibrium.
    !-----------------------------------------------------
    caseA = .true.
    call set_collisional_ionization_equilibrium( psys, caseA, CV%IsoTemp, DoHydrogen=.true., fit='hui' )
    call set_ye( psys )
            

    write (*,*) ' pareticles read in using read_Gpubhdf5_particles'
    write (*,*) ' rho  range (gadget units) ' ,minval(psys%par(:)%rho),maxval(psys%par(:)%rho)
    write (*,*) ' Temp range [K] ' ,minval(psys%par(:)%T),maxval(psys%par(:)%T)
    write (*,*) ' log(xHI) range ' ,log10(minval(psys%par(:)%xHI)),log10(maxval(psys%par(:)%xHI))

  end subroutine read_Gpubhdf5_particles
  
  

  
  
end module gadget_public_input_hdf5_mod
