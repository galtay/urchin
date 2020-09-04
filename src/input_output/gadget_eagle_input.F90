!> \file gadget_eagle_input.F90

!> \brief Handles readin of GADGET EAGLE HDF5 formatted files
!<

module gadget_eagle_input_mod

#ifdef EAGLE
use read_eagle
#endif

! libs
use myf90_mod
use hdf5_wrapper

! types
use gadget_general_class, only: gadget_units_type
use gadget_general_class, only: gadget_constants_type
use gadget_general_class, only: gadget_owls_parameters_type
use gadget_owls_header_class, only: gadget_owls_header_type
use gadget_sphray_header_class, only: gadget_sphray_header_type

! routines
use ion_table_class, only: return_gammaHI_at_z
use particle_system_mod, only: set_ye
use particle_system_mod, only: constrain_psys_to_selection
use gadget_general_class, only: form_gadget_snapshot_file_name
use gadget_general_class, only: gadget_constants_read_file
use gadget_general_class, only: gadget_units_read_file
use gadget_general_class, only: gadget_units_print_lun
use gadget_general_class, only: gadget_data_attributes_read_lun
use gadget_general_class, only: gadget_owls_parameters_read_file
use particle_returns_mod, only: return_nH
use gadget_owls_header_class, only: gadget_owls_header_read_file
use gadget_owls_header_class, only: gadget_owls_header_read_lun
use gadget_sphray_header_class, only: gadget_sphray_header_copy_owls
use analytic_ionization_solutions_mod, only: analytic_xHIeq

! variables
use config_mod, only: CV
use global_mod, only: GV
use global_mod, only: psys
use global_mod, only: saved_gheads
use global_mod, only: snap_file_names
use global_mod, only: h1_itab
use global_mod, only: global_owls_parameters
use gadget_general_class, only: gattrs

implicit none
private

real(r8b) :: zero = 0.0d0
real(r8b) :: one = 1.0d0

public :: get_planning_data_gadget_eagle
public :: read_Geagle_particles


contains






!>   gets run planning data from Gadget EAGLE Headers
!============================================================
subroutine get_planning_data_gadget_eagle()

  character(clen), parameter :: myname = 'get_planning_data_gadget_eagle'
  logical, parameter :: crash = .true.
  integer, parameter :: verb = 2

  type(gadget_owls_header_type) :: ghead
  type(gadget_units_type) :: gunits
  type(gadget_constants_type) :: gconst
  type(gadget_owls_parameters_type) :: gparam

  integer(i4b) :: iSnap           ! initial and final snapshot numbers
  integer(i4b) :: pfiles          ! files/snap for particles    
  integer(i4b) :: j               ! counters
  character(clen) :: snapfile     ! snapshot file name

  integer(i4b) :: fh
  character(clen) :: logfile
  integer(i4b) :: loglun
  
  ! these global variables are read from the config file
  !======================================================
  pfiles = CV%ParFilesPerSnap

  if ( allocated(saved_gheads) ) deallocate(saved_gheads)
  allocate( saved_gheads(0:pfiles-1) )

  if ( allocated(snap_file_names) ) deallocate(snap_file_names)
  allocate( snap_file_names(0:pfiles-1) )

  ! set global units and read constants
  !===================================================
  call form_gadget_snapshot_file_name(CV%SnapPath,CV%ParFileBase,0,snapfile,hdf5bool=.true.)
  call gadget_constants_read_file( gconst, snapfile )
  call gadget_units_read_file( gunits, snapfile )
  call gadget_owls_parameters_read_file( gparam, snapfile )

  GV%sf_gamma_eos = gparam%stellar_evolution_parameters%SF_EOSGammaEffective
  GV%nH_star = gparam%stellar_evolution_parameters%SF_THRESH_MinPhysDens_HpCM3

  global_owls_parameters = gparam

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
     snap_file_names(j) = snapfile

     call gadget_owls_header_read_file( ghead, snapfile )
     call gadget_sphray_header_copy_owls( saved_gheads(j), ghead )
     
     ! make sure there is gas in this snapshot
     if (.not. ghead%npar_all(0) > 0) then
        write(*,*) "Gadget snapshot does not contain any gas particles"
        write(*,*) "Urchin cannot read dark matter particles directly, "
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

  GV%UVB_gammaHI_cloudy = return_gammaHI_at_z( h1_itab, saved_gheads(0)%z )
  
  ! write units to log file
  !===================================================

  logfile = trim(CV%OutputDir) // "/" // "code_units.log"
  call open_formatted_file_w(logfile,loglun)
  call gadget_units_print_lun( gunits, loglun, saved_gheads(0)%h ) 
  close(loglun)  

  
end subroutine get_planning_data_gadget_eagle




subroutine read_Geagle_particles
  implicit none

#ifndef EAGLE

  stop ' EAGLE flag not on in Makefile '

#else
  
  character(clen), parameter :: myname="read_Geagle_particles" 
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  real(r4b), allocatable :: rblck(:)
  real(r8b), allocatable :: rblck3(:,:)
  integer(i8b), allocatable :: iblck(:)
 
  type(gadget_owls_header_type) :: ghead
  type(gadget_sphray_header_type) :: shead
 
  character(len=clen) :: fname
  character(clen) :: VarName, GroupName
  integer(i8b) :: ngas, np
  real(r8b) :: MB
  integer(i8b) :: ipar
  integer(i4b) :: err
  integer(i4b) :: fh
  type (eaglesnapshot) :: snap
  integer(i8b) :: i
  real(r4b) :: xmin, xmax, ymin, ymax, zmin, zmax

  real(r8b) :: meanweight
  logical :: caseA(2)

  real(r8b) :: ratio
  real(r8b) :: mu
  real(r8b) :: prs
  real(r8b) :: eos_prs_norm

  real(r8b) :: Tdum
  real(r8b) :: Hmf
  real(r8b) :: nH8
  real(r8b) :: T8
  real(r8b) :: xvec(5)

  integer(i4b), allocatable :: file_offset(:)
  integer(i4b), allocatable :: file_index(:)

  type(gadget_constants_type) :: gconst

  fname = trim(CV%SnapPath)//'/'//trim(CV%ParFileBase)//'.0.hdf5'


  shead = saved_gheads( 0 )


  ! read dataset attributes before anything else
  !-----------------------------------------------------------!
  GroupName = 'PartType0/'
  call hdf5_open_file(fh, fname, readonly=.true.)

  VarName = 'Coordinates'
  call gadget_data_attributes_read_lun( gattrs%pos, fh, GroupName, VarName )

  VarName = 'Velocity'
  call gadget_data_attributes_read_lun( gattrs%vel, fh, GroupName, VarName )

  VarName = 'ParticleIDs'
  call gadget_data_attributes_read_lun( gattrs%id, fh, GroupName, VarName )

  VarName = 'Mass'
  call gadget_data_attributes_read_lun( gattrs%mass, fh, GroupName, VarName )

  VarName = 'Temperature'
  call gadget_data_attributes_read_lun( gattrs%T, fh, GroupName, VarName )

  VarName = trim(CV%EOSvarName)
  call gadget_data_attributes_read_lun( gattrs%eos, fh, GroupName, VarName )

  VarName = 'Density'
  call gadget_data_attributes_read_lun( gattrs%rho, fh, GroupName, VarName )

  VarName = 'ElementAbundance/Hydrogen'
  call gadget_data_attributes_read_lun( gattrs%Hmf, fh, GroupName, VarName )

  VarName = 'ElementAbundance/Helium'
  call gadget_data_attributes_read_lun( gattrs%Hemf, fh, GroupName, VarName )

  
  ! create snapshot object using jch code
  !-----------------------------------------------------------!
  snap = open_snapshot(fname)

  write(*,*) "Boxsize   = ", snap%boxsize
  write(*,*) "Numfiles  = ", snap%numfiles
  write(*,*) "Hashbits  = ", snap%hashbits
  write(*,*) "NumPart   = ", snap%numpart_total

  if ( CV%DoSelection ) then
     xmin = CV%SelectXcoord - CV%SelectRadius
     xmax = CV%SelectXcoord + CV%SelectRadius
     ymin = CV%SelectYcoord - CV%SelectRadius
     ymax = CV%SelectYcoord + CV%SelectRadius
     zmin = CV%SelectZcoord - CV%SelectRadius
     zmax = CV%SelectZcoord + CV%SelectRadius
  else
     xmin = 0.0d0
     xmax = snap%boxsize
     ymin = 0.0d0
     ymax = snap%boxsize
     zmin = 0.0d0
     zmax = snap%boxsize
  end if

  write(*,*) 'Eagle Region: '
  write(*,*) '  xmin/max: ', xmin, xmax
  write(*,*) '  ymin/max: ', ymin, ymax
  write(*,*) '  zmin/max: ', zmin, zmax

  call select_region(snap, xmin, xmax, ymin, ymax, zmin, zmax )
  ngas = count_particles(snap, 0)

  write(*,*) '  ngas: ', ngas
  write(*,*) 

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

  ! store file and file location for each particle
  allocate( file_index(ngas), file_offset(ngas) )
  np = get_particle_locations(snap, 0, file_index, file_offset)
  forall(i=1:ngas) psys%par(i)%file_index  = file_index(i)
  forall(i=1:ngas) psys%par(i)%file_offset = file_offset(i)
  deallocate(file_index, file_offset)


  ! read positions 
  !-----------------------------------------------------------!
  VarName = 'Coordinates'
  allocate(rblck3(3,ngas), stat=err)
  if(err/=0) call myerr("allocating rblck3 for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), rblck3)
  forall(i=1:ngas) psys%par(i)%pos(1) = rblck3(1,i)
  forall(i=1:ngas) psys%par(i)%pos(2) = rblck3(2,i)
  forall(i=1:ngas) psys%par(i)%pos(3) = rblck3(3,i)
  deallocate(rblck3)

  ! read velocities 
  !-----------------------------------------------------------!  
#ifdef incVel
  VarName = 'Velocity'
  allocate( rblck3(3,ngas), stat=err )
  if(err/=0) call myerr("allocating rblck3 for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), rblck3)
  forall(i=1:ngas) psys%par(i)%vel(1) = rblck3(1,i)
  forall(i=1:ngas) psys%par(i)%vel(2) = rblck3(2,i)
  forall(i=1:ngas) psys%par(i)%vel(3) = rblck3(3,i)
  deallocate(rblck3)
#endif

  ! read id's
  !-----------------------------------------------------------!  
  VarName = 'ParticleIDs'
  allocate( iblck(ngas), stat=err )
  if(err/=0) call myerr("allocating iblck for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), iblck)
  forall(i=1:ngas) psys%par(i)%id = iblck(i)
  deallocate(iblck)

  ! read masses
  !-----------------------------------------------------------!  
  VarName = 'Mass'
  allocate( rblck(ngas), stat=err )
  if(err/=0) call myerr("allocating rblck for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), rblck)
  forall(i=1:ngas) psys%par(i)%mass = rblck(i)
  deallocate(rblck)

  ! read temperature
  !-----------------------------------------------------------!  
  VarName = 'Temperature'
  allocate( rblck(ngas), stat=err )
  if(err/=0) call myerr("allocating rblck for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), rblck)
  forall(i=1:ngas) psys%par(i)%T = rblck(i)
  deallocate(rblck)

  ! read EOS
  !-----------------------------------------------------------!  
#ifdef incEOS
  VarName = trim(CV%EOSvarName)
  allocate( rblck(ngas), stat=err )
  if(err/=0) call myerr("allocating rblck for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), rblck)
  forall(i=1:ngas) psys%par(i)%eos = rblck(i)
  deallocate(rblck)
#endif

  ! read density 
  !-----------------------------------------------------------!  
  VarName = 'Density'
  allocate(rblck(ngas), stat=err)  
  if(err/=0) call myerr("allocating rblck for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), rblck)
  forall(i=1:ngas) psys%par(i)%rho = rblck(i)
  deallocate(rblck)

  ! read smoothing lengths 
  !-----------------------------------------------------------!  
  VarName = 'SmoothingLength'
  allocate(rblck(ngas), stat=err)
  if(err/=0) call myerr("allocating rblck for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), rblck)
  forall(i=1:ngas) psys%par(i)%hsml = rblck(i)
  deallocate(rblck)


  ! read Hydrogen mass fractions
  !-----------------------------------------------------------!  
#ifdef incHmf
  VarName = 'ElementAbundance/Hydrogen'
  allocate(rblck(ngas), stat=err)
  if(err/=0) call myerr("allocating rblck for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), rblck)
  forall(i=1:ngas) psys%par(i)%Hmf = rblck(i)
  deallocate(rblck)
#endif


  ! read Helium mass fractions
  !-----------------------------------------------------------!  
#ifdef incHemf
  VarName = 'ElementAbundance/Helium'
  allocate(rblck(ngas), stat=err)
  if(err/=0) call myerr("allocating rblck for "//trim(VarName), myname, crash)
  np = read_dataset(snap, 0, trim(VarName), rblck)
  forall(i=1:ngas) psys%par(i)%Hemf = rblck(i)
  deallocate(rblck)
#endif



  ! jhc finalization routines 
  !-------------------------------------------------------------
  call clear_selection(snap)
  call close_snapshot(snap)




  ! adjust ISM properties if we want
  !==================================================
#ifdef incEOS

  ! Temperature
  if (CV%AdjustIsmPars) then
     fmt = "(T7, A, T40, ES15.5)"
     write(str,fmt) "ISM particles have f_hot: ", CV%IsmHotFrac
     call mywrite(str,verb)
     write(str,fmt) "setting ISM particles to T_ISM: ", CV%IsmTemp
     call mywrite(str,verb)
     do ipar = 1, size(psys%par)         
        if ( psys%par(ipar)%eos > zero ) then
           psys%par(ipar)%T = CV%IsmTemp
        endif
     end do
  endif
     

  ! Molecular Hydrogen
#ifdef incH2

  fmt = "(T7, A, T30, L1)"
  write(str,fmt) "Add H2 to ISM?: ", CV%IsmAddH2
  call mywrite(str,verb)

  if (CV%IsmAddH2) then
     call mywrite("   calculating fH2 for EAGLE input", verb)
     
     fmt = "(T7, A, T30, ES15.5)"
     write(str,fmt) "EOS effective gamma: ", GV%sf_gamma_eos
     call mywrite(str,verb)
     
     write(str,fmt) "EOS nH threshold: ", GV%nH_star
     call mywrite(str,verb)
     
     eos_prs_norm = 1.0d3 * gconst%BOLTZMANN / GV%nH_star**GV%sf_gamma_eos
     
     do ipar = 1, size(psys%par) 

        if ( psys%par(ipar)%eos > zero ) then
           nH8 = return_nH( psys%par(ipar), shead%h, shead%a, adjust_ism=CV%AdjustIsmPars) 
           prs = eos_prs_norm * nH8**GV%sf_gamma_eos
           ratio = ( prs / gconst%BOLTZMANN ) / CV%IsmH2Param1     
           psys%par(ipar)%fH2 = one / ( one + ratio**(-CV%IsmH2Param2) )
        else
           !        mu = 4.0d0 / ( 3.0d0 * H_mf + 1.0d0 + 4.0d0 * H_mf * psys%par(ipar)%ye )
           !        prs = gconst%BOLTZMANN / ( mu * gconst%PROTONMASS ) * &
           !             psys%par(ipar)%rho * (h**2 / a**3) * GV%cgs_rho * &
           !             psys%par(ipar)%T
           prs = zero
           psys%par(ipar)%fH2 = zero
        endif

     end do
     
  else
  
     do ipar = 1, size(psys%par)
        psys%par(ipar)%fH2 = zero
     end do

  endif

  fmt = "(T7, A, 2ES15.5)"
  write(str,fmt) "min/max fH2 = ", minval( psys%par%fH2 ), maxval( psys%par%fH2 )
  call mywrite(str,verb) 
  call mywrite('',verb)

#endif
#endif



  ! calculate xHI from CLOUDY iontables (USE ANALYTIC SOLUTION)
  !-----------------------------------------------------------!  
  call mywrite('',verb) 
  call mywrite("   calculating input xHI from THIN ANALYTIC SOLUTION ~ CLOUDY tables", verb)

  caseA = .true.

  ! loop through the gas particles and interpolate from the table
  !---------------------------------------------------------------
  do i = 1,size(psys%par)
     nH8 = return_nH( psys%par(i), h=shead%h, a=shead%a, adjust_ism=CV%AdjustIsmPars )
     T8  = psys%par(i)%T
     psys%par(i)%xHI = analytic_xHIeq( T8, GV%UVB_gammaHI_cloudy, &
          nH8, y=CV%NeBackground, tauHI_eff=zero )

     ! set xHII from xHI 
     !-----------------------------------------------------------!  
     psys%par(i)%xHII = one - psys%par(i)%xHI

     
     ! store in dedicated variable
     !-----------------------------------------------------------!  
#ifdef incCloudy
     psys%par(i)%xHI_cloudy = psys%par(i)%xHI
#endif

  end do


#ifdef incCloudy
  fmt = "(T7, A, 2ES15.5)"
  write(str,fmt) "min/max xHI_cloudy = ", minval( psys%par%xHI_cloudy ), maxval( psys%par%xHI_cloudy )
  call mywrite(str,verb) 
#endif
  str = ''
  call mywrite(str,verb)

  




  ! if Helium, initialize ionization fractions to collisional 
  ! equilibrium set caseA true or false for collisional equilibrium
  !------------------------------------------------------------
#ifdef incHe
  call set_collisional_ionization_equilibrium( psys, caseA, CV%IsoTemp, DoHydrogen=.false., fit='hui' )
#endif


  ! set the electron fractions from the ionization fractions
  !----------------------------------------------------------
  call set_ye( psys )



  !do i =1, np, 1
  !   write(*,'(i20,5f14.3)') psys%par(i)%id, psys%par(i)%pos(1:3), psys%par(i)%rho, psys%par(i)%T 
  !end do

  !stop 'debug'

#endif

end subroutine read_Geagle_particles






!!$!> reads a Gadget EAGLE HDF5 snapshot into a particle array  
!!$!========================================================================
!!$subroutine read_Geagle_particles()
!!$
!!$  character(clen), parameter :: myname="read_Geagle_particles" 
!!$  logical, parameter :: crash=.true.
!!$  integer, parameter :: verb=2
!!$  character(clen) :: str,fmt
!!$  
!!$  real(r4b), allocatable :: rblck(:)
!!$  real(r8b), allocatable :: rblck3(:,:)
!!$  integer(i8b), allocatable :: iblck(:)
!!$  integer(i8b) :: ngasread
!!$
!!$  type(gadget_owls_header_type) :: ghead
!!$  type(gadget_sphray_header_type) :: shead
!!$  character(clen) :: snapfile, VarName, GroupName
!!$  integer(i4b) :: fh
!!$  integer(i8b) :: i
!!$  integer(i4b) :: err
!!$
!!$  integer(i8b) :: npar, ngas, nmass, ipar
!!$  integer(i8b) :: npar1, ngas1, nmass1
!!$  logical :: varmass(0:5)
!!$  integer(i4b) :: fn
!!$
!!$  real(r8b) :: meanweight
!!$  logical :: caseA(2)
!!$  real(r8b) :: MB
!!$
!!$
!!$  real(r8b) :: ratio
!!$  real(r8b) :: mu
!!$  real(r8b) :: prs
!!$  real(r8b) :: eos_prs_norm
!!$
!!$  real(r8b) :: Tdum
!!$  real(r8b) :: Hmf
!!$  real(r8b) :: nH8
!!$  real(r8b) :: T8
!!$  real(r8b) :: xvec(5)
!!$
!!$  type(gadget_constants_type) :: gconst
!!$
!!$  logical :: hdf5bool
!!$
!!$
!!$  ! set hdf5 boolean
!!$  !======================================================
!!$  hdf5bool = .true.
!!$
!!$  ! determine maximum number of gas particles in any file
!!$  ! twice this value will be the maximum that can be read
!!$  ! in if Urchin Selection is chosen in the config file
!!$  !-----------------------------------------------------
!!$  psys%ngas_max = 2*maxval(saved_gheads(:)%npar_file(0))
!!$
!!$  ! set local particle numbers
!!$  !============================
!!$  shead = saved_gheads( 0 )
!!$  varmass = (shead%npar_all > 0 .and. shead%mass == 0)
!!$  npar = sum(shead%npar_all)
!!$  ngas = psys%ngas_max
!!$  nmass = sum(shead%npar_all, mask=varmass)
!!$
!!$
!!$  ! do Gadget dummy checks
!!$  !============================
!!$  if (ngas == 0) call myerr("snapshot has no gas particles",myname,crash)
!!$
!!$  ! calculate bytes per particle and allocate particle array
!!$  !===========================================================
!!$  MB = GV%bytesperpar * real(ngas) / 2.0d0**20
!!$  GV%MB = GV%MB + MB
!!$
!!$  fmt="(A,F10.4,A,I10,A)"
!!$  write(str,fmt) "   allocating ", MB, " MB for ", ngas, " particles"
!!$  call mywrite(str,verb) 
!!$
!!$  allocate (psys%par(ngas), stat=err)
!!$  if (err /= 0) call myerr("failed to allocate par",myname,crash)
!!$
!!$
!!$  ! now read all snapshot files
!!$  !==============================          
!!$  ngasread = 0
!!$  GroupName = 'PartType0/'
!!$  files: do fn = 0, shead%nfiles-1
!!$
!!$
!!$     ! recall the header info
!!$     !-----------------------------------------------------------!  
!!$     shead   = saved_gheads( fn )
!!$     varmass = (shead%npar_file > 0 .and. shead%mass == 0)
!!$     npar1   = sum(shead%npar_file)
!!$     ngas1   = shead%npar_file(0)
!!$     nmass1  = sum(shead%npar_file, mask=varmass)
!!$     if (ngas1 == 0) cycle
!!$     if(ngasread + ngas1 > psys%ngas_max) then
!!$        write (*,*) ' sorry: need to increase ngas_max '
!!$        stop
!!$     endif
!!$
!!$
!!$     ! begin read
!!$     !-----------------------------------------------------------!  
!!$     call form_gadget_snapshot_file_name(CV%SnapPath,CV%ParFileBase,fn,snapfile,hdf5bool)
!!$     call mywrite("   reading particle snapshot file: "//trim(snapfile), verb)
!!$     call hdf5_open_file(fh, snapfile, readonly=.true.)
!!$     call gadget_owls_header_read_lun(ghead,fh)
!!$
!!$
!!$
!!$     ! read positions 
!!$     !-----------------------------------------------------------!
!!$     allocate(rblck3(3,ngas1), stat=err)
!!$     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
!!$     VarName = 'Coordinates'
!!$     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%pos(1) = rblck3(1,i)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%pos(2) = rblck3(2,i)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%pos(3) = rblck3(3,i)
!!$     deallocate(rblck3)
!!$     call gadget_data_attributes_read_lun( gattrs%pos, fh, GroupName, VarName )
!!$
!!$
!!$
!!$     ! read velocities 
!!$     !-----------------------------------------------------------!  
!!$#ifdef incVel
!!$     allocate(rblck3(3,ngas1), stat=err)
!!$     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
!!$     VarName = 'Velocity'
!!$     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%vel(1) = rblck3(1,i)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%vel(2) = rblck3(2,i)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%vel(3) = rblck3(3,i)
!!$     deallocate(rblck3)
!!$     call gadget_data_attributes_read_lun( gattrs%vel, fh, GroupName, VarName )
!!$#endif
!!$
!!$     ! read id's 
!!$     !-----------------------------------------------------------!  
!!$     allocate(iblck(ngas1), stat=err )
!!$     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
!!$     VarName = 'ParticleIDs'
!!$     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),iblck)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%id = iblck(i)
!!$     deallocate(iblck)
!!$     call gadget_data_attributes_read_lun( gattrs%id, fh, GroupName, VarName ) 
!!$
!!$
!!$     ! read masses 
!!$     !-----------------------------------------------------------!  
!!$
!!$     ! if gas particles are variable mass
!!$     if (varmass(0)) then  
!!$        allocate(rblck(ngas1), stat=err)
!!$        if(err/=0) call myerr("allocating rblck for mass",myname,crash)
!!$        VarName = 'Mass'
!!$        call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
!!$        forall(i=1:ngas1) psys%par(ngasread+i)%mass = rblck(i)
!!$        deallocate(rblck)
!!$        call gadget_data_attributes_read_lun( gattrs%mass, fh, GroupName, VarName ) 
!!$
!!$     ! if gas particles are isomass
!!$     else
!!$        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(0)       
!!$     end if
!!$
!!$     ! read temperature
!!$     !-----------------------------------------------------------!  
!!$     allocate(rblck(ngas1), stat=err)
!!$     if(err/=0) call myerr("allocating rblck for T",myname,crash)
!!$     VarName = 'Temperature'
!!$     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%T = rblck(i)
!!$     deallocate(rblck)
!!$     call gadget_data_attributes_read_lun( gattrs%T, fh, GroupName, VarName )
!!$
!!$     ! read EOS
!!$     !-----------------------------------------------------------!  
!!$#ifdef incEOS
!!$     VarName = trim(CV%EOSvarName)
!!$     allocate(rblck(ngas1), stat=err)
!!$     if(err/=0) call myerr("allocating rblck for "//trim(VarName),myname,crash)
!!$     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%eos = rblck(i)
!!$     deallocate(rblck)
!!$     call gadget_data_attributes_read_lun( gattrs%eos, fh, GroupName, VarName )
!!$#endif
!!$
!!$     ! read density 
!!$     !-----------------------------------------------------------!  
!!$     allocate(rblck(ngas1), stat=err)
!!$     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
!!$     VarName = 'Density'
!!$     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%rho = rblck(i)
!!$     deallocate(rblck)
!!$     call gadget_data_attributes_read_lun( gattrs%rho, fh, GroupName, VarName ) 
!!$
!!$     ! read smoothing lengths 
!!$     !-----------------------------------------------------------!  
!!$     allocate(rblck(ngas1), stat=err)
!!$     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
!!$     VarName = 'SmoothingLength'
!!$     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%hsml = rblck(i)
!!$     deallocate(rblck)
!!$     call gadget_data_attributes_read_lun( gattrs%hsml, fh, GroupName, VarName )
!!$
!!$     ! read Hydrogen mass fractions
!!$     !-----------------------------------------------------------!  
!!$#ifdef incHmf
!!$     allocate(rblck(ngas1), stat=err)
!!$     if(err/=0) call myerr("allocating rblck for Hmf",myname,crash)
!!$     VarName = 'ElementAbundance/Hydrogen'
!!$     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%Hmf = rblck(i)
!!$     deallocate(rblck)
!!$     call gadget_data_attributes_read_lun( gattrs%Hmf, fh, GroupName, VarName )
!!$#endif
!!$
!!$
!!$     ! read Helium mass fractions
!!$     !-----------------------------------------------------------!  
!!$#ifdef incHemf
!!$     allocate(rblck(ngas1), stat=err)
!!$     if(err/=0) call myerr("allocating rblck for Hemf",myname,crash)
!!$     VarName = 'ElementAbundance/Helium'
!!$     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
!!$     forall(i=1:ngas1) psys%par(ngasread+i)%Hemf = rblck(i)
!!$     deallocate(rblck)
!!$     call gadget_data_attributes_read_lun( gattrs%Hemf, fh, GroupName, VarName )
!!$#endif
!!$
!!$
!!$
!!$     ! constrain_psys_to_selection updates ngasread using the 
!!$     ! number of particles in the selction from this file
!!$     ! otherwise update using all gas particles in the file
!!$     !-----------------------------------------------------------!
!!$     if (CV%DoSelection) then
!!$        call constrain_psys_to_selection(psys,ngasread)
!!$     else
!!$        ngasread = ngasread + ngas1
!!$     endif
!!$
!!$     !-----------------------------------------------------------!
!!$     call hdf5_close_file(fh)
!!$
!!$  end do files
!!$
!!$  ! final call of constrain fixes size of psys
!!$  ! to be the number of particles in the selection
!!$  !--------------------------------------------------
!!$  if (CV%DoSelection) then
!!$     call constrain_psys_to_selection(psys)
!!$  endif
!!$
!!$
!!$
!!$
!!$
!!$  ! adjust ISM properties if we want
!!$  !==================================================
!!$#ifdef incEOS
!!$
!!$  ! Temperature
!!$  if (CV%AdjustIsmPars) then
!!$     fmt = "(T7, A, T40, ES15.5)"
!!$     write(str,fmt) "ISM particles have f_hot: ", CV%IsmHotFrac
!!$     call mywrite(str,verb)
!!$     write(str,fmt) "setting ISM particles to T_ISM: ", CV%IsmTemp
!!$     call mywrite(str,verb)
!!$     do ipar = 1, size(psys%par)         
!!$        if ( psys%par(ipar)%eos > zero ) then
!!$           psys%par(ipar)%T = CV%IsmTemp
!!$        endif
!!$     end do
!!$  endif
!!$     
!!$
!!$  ! Molecular Hydrogen
!!$#ifdef incH2
!!$
!!$  fmt = "(T7, A, T30, L1)"
!!$  write(str,fmt) "Add H2 to ISM?: ", CV%IsmAddH2
!!$  call mywrite(str,verb)
!!$
!!$  if (CV%IsmAddH2) then
!!$     call mywrite("   calculating fH2 for EAGLE input", verb)
!!$     
!!$     fmt = "(T7, A, T30, ES15.5)"
!!$     write(str,fmt) "EOS effective gamma: ", GV%sf_gamma_eos
!!$     call mywrite(str,verb)
!!$     
!!$     write(str,fmt) "EOS nH threshold: ", GV%nH_star
!!$     call mywrite(str,verb)
!!$     
!!$     eos_prs_norm = 1.0d3 * gconst%BOLTZMANN / GV%nH_star**GV%sf_gamma_eos
!!$     
!!$     do ipar = 1, size(psys%par) 
!!$
!!$        if ( psys%par(ipar)%eos > zero ) then
!!$           nH8 = return_nH( psys%par(ipar), shead%h, shead%a, adjust_ism=CV%AdjustIsmPars) 
!!$           prs = eos_prs_norm * nH8**GV%sf_gamma_eos
!!$           ratio = ( prs / gconst%BOLTZMANN ) / CV%IsmH2Param1     
!!$           psys%par(ipar)%fH2 = one / ( one + ratio**(-CV%IsmH2Param2) )
!!$        else
!!$           !        mu = 4.0d0 / ( 3.0d0 * H_mf + 1.0d0 + 4.0d0 * H_mf * psys%par(ipar)%ye )
!!$           !        prs = gconst%BOLTZMANN / ( mu * gconst%PROTONMASS ) * &
!!$           !             psys%par(ipar)%rho * (h**2 / a**3) * GV%cgs_rho * &
!!$           !             psys%par(ipar)%T
!!$           prs = zero
!!$           psys%par(ipar)%fH2 = zero
!!$        endif
!!$
!!$     end do
!!$     
!!$  else
!!$  
!!$     do ipar = 1, size(psys%par)
!!$        psys%par(ipar)%fH2 = zero
!!$     end do
!!$
!!$  endif
!!$
!!$  fmt = "(T7, A, 2ES15.5)"
!!$  write(str,fmt) "min/max fH2 = ", minval( psys%par%fH2 ), maxval( psys%par%fH2 )
!!$  call mywrite(str,verb) 
!!$  call mywrite('',verb)
!!$
!!$#endif
!!$#endif
!!$
!!$
!!$
!!$  ! calculate xHI from CLOUDY iontables (USE ANALYTIC SOLUTION)
!!$  !-----------------------------------------------------------!  
!!$  call mywrite('',verb) 
!!$  call mywrite("   calculating input xHI from THIN ANALYTIC SOLUTION ~ CLOUDY tables", verb)
!!$
!!$  caseA = .true.
!!$
!!$  ! loop through the gas particles and interpolate from the table
!!$  !---------------------------------------------------------------
!!$  do i = 1,size(psys%par)
!!$     nH8 = return_nH( psys%par(i), h=shead%h, a=shead%a, adjust_ism=CV%AdjustIsmPars )
!!$     T8  = psys%par(i)%T
!!$     psys%par(i)%xHI = analytic_xHIeq( T8, GV%UVB_gammaHI_cloudy, &
!!$          nH8, y=CV%NeBackground, tauHI_eff=zero )
!!$
!!$     ! set xHII from xHI 
!!$     !-----------------------------------------------------------!  
!!$     psys%par(i)%xHII = one - psys%par(i)%xHI
!!$
!!$     
!!$     ! store in dedicated variable
!!$     !-----------------------------------------------------------!  
!!$#ifdef incCloudy
!!$     psys%par(i)%xHI_cloudy = psys%par(i)%xHI
!!$#endif
!!$
!!$  end do
!!$
!!$
!!$#ifdef incCloudy
!!$  fmt = "(T7, A, 2ES15.5)"
!!$  write(str,fmt) "min/max xHI_cloudy = ", minval( psys%par%xHI_cloudy ), maxval( psys%par%xHI_cloudy )
!!$  call mywrite(str,verb) 
!!$#endif
!!$  str = ''
!!$  call mywrite(str,verb)
!!$
!!$  
!!$
!!$
!!$
!!$
!!$  ! if Helium, initialize ionization fractions to collisional 
!!$  ! equilibrium set caseA true or false for collisional equilibrium
!!$  !------------------------------------------------------------
!!$#ifdef incHe
!!$  call set_collisional_ionization_equilibrium( psys, caseA, CV%IsoTemp, DoHydrogen=.false., fit='hui' )
!!$#endif
!!$
!!$
!!$  ! set the electron fractions from the ionization fractions
!!$  !----------------------------------------------------------
!!$  call set_ye( psys )
!!$
!!$
!!$
!!$end subroutine read_Geagle_particles







end module gadget_eagle_input_mod
