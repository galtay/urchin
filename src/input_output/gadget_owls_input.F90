!> \file gadget_owls_input.F90

!> \brief Handles readin of GADGET OWLS/GIMIC HDF5 formatted files
!<

module gadget_owls_input_mod

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

public :: get_planning_data_gadget_owls
public :: read_Gowls_particles


contains




!>   gets run planning data from Gadget OWLS/GIMIC Headers
!============================================================
subroutine get_planning_data_gadget_owls()

  character(clen), parameter :: myname = 'get_planning_data_gadget_owls'
  logical, parameter :: crash = .true.
  integer, parameter :: verb = 2

  type(gadget_owls_header_type) :: ghead
  type(gadget_units_type) :: gunits
  type(gadget_constants_type) :: gconst
  type(gadget_owls_parameters_type) :: gparam

  integer(i4b) :: iSnap           ! initial and final snapshot numbers
  integer(i4b) :: pfiles          ! files/snap for particles    
  integer(i4b) :: i,j             ! counters
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


end subroutine get_planning_data_gadget_owls







!> reads a Gadget OWLS/GIMIC HDF5 snapshot into a particle array  
!========================================================================
subroutine read_Gowls_particles()

  character(clen), parameter :: myname="read_Gowls_particles" 
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt
  
  real(r4b), allocatable :: rblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_owls_header_type) :: ghead
  type(gadget_sphray_header_type) :: shead
  character(clen) :: snapfile, VarName, GroupName
  integer(i4b) :: fh
  integer(i8b) :: i
  integer(i4b) :: err

  integer(i8b) :: npar, ngas, nmass, ipar
  integer(i8b) :: npar1, ngas1, nmass1
  logical :: varmass(0:5)
  integer(i4b) :: fn

  real(r8b) :: meanweight
  logical :: caseA(2)
  real(r8b) :: MB

  real(r8b) :: ratio
  real(r8b) :: mu
  real(r8b) :: prs
  real(r8b) :: eos_prs_norm

  real(r8b) :: Tdum
  real(r8b) :: Hmf
  real(r8b) :: nH8
  real(r8b) :: T8
  real(r8b) :: xvec(5)

  type(gadget_constants_type) :: gconst

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
     call gadget_owls_header_read_lun(ghead,fh)



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
     call gadget_data_attributes_read_lun( gattrs%pos, fh, GroupName, VarName )



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
     call gadget_data_attributes_read_lun( gattrs%vel, fh, GroupName, VarName )
#endif

     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
     VarName = 'ParticleIDs'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),iblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%id = iblck(i)
     deallocate(iblck)
     call gadget_data_attributes_read_lun( gattrs%id, fh, GroupName, VarName ) 


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
        call gadget_data_attributes_read_lun( gattrs%mass, fh, GroupName, VarName ) 

     ! if gas particles are isomass
     else
        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(0)       
     end if

     ! read temperature
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for T",myname,crash)
     VarName = 'Temperature'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%T = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( gattrs%T, fh, GroupName, VarName )

     ! read EOS
     !-----------------------------------------------------------!  
#ifdef incEOS
     VarName = trim(CV%EOSvarName)
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for "//trim(VarName),myname,crash)
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%eos = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( gattrs%eos, fh, GroupName, VarName )
#endif

     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
     VarName = 'Density'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%rho = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( gattrs%rho, fh, GroupName, VarName ) 

     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
     VarName = 'SmoothingLength'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%hsml = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( gattrs%hsml, fh, GroupName, VarName )

     ! read Hydrogen mass fractions
     !-----------------------------------------------------------!  
#ifdef incHmf
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for Hmf",myname,crash)
     VarName = 'ElementAbundance/Hydrogen'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%Hmf = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( gattrs%Hmf, fh, GroupName, VarName )
#endif


     ! read Helium mass fractions
     !-----------------------------------------------------------!  
#ifdef incHemf
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for Hemf",myname,crash)
     VarName = 'ElementAbundance/Helium'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%Hemf = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( gattrs%Hemf, fh, GroupName, VarName )
#endif


     !-----------------------------------------------------------!
     ngasread = ngasread + ngas1
     !-----------------------------------------------------------!

     call hdf5_close_file(fh)



  end do files






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
     call mywrite("   calculating fH2 for OWLS input", verb)
     
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



end subroutine read_Gowls_particles







end module gadget_owls_input_mod
