!> \file global.F90

!> \brief the module that handles global variables
!<
module global_mod
use myf90_mod 
use gadget_general_class, only: gadget_constants_type
use gadget_general_class, only: gadget_owls_parameters_type
use gadget_sphray_header_class, only: gadget_sphray_header_type
use gadget_sphray_header_class, only: gadget_sphray_header_print_lun
use particle_system_mod, only: particle_system_type
use oct_tree_mod, only: oct_tree_type
use raylist_mod, only: raylist_type
use atomic_rates_mod, only: atomic_rates_table_type
use atomic_rates_mod, only: atomic_rates_type
use ion_table_class, only: ion_table_type

use config_mod, only: CV

implicit none



! global variables
!=====================
type(particle_system_type) :: psys        !< particles + box
type(raylist_type) :: globalraylist       !< ray/particle intersections (for viewer)
type(oct_tree_type) :: tree               !< octree

type(atomic_rates_table_type) :: rtable   !< rates read in from file
type(atomic_rates_type) :: isoT_k         !< static rates for iso-temperature run
type(atomic_rates_type) :: cmbT_k         !< static rates for cmb-temperature
type(atomic_rates_type) :: xHII_k         !< static rates for xHII-temperature 

type(ion_table_type) :: h1_itab           !< H1   Cloudy Ionization Balance Table
type(ion_table_type) :: he1_itab          !< He1  Cloudy Ionization Balance Table
type(ion_table_type) :: he2_itab          !< He2  Cloudy Ionization Balance Table

type(gadget_constants_type), save :: gconst     !< gadget constants
type(gadget_owls_parameters_type), save :: global_owls_parameters !< OWLS/GIMIC run parameters

character(clen), allocatable :: snap_file_names(:)               !< all snapshot file names (nfiles)
type(gadget_sphray_header_type), allocatable :: saved_gheads(:) !< all headers (nfiles)




 
!> global variables type. 
!=========================
type global_variables_type


   ! these are all set in get_planning_data() called in initialize.f90
   !-----------------------------------------------------------------
   real(r8b) :: cgs_len   !< code length [cm/h]
   real(r8b) :: cgs_mass  !< code mass [g/h]
   real(r8b) :: cgs_vel   !< code velocity [cm/s] 

   real(r8b) :: cgs_time  !< code time [s/h]
   real(r8b) :: cgs_rho   !< code density [g/cm^3 h^2]
   real(r8b) :: cgs_prs   !< code pressure [dyne/cm^2 h^2]
   real(r8b) :: cgs_enrg  !< code energy [ergs/h]
   
   real(r8b) :: cgs_ilen2 !< 1 / (code length)^2 [h^2/cm^2]


   real(r8b) :: OmegaM    !< matter / critical density z=0
   real(r8b) :: OmegaL    !< lambda / critical density z=0
   real(r8b) :: OmegaB    !< baryon / critical density z=0
   real(r8b) :: TimeGyr   !< time at snapshot in Gyr

   real(r8b) :: LittleH   !< Hubble parameter @ z=0 in units of 100 km/s/Mpc
   real(r8b) :: ScaleFac  !< Scale factor at snapshot
   real(r8b) :: RhoConv   !< LittleH^2 / ScaleFac^3
   real(r8b) :: RedShift  !< RedShift at snapshot
   real(r8b) :: Tcmb_cur  !< CMB temperature for the current snapshot

   real(r8b) :: BoxLwrsComoh(3)   !< comoving h-mangled input coords of lower x,y,z corner
   real(r8b) :: BoxUprsComoh(3)   !< comoving h-mangled input coords of upper x,y,z corner

   integer(i4b) :: OutputIndx     !< labels the outputs


   ! RAM tracking
   !-----------------------------------------------------------
   integer(i4b) :: bytesperpar  !< bytes per particle
   real(r8b) :: MB              !< tracks memory consumption

   real(r8b) :: total_mass      !< summed mass of all particles in a snapshot
   real(r8b) :: total_H_atoms   !< sum of all H  atoms in computational volume
   real(r8b) :: total_He_atoms  !< sum of all He atoms in computational volume

   real(r8b) :: sf_gamma_eos       !< index for polytropic equation of state for star forming gas.
   real(r8b) :: nH_star            !< star formation nH threshold
   real(r8b) :: UVB_gammaHI_cloudy !< magnitude of UVB from cloudy ionization table

   
   ! these are updated continuosly while the code runs 
   ! and most are initialized in initialize.f90
   !----------------------------------------------------
   character(clen) :: raystat_file !< file where ray stats are put
   integer(i4b) :: raystatlun      !< lun for ray stat log file

   character(clen) :: pardata_file !< file with particle data summaries
   integer(i4b) :: pardatalun      !< lun for particle data log file
   
   integer(i8b) :: rayn            !< current ray number 

   real(r8b) :: nwionfrac          !< number weighted ionization fraction
   real(r8b) :: mwionfrac          !< mass weighted ionization fraction
   real(r8b) :: vwionfrac          !< volume weighted ionization fraction
   
end type global_variables_type
 

type(global_variables_type) :: GV           !< global variables          


contains



subroutine write_units_to_log_file()
  character(clen) :: logfile
  integer(i4b) :: loglun
  real(r8b) :: kpc2cm
  real(r8b) :: km2cm
  character(clen) :: fmt

  ! write units to log file
  !===================================================

  logfile = trim(CV%OutputDir) // "/" // "code_units.log"
  call open_formatted_file_w(logfile,loglun)

  kpc2cm = gconst%cm_per_mpc * 1.0d-3
  km2cm = 1.0d5

  fmt = '(A,ES12.5,A)'
  write(loglun,*) 
  write(loglun,'(A)') "setting code units ..."
  write(loglun,fmt) "  Hubble:     ", GV%LittleH, " = H0[km/s/Mpc] / 100 "
  write(loglun,fmt) "  length:     ", GV%cgs_len / kpc2cm, " [kpc/h]"
  write(loglun,fmt) "  mass:       ", GV%cgs_mass / gconst%solar_mass, " [Msun/h]"
  write(loglun,fmt) "  velocity    ", GV%cgs_vel / km2cm, " [km/s]"
  write(loglun,fmt) "  time:       ", GV%cgs_time / gconst%sec_per_megayear, " [Myr/h]"
  write(loglun,fmt) "  density:    ", GV%cgs_rho / (gconst%solar_mass / kpc2cm**3), " [Msun/kpc^3 h^2]"
  write(loglun,fmt) "  pressure:   ", GV%cgs_prs, " [dyne/cm^2 h^2]"
  write(loglun,fmt) "  energy:     ", GV%cgs_enrg, " [ergs/h]"
  write(loglun,*) 

  close(loglun)

end subroutine write_units_to_log_file





subroutine write_particle_headers_to_log_file()
  character(clen) :: logfile
  integer(i4b) :: loglun
  integer(i4b) :: ifile

  ! write all particle headers to log file
  !===================================================
  logfile = trim(CV%OutputDir) // "/" // "particle_headers.log"
  call open_formatted_file_w(logfile, loglun)

  write(loglun,'(A)') "logging GADGET particle header(s) ... "
  do ifile = 0, CV%ParFilesPerSnap-1
  
     write(loglun,'("     ",A)') trim(snap_file_names(ifile))    
     call gadget_sphray_header_print_lun(saved_gheads(ifile), loglun)
               
  end do

  close(loglun)


end subroutine write_particle_headers_to_log_file



!> puts read in header into a header for output
!---------------------------------------------------------
 
subroutine recall_saved_ghead(fnum, ghead)

  integer, intent(in) :: fnum                    !< which snapshot file
  type(gadget_sphray_header_type), intent(out) :: ghead !< gadget header
  
  
  ghead%npar_all = 0
  ghead%npar_all(0) = saved_gheads(fnum)%npar_all(0)     
  
  ghead%npar_file = 0
  ghead%npar_file(0) = saved_gheads(fnum)%npar_file(0)     
  
  ghead%mass = 0.0
  ghead%nfiles = saved_gheads(fnum)%nfiles
  
  ghead%a = saved_gheads(fnum)%a
  ghead%z = saved_gheads(fnum)%z

  ghead%boxlen   = saved_gheads(fnum)%boxlen
  ghead%OmegaM   = saved_gheads(fnum)%OmegaM
  ghead%OmegaL   = saved_gheads(fnum)%OmegaL
  ghead%OmegaB   = saved_gheads(fnum)%OmegaB
  ghead%time_gyr = saved_gheads(fnum)%time_gyr
  ghead%h        = saved_gheads(fnum)%h
  
  ghead%npar_hw = 0
  ghead%npar_hw(0) = saved_gheads(fnum)%npar_hw(0)
  
  ghead%flag_sfr      = saved_gheads(fnum)%flag_sfr
  ghead%flag_feedback = saved_gheads(fnum)%flag_feedback
  ghead%flag_cooling  = saved_gheads(fnum)%flag_cooling
  ghead%flag_age      = saved_gheads(fnum)%flag_age
  ghead%flag_metals   = saved_gheads(fnum)%flag_metals
  ghead%flag_entr_ics = saved_gheads(fnum)%flag_entr_ics
  
  ghead%rays_traced = GV%rayn
  
#ifdef incHmf
  ghead%flag_Hmf = 1
#else
  ghead%flag_Hmf = 0
#endif
  !-------------------------
#ifdef incHemf
  ghead%flag_Hemf = 1
#else
  ghead%flag_Hemf = 0
#endif
  !-------------------------
#ifdef incHe
  ghead%flag_helium = 1
#else
  ghead%flag_helium = 0
#endif
  !-------------------------
#ifdef outGammaHI
  ghead%flag_gammaHI = 1
#else
  ghead%flag_gammaHI = 0
#endif
  !-------------------------
#ifdef incCloudy
  ghead%flag_cloudy = 1
#else
  ghead%flag_cloudy = 0
#endif
  !-------------------------
#ifdef incEOS
  ghead%flag_eos = 1
#else
  ghead%flag_eos = 0
#endif
  !-------------------------
#ifdef incSFR
  ghead%flag_incsfr = 1
#else
  ghead%flag_incsfr = 0
#endif
  !-------------------------
  
  ghead%unused = 0

  
end subroutine recall_saved_ghead




end module global_mod
