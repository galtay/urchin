!> \file gadget_general_class.f90

!> \brief Handles universal GADGET stuff.  
!!
!! Provides types to handle units and constants in a 
!! standard way.  Note that where 6 element arrays are needed that
!! correspond to particle types, we have indexed them from 0-5 as 
!! opposed to the default Fortran 1-6. 
!< 

module gadget_general_class
use myf90_mod
use hdf5_wrapper

#ifdef useMPI
  use mpi
#endif





implicit none
private


public :: gadget_ptype_names

public :: gadget_units_type
public :: gadget_units_set
public :: gadget_units_print_lun
public :: gadget_units_read_lun
public :: gadget_units_read_file
public :: gadget_units_write_lun
#ifdef useMPI
public :: gadget_units_broadcast
#endif

public :: gadget_constants_type
public :: gadget_constants_read_lun
public :: gadget_constants_read_file
public :: gadget_constants_write_lun

public :: gadget_one_data_attributes_type
public :: gadget_data_attributes_type
public :: gadget_data_attributes_set
public :: gadget_data_attributes_read_lun
public :: gadget_data_attributes_write_lun
public :: gattrs

public :: gadget_owls_parameters_type
public :: gadget_owls_chemical_elements_type
public :: gadget_owls_numerical_parameters_type
public :: gadget_owls_stellar_evolution_parameters_type
public :: gadget_owls_wind_parameters_type

public :: gadget_owls_parameters_read_lun
public :: gadget_owls_parameters_read_file
public :: gadget_owls_parameters_write_lun

public :: form_gadget_snapshot_file_name




!> Particle type names
!-----------------------------------------
character(5), parameter :: gadget_ptype_names(0:5) = &
     (/"gas  ","halo ","disk ","bulge","stars","bndry"/)


!> Units type (default= kpc/h, 1.0e10 Msolar/h, km/s)
!-----------------------------------------------------
type gadget_units_type
   real(r8b) :: cgs_length   = 3.0856780d21  !<  [cm h^-1]                
   real(r8b) :: cgs_mass     = 1.989d43      !<  [g h^-1]                 
   real(r8b) :: cgs_velocity = 1.0d5         !<  [cm s^-1]                
   real(r8b) :: cgs_time     = 3.085678d16   !<  [s h^-1]                 
   real(r8b) :: cgs_density  = 6.769911d-22  !<  [g cm^-3 h^2]            
   real(r8b) :: cgs_pressure = 6.769911d-12  !<  [ba = g cm^-1 s^-2 h^2]  
   real(r8b) :: cgs_energy   = 1.989d53      !<  [erg = g cm^2 s^-2 h^-1] 
end type gadget_units_type


!> Constants type ( 23 doubles )
!-----------------------------------------
type gadget_constants_type
   real(r8b) :: PI = 3.1415927d0                 !< [pure]                  
   real(r8b) :: GAMMA = 1.6666667d0              !< [pure]
   real(r8b) :: GRAVITY = 6.6720000d-08          !< [cm^3 g^-1s^-2]    
   real(r8b) :: SOLAR_MASS = 1.9890000d+33       !< [g]                        
   real(r8b) :: SOLAR_LUM = 3.8260000d+33        !< [erg/s]
   real(r8b) :: RAD_CONST = 7.5650000d-15        !< [erg cm^-3 K^-4]  
   real(r8b) :: AVOGADRO = 6.0222000d+23         !< [pure] 
   real(r8b) :: BOLTZMANN = 1.3806000d-16        !< [cm^2 g s^-2 K^-1] 
   real(r8b) :: GAS_CONST = 83142500.d0          !< [erg K^-1 mol^-1]  
   real(r8b) :: C = 2.9979000d+10                !< [cm/s]   
   real(r8b) :: PLANCK = 6.6262000d-27           !< [cm^2 g s^-1]  
   real(r8b) :: CM_PER_MPC = 3.0856780d+24       !< [pure]    
   real(r8b) :: PROTONMASS = 1.6726000d-24       !< [g]
   real(r8b) :: ELECTRONMASS = 9.1095300d-28     !< [g] 
   real(r8b) :: ELECTRONCHARGE = 4.8032000d-10   !< [esu]
   real(r8b) :: HUBBLE = 3.2407789d-18           !< [s^-1 h]
   real(r8b) :: T_CMB0 = 2.7280000d0             !< [K]
   real(r8b) :: SEC_PER_MEGAYEAR = 3.1550000d+13 !< [pure]
   real(r8b) :: SEC_PER_YEAR = 31550000.d0       !< [pure]
   real(r8b) :: STEFAN = 7.5657000d-15     !< = rad. const. [erg cm^-3 K^-4]  
   real(r8b) :: THOMPSON = 6.6524587d-25   !< [cm^2]
   real(r8b) :: EV_TO_ERG = 1.6021765d-12  !< [pure]
   real(r8b) :: Z_SOLAR = 0.012663729d0    !< [Mass Fraction] 
end type gadget_constants_type


!> OWLS / GIMIC Parameters
!-------------------------------------
type gadget_owls_chemical_elements_type
   real(i4b) :: BG_NELEMENTS = 9
   character(20) :: ElementNames(9) = (/"Hydrogen            ",&
                                        "Helium              ",&
                                        "Carbon              ",&
                                        "Nitrogen            ",&
                                        "Oxygen              ",&
                                        "Neon                ",&
                                        "Magnesium           ",&
                                        "Silicon             ",&
                                        "Iron                "/)
   real(r8b) :: InitAbundance_Hydrogen = 0.752  
   real(r8b) :: InitAbundance_Helium = 0.248    
   real(r8b) :: InitAbundance_Carbon = 0.0      
   real(r8b) :: InitAbundance_Nitrogen = 0.0    
   real(r8b) :: InitAbundance_Oxygen = 0.0      
   real(r8b) :: InitAbundance_Neon = 0.0        
   real(r8b) :: InitAbundance_Magnesium = 0.0   
   real(r8b) :: InitAbundance_Silicon = 0.0     
   real(r8b) :: InitAbundance_Iron = 0.0        
   real(r8b) :: CalciumOverSilicon = 0.0941736  
   real(r8b) :: SulphurOverSilicon = 0.605416   
end type gadget_owls_chemical_elements_type



type gadget_owls_numerical_parameters_type
   integer(i4b) :: ComovingIntegrationOn = 1            
   integer(i4b) :: TypeOfTimestepCriterion = 0          
   real(r8b)    :: ErrTolIntAccuracy = 0.025            
   real(r8b)    :: MaxSizeTimestep = 0.025              
   real(r8b)    :: MinSizeTimestep = 1.0d-10            
   real(r8b)    :: ErrTolTheta = 0.0                    
   integer(i4b) :: TypeOfOpeningCriterion = 1           
   real(r8b)    :: ErrTolForceAcc = 0.0050              
   real(r8b)    :: TreeDomainUpdateFrequency = 0.0050   
   integer(i4b) :: DesNumNgb = 48                       
   real(r8b)    :: MaxNumNgbDeviation = 1.0             
   real(r8b)    :: ArtBulkViscConst = 1.0               
   real(r8b)    :: InitGasU_ERG = 2.70847d10            
   real(r8b)    :: MinGasU_ERG = 5.0d8                  
   real(r8b)    :: CourantFac = 0.15                    
   real(r8b)    :: PartAllocFactor = 1.3                
   real(r8b)    :: TreeAllocFactor = 0.9257499999999999 
   real(r8b)    :: BufferSize = 1.060998251294d-300
   real(r8b)    :: MinGasHsmlFractional = 0.01          
end type gadget_owls_numerical_parameters_type



type gadget_owls_stellar_evolution_parameters_type
   real(r8b) :: SF_EOSGammaEffective = 1.33333333
   real(r8b) :: SF_EOSEnergyAtThreshold_ERG = 1.008d12
   real(r8b) :: SF_THRESH_MinPhysDens_HpCM3 = 0.1
   real(r8b) :: SF_THRESH_MinOverDens = 57.7
   real(r8b) :: SF_THRESH_MetDepExponent = -0.64
   real(r8b) :: SF_THRESH_MetDepFiducialZ = 0.0020
   real(r8b) :: SF_THRESH_MetDepMaxPhysDens_HpCM3 = 10.0
   real(r8b) :: SF_THRESH_MaxTemp_K = 57.7
   real(r8b) :: SF_SchmidtLawCoeff_MSUNpYRpKPC2 = 1.515d-4
   real(r8b) :: SF_SchmidtLawExponent = 1.4
   real(r8b) :: SF_SchmidtLawCoeff_GpSpCM2 = 1.423770005086782d-16
   integer(i4b) :: SF_MetDepDensThresholdOn = 0
   character(100) :: IMF_Model =         "Chabrier                                         " // &
                                         "                                                 "
   character(100) :: IMF_LifetimeModel = "P98                                              " // &
                                         "                                                 "
   real(r8b) :: IMF_MinMass_MSUN = 0.1
   real(r8b) :: IMF_MaxMass_MSUN = 100.0
   character(100) :: SNIa_Model =        "Efolding                                         " // &
                                         "                                                 "
   real(r8b) :: SNIa_Efficiency_fracwd = 0.01
   real(r8b) :: SNIa_Energy_ERG = 1.0d51
   integer(i4b) :: SNIa_MassTransferOn = 1
   integer(i4b) :: SNIa_EnergyTransferOn = 1
   integer(i4b) :: SNII_MassTransferOn = 1
   real(r8b) :: SNII_MinMass_MSUN = 6.0
   real(r8b) :: SNII_MaxMass_MSUN = 100.0
   real(r8b) :: SNII_Factor_Hydrogen = 1.0
   real(r8b) :: SNII_Factor_Helium = 1.0
   real(r8b) :: SNII_Factor_Carbon = 0.5
   real(r8b) :: SNII_Factor_Nitrogen = 1.0
   real(r8b) :: SNII_Factor_Oxygen = 1.0
   real(r8b) :: SNII_Factor_Neon = 1.0
   real(r8b) :: SNII_Factor_Magnesium = 2.0
   real(r8b) :: SNII_Factor_Silicon = 1.0
   real(r8b) :: SNII_Factor_Iron = 0.5
   integer(i4b) :: AGB_MassTransferOn = 1
   integer(i4b) :: POPIII_MassTransferOn = 0
   integer(i4b) :: POPIII_EnergyTransferOn = 0
   real(r8b) :: POPIII_Energy_ERG = 0.0
   real(r8b) :: POPIII_NumPerMsun = 0.0
end type gadget_owls_stellar_evolution_parameters_type



type gadget_owls_wind_parameters_type
   integer(i4b) :: SNII_WindOn = 1
   integer(i4b) :: SNII_WindIsotropicOn = 1
   real(r8b) :: SNII_WindSpeed_KMpS = 600.0
   real(r8b) :: SNII_WindMassLoading = 2.0
   real(r8b) :: SNII_WindDelay_YR = 3.0d7
end type gadget_owls_wind_parameters_type


type gadget_owls_parameters_type
   type(gadget_owls_chemical_elements_type) :: chemical_elements
   type(gadget_owls_numerical_parameters_type) :: numerical_parameters
   type(gadget_owls_stellar_evolution_parameters_type) :: stellar_evolution_parameters
   type(gadget_owls_wind_parameters_type) :: wind_parameters
end type gadget_owls_parameters_type


!> Attributes for data in HDF5 files.
!-----------------------------------------------------
type gadget_one_data_attributes_type
   real(r8b) :: CGSConversionFactor 
   real(r4b) :: h_scale_exponent
   real(r4b) :: aexp_scale_exponent
   character(200) :: VarDescription
end type gadget_one_data_attributes_type


type gadget_data_attributes_type

   type(gadget_one_data_attributes_type) :: pos
   type(gadget_one_data_attributes_type) :: id
   type(gadget_one_data_attributes_type) :: mass
   type(gadget_one_data_attributes_type) :: T
   type(gadget_one_data_attributes_type) :: rho
   type(gadget_one_data_attributes_type) :: ye
   type(gadget_one_data_attributes_type) :: xHI
   type(gadget_one_data_attributes_type) :: hsml
   type(gadget_one_data_attributes_type) :: lasthit
   
   type(gadget_one_data_attributes_type) :: vel
   type(gadget_one_data_attributes_type) :: xHI_cloudy
   type(gadget_one_data_attributes_type) :: Hmf
   type(gadget_one_data_attributes_type) :: xHeI
   type(gadget_one_data_attributes_type) :: xHeII
   type(gadget_one_data_attributes_type) :: Hemf
   type(gadget_one_data_attributes_type) :: gammaHI
   type(gadget_one_data_attributes_type) :: eos
   type(gadget_one_data_attributes_type) :: sfr
   type(gadget_one_data_attributes_type) :: fh2
   type(gadget_one_data_attributes_type) :: nH

end type gadget_data_attributes_type


type(gadget_data_attributes_type) :: gattrs


contains


!!============================================
!!
!!    UNITS    
!!
!!============================================



!> sets user supplied units
!--------------------------------------------------------------
subroutine gadget_units_set(this, cgs_length, cgs_mass, cgs_velocity)
  type(gadget_units_type) :: this
  real(r8b) :: cgs_length
  real(r8b) :: cgs_mass
  real(r8b) :: cgs_velocity

  this%cgs_length   = cgs_length
  this%cgs_mass     = cgs_mass
  this%cgs_velocity = cgs_velocity

  this%cgs_density  = this%cgs_mass / this%cgs_length**3
  this%cgs_energy   = this%cgs_mass * this%cgs_velocity**2
  this%cgs_time     = this%cgs_length / this%cgs_velocity
  this%cgs_pressure = this%cgs_mass / &
       (this%cgs_length**3 / this%cgs_velocity**2)

end subroutine gadget_units_set


!> formatted print of units to lun (including standard out)
!---------------------------------------------------------------
subroutine gadget_units_print_lun(this, lun, h)
  type(gadget_units_type), intent(in) :: this 
  integer(i4b), intent(in) :: lun
  real(r8b), intent(in) :: h

  type(gadget_constants_type) :: gconst

  character(clen) :: n1
  character(clen) :: star_fmt
  character(clen) :: unit_fmt
  character(clen) :: line_fmt

  real(r8b), parameter :: cm_per_km = 1.0d5

  ! binding energy of 1.0d8 [Msolar/h] halo @ z=0 
  ! Delta_c = 18 pi^2, units = erg/h
  ! http://arxiv.org/pdf/astro-ph/0010468v3
  real(r8b), parameter :: E8 = 5.45d53 * 1.0d-1

  real(r8b) :: rho_crit_0

  rho_crit_0 = 3.0d0 * gconst%hubble**2 / (8.0d0 * gconst%pi * gconst%gravity)

  star_fmt = "(78('='))"
  line_fmt = "(78('-'))"
  unit_fmt = "(A,T30,A)"

  write(lun,*)
  write(lun,star_fmt)
  write(lun,*) "   Units"
  write(lun,line_fmt)
  write(lun,*) "  h = ", h, " = H0[km/s/Mpc] / 100 "
  write(n1,'(ES20.6)') this%cgs_length
  write(lun,unit_fmt) "length [cm/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_mass
  write(lun,unit_fmt) "mass [g/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_velocity
  write(lun,unit_fmt) "velocity [cm/s]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_density
  write(lun,unit_fmt) "density [g/cm^3 h^2]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_energy
  write(lun,unit_fmt) "energy [erg/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_time
  write(lun,unit_fmt) "time [s/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_pressure
  write(lun,unit_fmt) "pressure [ba h^2]:", trim(adjustl(n1)) 

  write(lun,*) 

  write(n1,'(ES20.6)') rho_crit_0
  write(lun,unit_fmt) "rho_crit_0 [g/cm^3 h^2]:", trim(adjustl(n1))
  write(n1,'(ES20.6)') E8
  write(lun,unit_fmt) "E8 [erg/h]:", trim(adjustl(n1))

  write(lun,*)

  write(n1,'(ES20.6)') this%cgs_length / gconst%cm_per_mpc
  write(lun,unit_fmt) "length [Mpc/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_mass / gconst%solar_mass
  write(lun,unit_fmt) "mass [Msolar/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_velocity / cm_per_km
  write(lun,unit_fmt) "velocity [km/s]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_density / rho_crit_0
  write(lun,unit_fmt) "density [rho_crit_0 h^2]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_energy / E8
  write(lun,unit_fmt) "energy [E8]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_time / gconst%sec_per_megayear
  write(lun,unit_fmt) "time [Myr/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_pressure / gconst%boltzmann 
  write(lun,unit_fmt) "pressure/k_b [cm^-3 K h^2]:", trim(adjustl(n1)) 

  write(lun,star_fmt)

end subroutine gadget_units_print_lun


!> reads units from an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_units_read_lun(this, fh)
  type(gadget_units_type), intent(out) :: this
  integer, intent(in) :: fh
  
  call hdf5_read_attribute(fh, 'Units/UnitLength_in_cm', this%cgs_length)
  call hdf5_read_attribute(fh, 'Units/UnitMass_in_g', this%cgs_mass)
  call hdf5_read_attribute(fh, 'Units/UnitVelocity_in_cm_per_s', this%cgs_velocity)
  call hdf5_read_attribute(fh, 'Units/UnitDensity_in_cgs', this%cgs_density)
  call hdf5_read_attribute(fh, 'Units/UnitEnergy_in_cgs', this%cgs_energy)
  call hdf5_read_attribute(fh, 'Units/UnitPressure_in_cgs', this%cgs_pressure)
  call hdf5_read_attribute(fh, 'Units/UnitTime_in_s', this%cgs_time)
  
end subroutine gadget_units_read_lun


!> reads units from an hdf5 file
!--------------------------------------------------------------
subroutine gadget_units_read_file(this, snapfile)
  type(gadget_units_type), intent(out) :: this
  character(*), intent(in) :: snapfile
  integer :: fh
  
  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call gadget_units_read_lun(this, fh)
  call hdf5_close_file( fh )

end subroutine gadget_units_read_file


!> writes units to an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_units_write_lun(this, fh)
  type(gadget_units_type), intent(in) :: this
  integer, intent(in) :: fh

  call hdf5_write_attribute(fh, 'Units/UnitLength_in_cm', this%cgs_length)
  call hdf5_write_attribute(fh, 'Units/UnitMass_in_g', this%cgs_mass)
  call hdf5_write_attribute(fh, 'Units/UnitVelocity_in_cm_per_s', this%cgs_velocity)
  call hdf5_write_attribute(fh, 'Units/UnitDensity_in_cgs', this%cgs_density)
  call hdf5_write_attribute(fh, 'Units/UnitEnergy_in_cgs', this%cgs_energy)
  call hdf5_write_attribute(fh, 'Units/UnitPressure_in_cgs', this%cgs_pressure)
  call hdf5_write_attribute(fh, 'Units/UnitTime_in_s', this%cgs_time)


end subroutine gadget_units_write_lun


!> broadcasts units
!-----------------------------------------
#ifdef useMPI
subroutine gadget_units_broadcast( units )
  type(gadget_units_type) :: units
  integer :: count
  integer :: root
  integer :: ierr

  count = 7
  root = 0
  call mpi_bcast( units, count, mpi_double_precision, root, mpi_comm_world, ierr )
  if (ierr /= 0) then
     call mpi_barrier(mpi_comm_world, ierr)
     call mpi_finalize(ierr)
  endif

!  write(*,*) 'useMPI macro not defined in Makefile'
!  stop

end subroutine gadget_units_broadcast
#endif


!!============================================
!!
!!    OWLS Parameters
!!
!!============================================

!! applying the principle of read the least amount of stuff


!> reads OWLS parameters from an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_owls_parameters_read_lun(this, fh)
  type(gadget_owls_parameters_type), intent(out) :: this
  integer, intent(in) :: fh

  type(gadget_owls_chemical_elements_type) :: ce
  type(gadget_owls_numerical_parameters_type) :: np
  type(gadget_owls_stellar_evolution_parameters_type) :: sep
  type(gadget_owls_wind_parameters_type) :: wp

  character(clen) :: grp
  
  grp = 'Parameters/ChemicalElements/'
  call hdf5_read_attribute(fh,trim(grp)//'BG_NELEMENTS', ce%BG_NELEMENTS )  
  call hdf5_read_attribute(fh,trim(grp)//'ElementNames', ce%ElementNames )
  call hdf5_read_attribute(fh,trim(grp)//'InitAbundance_Hydrogen', ce%InitAbundance_Hydrogen )  
  call hdf5_read_attribute(fh,trim(grp)//'InitAbundance_Helium', ce%InitAbundance_Helium )
  call hdf5_read_attribute(fh,trim(grp)//'InitAbundance_Carbon', ce%InitAbundance_Carbon )      
  call hdf5_read_attribute(fh,trim(grp)//'InitAbundance_Nitrogen', ce%InitAbundance_Nitrogen )
  call hdf5_read_attribute(fh,trim(grp)//'InitAbundance_Oxygen', ce%InitAbundance_Oxygen )      
  call hdf5_read_attribute(fh,trim(grp)//'InitAbundance_Neon', ce%InitAbundance_Neon )        
  call hdf5_read_attribute(fh,trim(grp)//'InitAbundance_Magnesium', ce%InitAbundance_Magnesium )   
  call hdf5_read_attribute(fh,trim(grp)//'InitAbundance_Silicon', ce%InitAbundance_Silicon )     
  call hdf5_read_attribute(fh,trim(grp)//'InitAbundance_Iron', ce%InitAbundance_Iron )        
  call hdf5_read_attribute(fh,trim(grp)//'CalciumOverSilicon', ce%CalciumOverSilicon )  
  call hdf5_read_attribute(fh,trim(grp)//'SulphurOverSilicon', ce%SulphurOverSilicon )   

#ifndef EAGLE
  grp = 'Parameters/NumericalParameters/'
  call hdf5_read_attribute(fh,trim(grp)//'ComovingIntegrationOn', np%ComovingIntegrationOn )
  call hdf5_read_attribute(fh,trim(grp)//'TypeOfTimestepCriterion', np%TypeOfTimestepCriterion )
  call hdf5_read_attribute(fh,trim(grp)//'ErrTolIntAccuracy', np%ErrTolIntAccuracy )
  call hdf5_read_attribute(fh,trim(grp)//'MaxSizeTimestep', np%MaxSizeTimestep )
  call hdf5_read_attribute(fh,trim(grp)//'MinSizeTimestep', np%MinSizeTimestep )
  call hdf5_read_attribute(fh,trim(grp)//'ErrTolTheta', np%ErrTolTheta )
  call hdf5_read_attribute(fh,trim(grp)//'TypeOfOpeningCriterion', np%TypeOfOpeningCriterion )
  call hdf5_read_attribute(fh,trim(grp)//'ErrTolForceAcc', np%ErrTolForceAcc )
  call hdf5_read_attribute(fh,trim(grp)//'TreeDomainUpdateFrequency', np%TreeDomainUpdateFrequency )
  call hdf5_read_attribute(fh,trim(grp)//'DesNumNgb', np%DesNumNgb )
  call hdf5_read_attribute(fh,trim(grp)//'MaxNumNgbDeviation', np%MaxNumNgbDeviation )
  call hdf5_read_attribute(fh,trim(grp)//'ArtBulkViscConst', np%ArtBulkViscConst )
  call hdf5_read_attribute(fh,trim(grp)//'InitGasU_ERG', np%InitGasU_ERG )
  call hdf5_read_attribute(fh,trim(grp)//'MinGasU_ERG', np%MinGasU_ERG )
  call hdf5_read_attribute(fh,trim(grp)//'CourantFac', np%CourantFac )
  call hdf5_read_attribute(fh,trim(grp)//'PartAllocFactor', np%PartAllocFactor )
  call hdf5_read_attribute(fh,trim(grp)//'TreeAllocFactor', np%TreeAllocFactor )
  call hdf5_read_attribute(fh,trim(grp)//'BufferSize', np%BufferSize )
  call hdf5_read_attribute(fh,trim(grp)//'MinGasHsmlFractional', np%MinGasHsmlFractional )

  grp = 'Parameters/StellarEvolutionParameters/'
  call hdf5_read_attribute(fh,trim(grp)//'SF_EOSGammaEffective', sep%SF_EOSGammaEffective )
  call hdf5_read_attribute(fh,trim(grp)//'SF_EOSEnergyAtThreshold_ERG', sep%SF_EOSEnergyAtThreshold_ERG )
  call hdf5_read_attribute(fh,trim(grp)//'SF_THRESH_MinPhysDens_HpCM3', sep%SF_THRESH_MinPhysDens_HpCM3 ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SF_THRESH_MinOverDens', sep%SF_THRESH_MinOverDens )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SF_THRESH_MetDepExponent', sep%SF_THRESH_MetDepExponent ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SF_THRESH_MetDepFiducialZ', sep%SF_THRESH_MetDepFiducialZ )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SF_THRESH_MetDepMaxPhysDens_HpCM3', sep%SF_THRESH_MetDepMaxPhysDens_HpCM3 ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SF_THRESH_MaxTemp_K', sep%SF_THRESH_MaxTemp_K )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SF_SchmidtLawCoeff_MSUNpYRpKPC2', sep%SF_SchmidtLawCoeff_MSUNpYRpKPC2 )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SF_SchmidtLawExponent', sep%SF_SchmidtLawExponent )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SF_SchmidtLawCoeff_GpSpCM2', sep%SF_SchmidtLawCoeff_GpSpCM2 ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SF_MetDepDensThresholdOn', sep%SF_MetDepDensThresholdOn ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'IMF_Model', sep%IMF_Model )
!!$  call hdf5_read_attribute(fh,trim(grp)//'IMF_LifetimeModel', sep%IMF_LifetimeModel ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'IMF_MinMass_MSUN', sep%IMF_MinMass_MSUN )
!!$  call hdf5_read_attribute(fh,trim(grp)//'IMF_MaxMass_MSUN', sep%IMF_MaxMass_MSUN ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNIa_Model', sep%SNIa_Model )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNIa_Efficiency_fracwd', sep%SNIa_Efficiency_fracwd ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNIa_Energy_ERG', sep%SNIa_Energy_ERG )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNIa_MassTransferOn', sep%SNIa_MassTransferOn )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNIa_EnergyTransferOn', sep%SNIa_EnergyTransferOn ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_MassTransferOn', sep%SNII_MassTransferOn ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_MinMass_MSUN', sep%SNII_MinMass_MSUN )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_MaxMass_MSUN', sep%SNII_MaxMass_MSUN )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_Factor_Hydrogen', sep%SNII_Factor_Hydrogen ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_Factor_Helium', sep%SNII_Factor_Helium )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_Factor_Carbon', sep%SNII_Factor_Carbon )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_Factor_Nitrogen', sep%SNII_Factor_Nitrogen ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_Factor_Oxygen', sep%SNII_Factor_Oxygen ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_Factor_Neon', sep%SNII_Factor_Neon )
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_Factor_Magnesium', sep%SNII_Factor_Magnesium ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_Factor_Silicon', sep%SNII_Factor_Silicon ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'SNII_Factor_Iron', sep%SNII_Factor_Iron )
!!$  call hdf5_read_attribute(fh,trim(grp)//'AGB_MassTransferOn', sep%AGB_MassTransferOn )
!!$  call hdf5_read_attribute(fh,trim(grp)//'POPIII_MassTransferOn', sep%POPIII_MassTransferOn )
!!$  call hdf5_read_attribute(fh,trim(grp)//'POPIII_EnergyTransferOn', sep%POPIII_EnergyTransferOn ) 
!!$  call hdf5_read_attribute(fh,trim(grp)//'POPIII_Energy_ERG', sep%POPIII_Energy_ERG )
!!$  call hdf5_read_attribute(fh,trim(grp)//'POPIII_NumPerMsun', sep%POPIII_NumPerMsun ) 


  grp = 'Parameters/WindParameters/'
  call hdf5_read_attribute(fh,trim(grp)//'SNII_WindOn', wp%SNII_WindOn )
  call hdf5_read_attribute(fh,trim(grp)//'SNII_WindIsotropicOn', wp%SNII_WindIsotropicOn )
  call hdf5_read_attribute(fh,trim(grp)//'SNII_WindSpeed_KMpS', wp%SNII_WindSpeed_KMps )
  call hdf5_read_attribute(fh,trim(grp)//'SNII_WindMassLoading', wp%SNII_WindMassLoading )
  call hdf5_read_attribute(fh,trim(grp)//'SNII_WindDelay_YR', wp%SNII_WindDelay_YR )
#endif


  this%chemical_elements = ce
  this%numerical_parameters = np
  this%stellar_evolution_parameters = sep
  this%wind_parameters = wp



end subroutine gadget_owls_parameters_read_lun

  
!> reads OWLS parameters from an hdf5 file
!--------------------------------------------------------------
subroutine gadget_owls_parameters_read_file(this, snapfile)
  type(gadget_owls_parameters_type), intent(out) :: this
  character(*), intent(in) :: snapfile
  integer :: fh

  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call gadget_owls_parameters_read_lun(this, fh)
  call hdf5_close_file( fh )


end subroutine gadget_owls_parameters_read_file


!> writes owls parameters to an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_owls_parameters_write_lun(this, fh)
  type(gadget_owls_parameters_type), intent(in) :: this
  integer, intent(in) :: fh

  type(gadget_owls_chemical_elements_type) :: ce
  type(gadget_owls_numerical_parameters_type) :: np
  type(gadget_owls_stellar_evolution_parameters_type) :: sep
  type(gadget_owls_wind_parameters_type) :: wp
  character(clen) :: grp

  ce = this%chemical_elements 
  np = this%numerical_parameters 
  sep = this%stellar_evolution_parameters 
  wp = this%wind_parameters 


  grp = 'Parameters/ChemicalElements/'
  call hdf5_write_attribute(fh,trim(grp)//'BG_NELEMENTS', ce%BG_NELEMENTS )  
  call hdf5_write_attribute(fh,trim(grp)//'ElementNames', ce%ElementNames )
  call hdf5_write_attribute(fh,trim(grp)//'InitAbundance_Hydrogen', ce%InitAbundance_Hydrogen )  
  call hdf5_write_attribute(fh,trim(grp)//'InitAbundance_Helium', ce%InitAbundance_Helium )
  call hdf5_write_attribute(fh,trim(grp)//'InitAbundance_Carbon', ce%InitAbundance_Carbon )      
  call hdf5_write_attribute(fh,trim(grp)//'InitAbundance_Nitrogen', ce%InitAbundance_Nitrogen )
  call hdf5_write_attribute(fh,trim(grp)//'InitAbundance_Oxygen', ce%InitAbundance_Oxygen )      
  call hdf5_write_attribute(fh,trim(grp)//'InitAbundance_Neon', ce%InitAbundance_Neon )        
  call hdf5_write_attribute(fh,trim(grp)//'InitAbundance_Magnesium', ce%InitAbundance_Magnesium )   
  call hdf5_write_attribute(fh,trim(grp)//'InitAbundance_Silicon', ce%InitAbundance_Silicon )     
  call hdf5_write_attribute(fh,trim(grp)//'InitAbundance_Iron', ce%InitAbundance_Iron )        
  call hdf5_write_attribute(fh,trim(grp)//'CalciumOverSilicon', ce%CalciumOverSilicon )  
  call hdf5_write_attribute(fh,trim(grp)//'SulphurOverSilicon', ce%SulphurOverSilicon )   

  grp = 'Parameters/NumericalParameters/'
  call hdf5_write_attribute(fh,trim(grp)//'ComovingIntegrationOn', np%ComovingIntegrationOn )
  call hdf5_write_attribute(fh,trim(grp)//'TypeOfTimestepCriterion', np%TypeOfTimestepCriterion )
  call hdf5_write_attribute(fh,trim(grp)//'ErrTolIntAccuracy', np%ErrTolIntAccuracy )
  call hdf5_write_attribute(fh,trim(grp)//'MaxSizeTimestep', np%MaxSizeTimestep )
  call hdf5_write_attribute(fh,trim(grp)//'MinSizeTimestep', np%MinSizeTimestep )
  call hdf5_write_attribute(fh,trim(grp)//'ErrTolTheta', np%ErrTolTheta )
  call hdf5_write_attribute(fh,trim(grp)//'TypeOfOpeningCriterion', np%TypeOfOpeningCriterion )
  call hdf5_write_attribute(fh,trim(grp)//'ErrTolForceAcc', np%ErrTolForceAcc )
  call hdf5_write_attribute(fh,trim(grp)//'TreeDomainUpdateFrequency', np%TreeDomainUpdateFrequency )
  call hdf5_write_attribute(fh,trim(grp)//'DesNumNgb', np%DesNumNgb )
  call hdf5_write_attribute(fh,trim(grp)//'MaxNumNgbDeviation', np%MaxNumNgbDeviation )
  call hdf5_write_attribute(fh,trim(grp)//'ArtBulkViscConst', np%ArtBulkViscConst )
  call hdf5_write_attribute(fh,trim(grp)//'InitGasU_ERG', np%InitGasU_ERG )
  call hdf5_write_attribute(fh,trim(grp)//'MinGasU_ERG', np%MinGasU_ERG )
  call hdf5_write_attribute(fh,trim(grp)//'CourantFac', np%CourantFac )
  call hdf5_write_attribute(fh,trim(grp)//'PartAllocFactor', np%PartAllocFactor )
  call hdf5_write_attribute(fh,trim(grp)//'TreeAllocFactor', np%TreeAllocFactor )
  call hdf5_write_attribute(fh,trim(grp)//'BufferSize', np%BufferSize )
  call hdf5_write_attribute(fh,trim(grp)//'MinGasHsmlFractional', np%MinGasHsmlFractional )


  grp = 'Parameters/StellarEvolutionParameters/'
  call hdf5_write_attribute(fh,trim(grp)//'SF_EOSGammaEffective', sep%SF_EOSGammaEffective )
  call hdf5_write_attribute(fh,trim(grp)//'SF_EOSEnergyAtThreshold_ERG', sep%SF_EOSEnergyAtThreshold_ERG )
  call hdf5_write_attribute(fh,trim(grp)//'SF_THRESH_MinPhysDens_HpCM3', sep%SF_THRESH_MinPhysDens_HpCM3 ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SF_THRESH_MinOverDens', sep%SF_THRESH_MinOverDens )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SF_THRESH_MetDepExponent', sep%SF_THRESH_MetDepExponent ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SF_THRESH_MetDepFiducialZ', sep%SF_THRESH_MetDepFiducialZ )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SF_THRESH_MetDepMaxPhysDens_HpCM3', sep%SF_THRESH_MetDepMaxPhysDens_HpCM3 ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SF_THRESH_MaxTemp_K', sep%SF_THRESH_MaxTemp_K )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SF_SchmidtLawCoeff_MSUNpYRpKPC2', sep%SF_SchmidtLawCoeff_MSUNpYRpKPC2 )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SF_SchmidtLawExponent', sep%SF_SchmidtLawExponent )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SF_SchmidtLawCoeff_GpSpCM2', sep%SF_SchmidtLawCoeff_GpSpCM2 ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SF_MetDepDensThresholdOn', sep%SF_MetDepDensThresholdOn ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'IMF_Model', sep%IMF_Model )
!!$  call hdf5_write_attribute(fh,trim(grp)//'IMF_LifetimeModel', sep%IMF_LifetimeModel ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'IMF_MinMass_MSUN', sep%IMF_MinMass_MSUN )
!!$  call hdf5_write_attribute(fh,trim(grp)//'IMF_MaxMass_MSUN', sep%IMF_MaxMass_MSUN ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNIa_Model', sep%SNIa_Model )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNIa_Efficiency_fracwd', sep%SNIa_Efficiency_fracwd ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNIa_Energy_ERG', sep%SNIa_Energy_ERG )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNIa_MassTransferOn', sep%SNIa_MassTransferOn )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNIa_EnergyTransferOn', sep%SNIa_EnergyTransferOn ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_MassTransferOn', sep%SNII_MassTransferOn ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_MinMass_MSUN', sep%SNII_MinMass_MSUN )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_MaxMass_MSUN', sep%SNII_MaxMass_MSUN )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_Factor_Hydrogen', sep%SNII_Factor_Hydrogen ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_Factor_Helium', sep%SNII_Factor_Helium )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_Factor_Carbon', sep%SNII_Factor_Carbon )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_Factor_Nitrogen', sep%SNII_Factor_Nitrogen ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_Factor_Oxygen', sep%SNII_Factor_Oxygen ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_Factor_Neon', sep%SNII_Factor_Neon )
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_Factor_Magnesium', sep%SNII_Factor_Magnesium ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_Factor_Silicon', sep%SNII_Factor_Silicon ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'SNII_Factor_Iron', sep%SNII_Factor_Iron )
!!$  call hdf5_write_attribute(fh,trim(grp)//'AGB_MassTransferOn', sep%AGB_MassTransferOn )
!!$  call hdf5_write_attribute(fh,trim(grp)//'POPIII_MassTransferOn', sep%POPIII_MassTransferOn )
!!$  call hdf5_write_attribute(fh,trim(grp)//'POPIII_EnergyTransferOn', sep%POPIII_EnergyTransferOn ) 
!!$  call hdf5_write_attribute(fh,trim(grp)//'POPIII_Energy_ERG', sep%POPIII_Energy_ERG )
!!$  call hdf5_write_attribute(fh,trim(grp)//'POPIII_NumPerMsun', sep%POPIII_NumPerMsun ) 

  grp = 'Parameters/WindParameters/'
  call hdf5_write_attribute(fh,trim(grp)//'SNII_WindOn', wp%SNII_WindOn )
  call hdf5_write_attribute(fh,trim(grp)//'SNII_WindIsotropicOn', wp%SNII_WindIsotropicOn )
  call hdf5_write_attribute(fh,trim(grp)//'SNII_WindSpeed_KMpS', wp%SNII_WindSpeed_KMps )
  call hdf5_write_attribute(fh,trim(grp)//'SNII_WindMassLoading', wp%SNII_WindMassLoading )
  call hdf5_write_attribute(fh,trim(grp)//'SNII_WindDelay_YR', wp%SNII_WindDelay_YR )


end subroutine gadget_owls_parameters_write_lun





!!============================================
!!
!!    CONSTANTS
!!
!!============================================


!> reads constants from an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_constants_read_lun(this, fh)
  type(gadget_constants_type), intent(out) :: this
  integer, intent(in) :: fh
  
  call hdf5_read_attribute(fh,'Constants/PI',this%pi)
  call hdf5_read_attribute(fh,'Constants/GAMMA',this%gamma)
  call hdf5_read_attribute(fh,'Constants/GRAVITY',this%gravity)
  call hdf5_read_attribute(fh,'Constants/SOLAR_MASS',this%solar_mass)
  call hdf5_read_attribute(fh,'Constants/SOLAR_LUM',this%solar_lum)
  call hdf5_read_attribute(fh,'Constants/RAD_CONST',this%rad_const)
  call hdf5_read_attribute(fh,'Constants/AVOGADRO',this%avogadro)
  call hdf5_read_attribute(fh,'Constants/BOLTZMANN',this%boltzmann)
  call hdf5_read_attribute(fh,'Constants/GAS_CONST',this%gas_const)
  call hdf5_read_attribute(fh,'Constants/C',this%c)
  call hdf5_read_attribute(fh,'Constants/PLANCK',this%planck)
  call hdf5_read_attribute(fh,'Constants/CM_PER_MPC',this%cm_per_mpc)
  call hdf5_read_attribute(fh,'Constants/PROTONMASS',this%protonmass)
  call hdf5_read_attribute(fh,'Constants/ELECTRONMASS',this%electronmass)
  call hdf5_read_attribute(fh,'Constants/ELECTRONCHARGE',this%electroncharge)
  call hdf5_read_attribute(fh,'Constants/HUBBLE',this%hubble)
  call hdf5_read_attribute(fh,'Constants/T_CMB0',this%t_cmb0)
  call hdf5_read_attribute(fh,'Constants/SEC_PER_MEGAYEAR',this%sec_per_megayear)
  call hdf5_read_attribute(fh,'Constants/SEC_PER_YEAR',this%sec_per_year)
  call hdf5_read_attribute(fh,'Constants/STEFAN',this%stefan)
  call hdf5_read_attribute(fh,'Constants/THOMPSON',this%thompson)
  call hdf5_read_attribute(fh,'Constants/EV_TO_ERG',this%ev_to_erg)
  call hdf5_read_attribute(fh,'Constants/Z_Solar',this%z_solar)
  

end subroutine gadget_constants_read_lun

  
!> reads constants from an hdf5 file
!--------------------------------------------------------------
subroutine gadget_constants_read_file(this, snapfile)
  type(gadget_constants_type), intent(out) :: this
  character(*), intent(in) :: snapfile
  integer :: fh

  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call gadget_constants_read_lun(this, fh)
  call hdf5_close_file( fh )


end subroutine gadget_constants_read_file


!> writes constants to an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_constants_write_lun(this, fh)
  type(gadget_constants_type), intent(in) :: this
  integer, intent(in) :: fh

  call hdf5_write_attribute(fh,'Constants/PI',this%pi)
  call hdf5_write_attribute(fh,'Constants/GAMMA',this%gamma)
  call hdf5_write_attribute(fh,'Constants/GRAVITY',this%gravity)
  call hdf5_write_attribute(fh,'Constants/SOLAR_MASS',this%solar_mass)
  call hdf5_write_attribute(fh,'Constants/SOLAR_LUM',this%solar_lum)
  call hdf5_write_attribute(fh,'Constants/RAD_CONST',this%rad_const)
  call hdf5_write_attribute(fh,'Constants/AVOGADRO',this%avogadro)
  call hdf5_write_attribute(fh,'Constants/BOLTZMANN',this%boltzmann)
  call hdf5_write_attribute(fh,'Constants/GAS_CONST',this%gas_const)
  call hdf5_write_attribute(fh,'Constants/C',this%c)
  call hdf5_write_attribute(fh,'Constants/PLANCK',this%planck)
  call hdf5_write_attribute(fh,'Constants/CM_PER_MPC',this%cm_per_mpc)
  call hdf5_write_attribute(fh,'Constants/PROTONMASS',this%protonmass)
  call hdf5_write_attribute(fh,'Constants/ELECTRONMASS',this%electronmass)
  call hdf5_write_attribute(fh,'Constants/ELECTRONCHARGE',this%electroncharge)
  call hdf5_write_attribute(fh,'Constants/HUBBLE',this%hubble)
  call hdf5_write_attribute(fh,'Constants/T_CMB0',this%t_cmb0)
  call hdf5_write_attribute(fh,'Constants/SEC_PER_MEGAYEAR',this%sec_per_megayear)
  call hdf5_write_attribute(fh,'Constants/SEC_PER_YEAR',this%sec_per_year)
  call hdf5_write_attribute(fh,'Constants/STEFAN',this%stefan)
  call hdf5_write_attribute(fh,'Constants/THOMPSON',this%thompson)
  call hdf5_write_attribute(fh,'Constants/EV_TO_ERG',this%ev_to_erg)
  call hdf5_write_attribute(fh,'Constants/Z_Solar',this%z_solar)


end subroutine gadget_constants_write_lun





!!============================================
!!
!!    DATA ATTRIBUTES
!!
!!============================================



!> writes data attributes to an hdf5 file  
!------------------------------------------------    
subroutine gadget_data_attributes_write_lun( attr, fh, grp_tag, dat_tag  )
  type(gadget_one_data_attributes_type) :: attr
  integer(i4b) :: fh
  character(*) :: grp_tag
  character(*) :: dat_tag

  character(clen) :: tag

  tag = trim(grp_tag) // trim(dat_tag) // '/CGSConversionFactor'
  call hdf5_write_attribute( fh, tag, attr%CGSConversionFactor )

  tag = trim(grp_tag) // trim(dat_tag) // '/h-scale-exponent'
  call hdf5_write_attribute( fh, tag, attr%h_scale_exponent )

  tag = trim(grp_tag) // trim(dat_tag) // '/aexp-scale-exponent'
  call hdf5_write_attribute( fh, tag, attr%aexp_scale_exponent )

  tag = trim(grp_tag) // trim(dat_tag) // '/VarDescription'
  call hdf5_write_attribute( fh, tag, attr%VarDescription )


end subroutine gadget_data_attributes_write_lun



!> reads data attributes from an hdf5 file  
!------------------------------------------------    
subroutine gadget_data_attributes_read_lun( attr, fh, grp_tag, dat_tag  )
  type(gadget_one_data_attributes_type) :: attr
  integer(i4b) :: fh
  character(*) :: grp_tag
  character(*) :: dat_tag

  character(clen) :: tag

  tag = trim(grp_tag) // trim(dat_tag) // '/CGSConversionFactor'
  call hdf5_read_attribute( fh, tag, attr%CGSConversionFactor )

  tag = trim(grp_tag) // trim(dat_tag) // '/h-scale-exponent'
  call hdf5_read_attribute( fh, tag, attr%h_scale_exponent )

  tag = trim(grp_tag) // trim(dat_tag) // '/aexp-scale-exponent'
  call hdf5_read_attribute( fh, tag, attr%aexp_scale_exponent )

  tag = trim(grp_tag) // trim(dat_tag) // '/VarDescription'
  call hdf5_read_attribute( fh, tag, attr%VarDescription )

end subroutine gadget_data_attributes_read_lun





!> sets data attributes using supplied units
!------------------------------------------------    
subroutine gadget_data_attributes_set(gunits)
  type(gadget_units_type), intent(in) :: gunits

  gattrs%pos%CGSConversionFactor = gunits%cgs_length
  gattrs%vel%CGSConversionFactor = gunits%cgs_velocity
  gattrs%id%CGSConversionFactor = 1.0d0
  gattrs%mass%CGSConversionFactor = gunits%cgs_mass
  gattrs%T%CGSConversionFactor = 1.0d0
  gattrs%rho%CGSConversionFactor = gunits%cgs_density
  gattrs%ye%CGSConversionFactor = 1.0d0
  gattrs%xHI%CGSConversionFactor = 1.0d0
  gattrs%hsml%CGSConversionFactor = gunits%cgs_length
  gattrs%lasthit%CGSConversionFactor = 1.0d0

  gattrs%xHI_cloudy%CGSConversionFactor = 1.0d0
  gattrs%Hmf%CGSConversionFactor = 1.0d0
  gattrs%xHeI%CGSConversionFactor = 1.0d0
  gattrs%xHeII%CGSConversionFactor = 1.0d0
  gattrs%Hemf%CGSConversionFactor = 1.0d0
  gattrs%gammaHI%CGSConversionFactor = 1.0d0
  gattrs%eos%CGSConversionFactor = 1.0d0
  gattrs%sfr%CGSConversionFactor = 1.0d0
  gattrs%fh2%CGSConversionFactor = 1.0d0
  gattrs%nH%CGSConversionFactor = 1.0d0



  gattrs%pos%VarDescription = "Co-moving coordinates. Physical: r = a x = Coordinates h^-1 a U_L [cm]"
  gattrs%vel%VarDescription = "Co-moving velocities. Physical v_p = a dx/dt  = Velocities a^1/2 U_V [cm/s]"
  gattrs%id%VarDescription = "Unique particle identifier"
  gattrs%mass%VarDescription = "Particle mass. Physical m = Mass h^-1 U_M [g]"
  gattrs%T%VarDescription = "Temperature [K]"
  gattrs%rho%VarDescription = "Co-moving mass densities. Physical rho = Densities h^2 a^-3 U_M U_L^-3 [g/cm^3]"
  gattrs%ye%VarDescription = "Electron abundance = n_e / n_H"
  gattrs%xHI%VarDescription = "Hydrogen neutral fraction = n_HI / n_H"
  gattrs%hsml%VarDescription = "Co-moving smoothing length. Physical h = SmoothingLength h^-1 a U_L [cm]"
  gattrs%lasthit%VarDescription = "Index of last ray to cross this particle"

  gattrs%xHI_cloudy%VarDescription = "Hydrogen neutral fraction from Cloudy Tables xHI = n_HI / n_H"
  gattrs%Hmf%VarDescription = "Hydrogen mass fraction = mass in atomic Hydrogen / particle mass"
  gattrs%xHeI%VarDescription = "Helium I fraction = n_HeI / n_He"
  gattrs%xHeII%VarDescription = "Helium II fraction = n_HeII / n_He"
  gattrs%Hemf%VarDescription = "Helium mass fraction = mass in atomic Helium / particle mass"
  gattrs%gammaHI%VarDescription = "Hydrogen photoionization rate [1/s]"
  gattrs%eos%VarDescription = "Flag for effective equation of state. 1 if currently on EoS" // &
       ", 0 if has never been on EoS, -ExpansionFactor if left the EoS at ExpansionFactor"
  gattrs%sfr%VarDescription = "Star formation rate. Physical sfr = StarformationRate SOLAR_MASS SEC_PER_YEAR^-1 [g/s]"
  gattrs%fh2%VarDescription = "Molecular Hydrogen Mass Fraction. Mass in molecules over total Hydrogen mass"
  gattrs%nH%VarDescription = "Physical Atomic Hydrogen Number Density [cm^-3]. No factors of h, no factors of a, nHI = nH * xHI"



  gattrs%pos%h_scale_exponent = -1.0d0
  gattrs%vel%h_scale_exponent = 0.0d0
  gattrs%id%h_scale_exponent = 0.0d0
  gattrs%mass%h_scale_exponent = -1.0d0
  gattrs%T%h_scale_exponent = 0.0d0
  gattrs%rho%h_scale_exponent = 2.0d0
  gattrs%ye%h_scale_exponent = 0.0d0
  gattrs%xHI%h_scale_exponent = 0.0d0
  gattrs%hsml%h_scale_exponent = -1.0d0
  gattrs%lasthit%h_scale_exponent = 0.0d0

  gattrs%xHI_cloudy%h_scale_exponent = 0.0d0
  gattrs%Hmf%h_scale_exponent = 0.0d0
  gattrs%xHeI%h_scale_exponent = 0.0d0
  gattrs%xHeII%h_scale_exponent = 0.0d0
  gattrs%Hemf%h_scale_exponent = 0.0d0
  gattrs%gammaHI%h_scale_exponent = 0.0d0
  gattrs%eos%h_scale_exponent = 0.0d0
  gattrs%sfr%h_scale_exponent = 0.0d0
  gattrs%fh2%h_scale_exponent = 0.0d0
  gattrs%nH%h_scale_exponent = 0.0d0



  gattrs%pos%aexp_scale_exponent = 1.0d0
  gattrs%vel%aexp_scale_exponent = 0.5d0
  gattrs%id%aexp_scale_exponent = 0.0d0
  gattrs%mass%aexp_scale_exponent = 0.0d0
  gattrs%T%aexp_scale_exponent = 0.0d0
  gattrs%rho%aexp_scale_exponent = -3.0d0
  gattrs%ye%aexp_scale_exponent = 0.0d0
  gattrs%xHI%aexp_scale_exponent = 0.0d0
  gattrs%hsml%aexp_scale_exponent = 1.0d0
  gattrs%lasthit%aexp_scale_exponent = 0.0d0

  gattrs%xHI_cloudy%aexp_scale_exponent = 0.0d0
  gattrs%Hmf%aexp_scale_exponent = 0.0d0
  gattrs%xHeI%aexp_scale_exponent = 0.0d0
  gattrs%xHeII%aexp_scale_exponent = 0.0d0
  gattrs%Hemf%aexp_scale_exponent = 0.0d0
  gattrs%gammaHI%aexp_scale_exponent = 0.0d0
  gattrs%eos%aexp_scale_exponent = 0.0d0
  gattrs%sfr%aexp_scale_exponent = 0.0d0
  gattrs%fh2%aexp_scale_exponent = 0.0d0
  gattrs%nH%aexp_scale_exponent = 0.0d0


end subroutine gadget_data_attributes_set











!> forms a snapshot file name from, path, base, SnapNum, and FileNum.
!! For example, path = "/home/galtay/data/snapshots", base = "snap",
!! SnapNum = 2 and FileNum = 15 would return in SnapFile, 
!! "/home/galtay/data/snapshots/snap_002.15".  Setting the hdf5bool to 
!! true appends ".hdf5" to the file. 
!------------------------------------------------------------------------
subroutine form_gadget_snapshot_file_name(path, base, FileNum, &
     SnapFile, hdf5bool)

  character(*), intent(in)     :: path        !< path to snapshot dir
  character(*), intent(in)     :: base        !< base snapshot name
  integer(i4b), intent(in)     :: FileNum     !< file number of snapshot
  character(clen), intent(out) :: SnapFile    !< snapshot filename
  logical, intent(in)          :: hdf5bool    !< hdf5 file? 

  character(clen) :: SnapFileTmp
  character(10) :: FileNumChar
  character(clen) :: fmt
  logical :: Fthere

  write(FileNumChar,"(I6)") FileNum
  fmt = "(A,'/',A)"

  ! first write a file with no extension
  !--------------------------------------
  write(SnapFileTmp,fmt) trim(path), trim(base)
  SnapFile = trim(SnapFileTmp)
  if (hdf5bool) SnapFile = trim(SnapFile) // ".hdf5"
  inquire( file=SnapFile, exist=Fthere )


  ! if the file number is 0 and a file with no extension exists then return
  ! otherwise append the FileNum to the file name
  !-------------------------------------------------------------------------
  if (FileNum == 0 .and. Fthere) then
     return
  else
     SnapFile = trim(SnapFileTmp) // "." // trim(adjustl(FileNumChar))
     if (hdf5bool) SnapFile = trim(SnapFile) // ".hdf5"
  end if


end subroutine form_gadget_snapshot_file_name










end module gadget_general_class


