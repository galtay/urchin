!> \file config.F90

!> \brief the module that handles the config variables.
!<
module config_mod
use myf90_mod
use hdf5_wrapper
implicit none
private


public :: read_config_file
public :: write_config_hdf5_lun
public :: CV


! these variables are read in directly from the config file
!----------------------------------------------------------
type config_variables_type

   character(clen) :: config_file   !< name of the config file

   integer(i4b)    :: Verbosity !< 0=silent, 1=whisper, 2=talk, 3=debug
   logical         :: JustInit  !< set true to stop after initialization
   logical         :: Comoving  !< convert input from co-moving to physical? 
   integer(i8b)    :: IntSeed   !< seed for mersenne twister

   logical         :: UseIsoTemp !< fix all particles to single temperature? 
   real(r8b)       :: IsoTemp    !< if UseIsoTemp, all pars fixed @ T=IsoTemp

   logical         :: AdjustIsmPars       !< adjust temperature and hot fraction of star forming particles?
   real(r8b)       :: IsmTemp             !< if AdjustIsmPars, Cold Phase T of ISM particles = IsmTemp 
   real(r8b)       :: IsmHotFrac          !< if AdjustIsmPars, mass fraction of ISM gas in hot phase
   logical         :: IsmAddH2            !< add molecular hydrogen to ISM particles? 
   real(r8b)       :: IsmH2Param1         !< P0 from Blitz 06
   real(r8b)       :: IsmH2Param2         !< alpha from Blitz 06

   logical         :: AdjustShldPars      !< adjust temperature of self-shielded particles?
   real(r8b)       :: ShldTemp            !< if AdjustShldPars, particles with tau_eff > ShldTauThresh get T=ShldTemp
   real(r8b)       :: ShldTauThresh       !< tau_eff threshold to use ShldTemp

   character(clen) :: InputType           !< type of Gadget, one of {"public", "cosmobh", "owls",
                                          !< "vbromm", "publichdf5", "eagle"}
   character(clen) :: EOSvarName          !< for owls/eagle sims one of {"SfFlag", "OnEquationOfState"}
   character(clen) :: SnapPath            !< dir containing particle snapshots
   character(clen) :: ParFileBase         !< particle snapshot file base
   integer(i4b)    :: ParFilesPerSnap     !< files per particle snapshot

   character(clen) :: ItabDir             !< file containing background spectrum
   character(clen) :: AtomicRatesFile     !< file containing atomic rates tables
   character(clen) :: SlabLookupFile      !< file containing analytic slab solutions

   logical         :: Monochromatic       !< monochromatic radiation? 
   real(r8b)       :: MonoRydbergs        !< energy of monochromatic photons
   character(clen) :: MonoNormMethod      !< method to normalize monochromatic flux to HM01 spectrum
                                          !< one of ['ngamma', 'gammaHI']

   logical         :: DoSelection         !< trace a sub-selection and output full data on it
   real(r8b)       :: SelectRadius        !< radius of spehre to select particles in
   real(r8b)       :: SelectXcoord        !< x-posiiton of sphere
   real(r8b)       :: SelectYcoord        !< y-posiiton of sphere
   real(r8b)       :: SelectZcoord        !< z-posiiton of sphere
 
   integer(i4b)    :: XcoordBndryCond     !< one of {-1:reflecting 0:vacuum 1:periodic}
   integer(i4b)    :: YcoordBndryCond     !< one of {-1:reflecting 0:vacuum 1:periodic}
   integer(i4b)    :: ZcoordBndryCond     !< one of {-1:reflecting 0:vacuum 1:periodic}
 
   logical         :: RayStats            !< T = mini output file on ray statistics in raystats.dat
   real(r8b)       :: MaxRayDist          !< if > 0, max ray distance in physical code units

   integer(i4b)    :: HealpixNside        !< n pixels = 12 * nside^2
   integer(i4b)    :: NumIterations       !< number of iterations over particles
   real(r8b)       :: ConvTol             !< convergence tolerance

   real(r8b)       :: SpecLoRyd           !< lowest energy photon to include in spectrum [Rydbergs]
   real(r8b)       :: SpecHiRyd           !< highest energy photon to include in spectrum [Rydbergs]
   real(r8b)       :: GammaMultiplier     !< factor by which to multiply background flux
   integer(i4b)    :: FluxDistribution    !< 0=isotropic, +/-3 = incident from +/- z direction
   real(r8b)       :: RhoMultiplier       !< factor by which to multiply rho (and mass) of particles
   real(r8b)       :: NeBackground        !< constant electron number density from helium and metals

   real(r8b)       :: H_mf                !< if not in snapshots, use this hydrogen mass fraction
   real(r8b)       :: He_mf               !< if not in snapshots, use this helium mass fraction
   
   character(clen) :: OutputDir           !< path to output directory
   character(clen) :: OutputFileBase      !< output file base
   character(clen) :: OutputType          !< one of {"binary", "hdf5"}

   integer(i4b)    :: PartPerCell         !< minimum particles in a tree leaf

end type config_variables_type


type(config_variables_type) :: CV           !< config variables          




contains


!> reads the config file 
!==============================
subroutine read_config_file(config_file)
  character(clen), parameter :: myname="read_config_file"
  logical, parameter :: crash = .true.
  character(clen), intent(in) :: config_file  !< file to read config vars from

  character(clen) :: keyword
  character(clen) :: str
  logical :: file_exists


    write(str,'(A,T27,A)') 'using configuration file: ', trim(config_file)
    call mywrite(str,verb=0) 

    inquire( file=config_file, exist=file_exists )
    if (.not. file_exists) then
       call myerr("config file does not exist",myname,crash)
    end if


    !   initialization
    !------------------------------
    keyword = "Verbosity:"
    call scanfile(config_file,keyword,CV%Verbosity)
    myf90_verbosity = CV%verbosity

    keyword = "JustInit:"
    call scanfile(config_file,keyword,CV%JustInit)

    keyword = "Comoving:"
    call scanfile(config_file,keyword,CV%Comoving)

    keyword = "IntSeed:"
    call scanfile(config_file,keyword,CV%IntSeed)


    !   special treatment
    !------------------------------

    keyword = "UseIsoTemp:"
    call scanfile(config_file,keyword,CV%UseIsoTemp)

    keyword = "IsoTemp:"
    call scanfile(config_file,keyword,CV%IsoTemp)


    keyword = "AdjustIsmPars:"
    call scanfile(config_file,keyword,CV%AdjustIsmPars)

    keyword = "IsmTemp:"
    call scanfile(config_file,keyword,CV%IsmTemp)

    keyword = "IsmHotFrac:"
    call scanfile(config_file,keyword,CV%IsmHotFrac)

    keyword = "IsmAddH2:"
    call scanfile(config_file,keyword,CV%IsmAddH2)

    keyword = "IsmH2Param1:"
    call scanfile(config_file,keyword,CV%IsmH2Param1)

    keyword = "IsmH2Param2:"
    call scanfile(config_file,keyword,CV%IsmH2Param2)


    keyword = "AdjustShldPars:"
    call scanfile(config_file,keyword,CV%AdjustShldPars)

    keyword = "ShldTemp:"
    call scanfile(config_file,keyword,CV%ShldTemp)

    keyword = "ShldTauThresh:"
    call scanfile(config_file,keyword,CV%ShldTauThresh)



    !   input snapshot information
    !------------------------------
    keyword = "InputType:"
    call scanfile(config_file,keyword,CV%InputType)

    keyword = "EOSvarName:"
    call scanfile(config_file,keyword,CV%EOSvarName)

    keyword = "SnapPath:"
    call scanfile(config_file,keyword,CV%SnapPath)

    keyword = "ParFileBase:"
    call scanfile(config_file,keyword,CV%ParFileBase)

    keyword = "ParFilesPerSnap:"
    call scanfile(config_file,keyword,CV%ParFilesPerSnap)


    keyword = "ItabDir:"
    call scanfile(config_file,keyword,CV%ItabDir)

    keyword = "AtomicRatesFile:"
    call scanfile(config_file,keyword,CV%AtomicRatesFile)

    keyword = "SlabLookupFile:"
    call scanfile(config_file,keyword,CV%SlabLookupFile)

    !   monochromatic spectrum
    !----------------------------
    keyword = "Monochromatic:"
    call scanfile(config_file,keyword,CV%Monochromatic)

    keyword = "MonoRydbergs:"
    call scanfile(config_file,keyword,CV%MonoRydbergs)

    keyword = "MonoNormMethod:"
    call scanfile(config_file,keyword,CV%MonoNormMethod)


    !   if extracting particles in spehre
    !----------------------------
    keyword = "DoSelection:"
    call scanfile(config_file,keyword,CV%DoSelection)

    keyword = "SelectRadius:"
    call scanfile(config_file,keyword,CV%SelectRadius)

    keyword = "SelectXcoord:"
    call scanfile(config_file,keyword,CV%SelectXcoord)

    keyword = "SelectYcoord:"
    call scanfile(config_file,keyword,CV%SelectYcoord)

    keyword = "SelectZcoord:"
    call scanfile(config_file,keyword,CV%SelectZcoord)


    !   ray tracing
    !---------------------------- 
    keyword = "XcoordBndryCond:"
    call scanfile(config_file,keyword,CV%XcoordBndryCond)

    keyword = "YcoordBndryCond:"
    call scanfile(config_file,keyword,CV%YcoordBndryCond)

    keyword = "ZcoordBndryCond:"
    call scanfile(config_file,keyword,CV%ZcoordBndryCond)

    keyword = "RayStats:"
    call scanfile(config_file,keyword,CV%RayStats)

    keyword = "MaxRayDist:"
    call scanfile(config_file,keyword,CV%MaxRayDist)

    keyword = "HealpixNside:"
    call scanfile(config_file,keyword,CV%HealpixNside)

    keyword = "NumIterations:"
    call scanfile(config_file,keyword,CV%NumIterations)

    keyword = "ConvTol:"
    call scanfile(config_file,keyword,CV%ConvTol)



    !   ion/temp solving
    !----------------------------
    keyword = "SpecLoRyd:"
    call scanfile(config_file,keyword,CV%SpecLoRyd)

    keyword = "SpecHiRyd:"
    call scanfile(config_file,keyword,CV%SpecHiRyd)

    keyword = "GammaMultiplier:"
    call scanfile(config_file,keyword,CV%GammaMultiplier)

    keyword = "FluxDistribution:"
    call scanfile(config_file,keyword,CV%FluxDistribution)

    keyword = "RhoMultiplier:"
    call scanfile(config_file,keyword,CV%RhoMultiplier)

    keyword = "NeBackground:"
    call scanfile(config_file,keyword,CV%NeBackground)



    !   gas composition
    !----------------------------
    keyword = "H_mf:"
    call scanfile(config_file,keyword,CV%H_mf)

    keyword = "He_mf:"
    call scanfile(config_file,keyword,CV%He_mf)


    !   output
    !-------------
    keyword = "OutputDir:"
    call scanfile(config_file,keyword,CV%OutputDir)

    keyword = "OutputFileBase:"
    call scanfile(config_file,keyword,CV%OutputFileBase)

    keyword = "OutputType:"
    call scanfile(config_file,keyword,CV%OutputType)


    !  tree building
    !--------------------
    keyword = "PartPerCell:"
    call scanfile(config_file,keyword,CV%PartPerCell)


    CV%config_file = config_file

    call dummy_check_config_variables()
    call config_info_to_file()


end subroutine read_config_file


!> run dummy checks on config variables
!========================================
subroutine dummy_check_config_variables()

  logical :: config_good
!  logical :: charmatch

  config_good = .true. 


  if (CV%UseIsoTemp) then
     if (CV%IsoTemp <= 0.0d0) then
        write(*,*)
        write(*,*) "You have set an IsoTemp run, therefore IsoTemp must be > 0.0"
        write(*,*) 
        config_good = .false. 
     endif
  endif
       
  
  if (CV%Monochromatic) then
     select case ( trim(adjustl(CV%MonoNormMethod)) )
     case( 'ngamma' )
     case( 'gammaHI' )
     case default
        write(*,*) 
        write(*,*) "You have set Monochromatic = True in the config file"
        write(*,*) "For Monochromatic runs, "
        write(*,*) "MonoNormMethod must be set to either 'ngamma' or 'gammaHI'"
        write(*,*) "MonoNormMethod = ", trim(CV%MonoNormMethod)
        config_good = .false. 
        write(*,*) 
     end select
  endif

#ifdef incEOS
  if ( trim(CV%InputType) /= 'owls' .and. trim(CV%InputType) /= 'eagle') then
     write(*,*) 
     write(*,*) "Special treatment of ISM not supported with the current input type: "
     write(*,*) trim(CV%InputType)
     write(*,*) 
     config_good = .false.
  endif

#else
  if (CV%AdjustIsmPars) then
     write(*,*) 
     write(*,*) "The Makefile macro incEOS must be defined"
     write(*,*) "to use special treatment of ISM particles"
     write(*,*) 
     config_good = .false.
  endif

  if (CV%IsmAddH2) then
     write(*,*) 
     write(*,*) "The Makefile macro incEOS must be defined"
     write(*,*) "to add H2 treatment to the ISM"
     write(*,*) 
     config_good = .false.
  endif

#endif
  



  if (CV%AdjustIsmPars) then
     if (CV%IsmTemp <= 0.0d0) then
        write(*,*)
        write(*,*) "You have chosen special treatment of ISM particles. "
        write(*,*) "IsmTemp must be > 0.0"
        write(*,*) 
        config_good = .false. 
     endif
     if (CV%IsmHotFrac < 0.0d0) then
        write(*,*)
        write(*,*) "You have chosen special treatment of ISM particles. "
        write(*,*) "IsmHotFrac must be >= 0.0"
        write(*,*) 
        config_good = .false. 
     endif
  endif


  select case ( trim(adjustl(CV%InputType)) )
  case( 'public' )
  case( 'cosmobh' )
  case( 'owls' )
  case( 'vbromm' )
  case( 'publichdf5' )
  case( 'eagle' )
  case default
     write(*,*) "Input Type ", trim(CV%InputType), " not recognized"
     write(*,*) "must be one of"
     write(*,*) "'public':     (Gadget 2.0 Public Standard)"
     write(*,*) "'cosmobh':    (Gadget CosmoBH)"
     write(*,*) "'owls':       (Gadget OWLS/GIMIC)"
     write(*,*) "'vbromm':     (Gadget Bromm)"
     write(*,*) "'publichdf5': (Gadget Public HDF5)"
     write(*,*) "'eagle':      (Gadget EAGLE)"
     config_good = .false. 
  end select
  

  select case ( trim(adjustl(CV%EOSvarName)) )
  case( 'SfFlag' )
  case( 'OnEquationOfState' )
  case default
     write(*,*) "EOS variable name ", trim(CV%EOSvarName), " not recognized"
     write(*,*) "must be one of"
     write(*,*) "'SfFlag'"
     write(*,*) "'OnEquationOfState'"
     config_good = .false. 
  end select


  select case ( trim(adjustl(CV%OutputType)) )
  case( 'binary' )
  case( 'hdf5' )
  case default
     write(*,*) "Output Type ", CV%OutputType, " not recognized"
     write(*,*) "must be 1 (Standard Gadget) or 2 (HDF5)"
     config_good = .false. 
  end select
  
  
#ifdef incH2
  if ( trim(CV%InputType) /= 'owls' .and. trim(CV%InputType) /= 'eagle' ) then
     write(*,*) 
     write(*,*) "Addition of H_2 not supported with the current input type: "
     write(*,*) trim(CV%InputType)
     write(*,*) 
     config_good = .false.
  endif
     

  if (.not. CV%IsmAddH2) then
     write(*,*)
     write(*,*) "You have defined the incH2 macro in the Makefile, but"
     write(*,*) "IsmAddH2 is False in the config file.  Please make these"
     write(*,*) "consistent by commenting out incH2 in the Makefile or "
     write(*,*) "setting IsmAddH2 True."
     write(*,*) 
     config_good = .false. 
  endif
#endif
  

#ifndef incH2
  if (CV%IsmAddH2) then
     write(*,*)
     write(*,*) "You have not defined the incH2 macro in the Makefile, but"
     write(*,*) "IsmAddH2 is True in the config file.  Please make these"
     write(*,*) "consistent by un-commenting incH2 in the Makefile or "
     write(*,*) "setting IsmAddH2 False."
     write(*,*) 
     config_good = .false. 
  endif
#endif




#ifdef incHmf
  if (CV%H_mf < 0.0) then
     write(*,*)
     write(*,*) "You have defined the incHmf macro in the Makefile, and"
     write(*,*) "the Hydrogen mass fraction is set < 0.0 in the config file."
     write(*,*) "This means the H mass fractions should be present in the "
     write(*,*) "particle snapshot your are using! Just a dummy check."
     write(*,*) 
  end if
#else
  if (CV%H_mf < 0.0) then
     write(*,*) "You have not defined the incHmf macro in the Makefile, and"
     write(*,*) "the Hydrogen mass fraction is < 0.0 in the config file."
     write(*,*) "I don't know what the Hydrogen mass fraction should be" 
     write(*,*) "If they are in the snapshot, uncomment incHmf.  If you want"
     write(*,*) "to set a constant value, do so in the config file" 
     config_good = .false.
  end if
#endif




#ifdef incHemf
#ifndef incHe
  write(*,*) "You have included He mass fractions but not included He."
  write(*,*) "Either comment out incHemf or uncomment incHe"
  config_good = .false.
#endif
#endif

#ifdef incHemf
  if (CV%He_mf < 0.0) then
     write(*,*)
     write(*,*) "You have defined the incHemf macro in the Makefile, and"
     write(*,*) "the Helium mass fraction is set < 0.0 in the config file."
     write(*,*) "This means the He mass fractions should be present in the "
     write(*,*) "particle snapshot your are using! Just a dummy check."
     write(*,*) 
  end if
#else
  if (CV%He_mf > 0.0) then
     write(*,*) "You have not defined the incHe macro in the Makefile, but"
     write(*,*) "the Helium mass fraction is > 0.0 in the config file."
     write(*,*) "Please uncomment the incHe line in the Makefile or set the"
     write(*,*) "Helium mass fraction to a negative number in the config file."
     config_good = .false.
  end if
#endif


#ifdef EAGLE
  if(trim(CV%InputType) /= 'eagle') then
     write (*,*) ' EAGLE Makefile macro set but using InputType = ',trim(CV%InputType)
     write (*,*) ' .... this is not alllowed'
     config_good = .false.
  endif
#else
  if ( trim(CV%InputType) == 'eagle' ) then
     write(*,*) ' EAGLE Makefile macro not set, yet using InputType = eagle'
     write(*,*) ' .... this is not allowed'
     config_good = .false.
  endif
#endif


  if (.not. config_good) then
     write(*,*) 'please edit'
     write(*,*) trim(CV%config_file)
     write(*,*) 'or recompile with different Makefile macros'
     stop
  end if


end subroutine dummy_check_config_variables




!> writes the configuration file information to the output directory
!! log file
!====================================================================
subroutine config_info_to_file()

  character(clen) :: config_log_file
  integer(i4b) :: lun
  character(clen) :: i1fmt
  character(clen) :: i3fmt
  character(clen) :: lfmt
  character(clen) :: dfmt
  character(clen) :: sfmt
  character(clen) :: str

  config_log_file = trim(CV%OutputDir) // "/config_values_used.log"
  
  call open_formatted_file_w(config_log_file, lun)

!  105 format(T2,A,I3.3)
!  111 format(T2,A,I10)
!  120 format(T2,A,F10.3)


  i3fmt = '(I3.3)'
  i1fmt = '(I8)' 
  lfmt = '(L1)'
  dfmt = '(ES15.4)'
  sfmt = '(A,T25,A)'



  write(lun,*)"====================================================="
  write(lun,*)"The URCHIN configuration variables used for this run "
  write(lun,*)"====================================================="
  write(lun,*) 

  write(str,i1fmt) CV%Verbosity
  write(lun,sfmt)"Verbosity:", trim(adjustl(str))
  write(lun,*)

  write(str,lfmt) CV%JustInit
  write(lun,sfmt)"JustInit:", trim(adjustl(str))
  write(str,lfmt) CV%Comoving
  write(lun,sfmt)"Comoving:", trim(adjustl(str))
  write(str,'(I20)') CV%IntSeed
  write(lun,sfmt)"IntSeed:", trim(adjustl(str))
  write(lun,*) 

  write(str,lfmt) CV%UseIsoTemp
  write(lun,sfmt)"UseIsoTemp:", trim(adjustl(str))
  write(str,dfmt) CV%IsoTemp
  write(lun,sfmt)"IsoTemp:", trim(adjustl(str))
  write(lun,*)

  write(str,lfmt) CV%AdjustIsmPars
  write(lun,sfmt)"AdjustIsmPars:", trim(adjustl(str))
  write(str,dfmt) CV%IsmTemp
  write(lun,sfmt)"IsmTemp:", trim(adjustl(str))
  write(str,dfmt) CV%IsmHotFrac
  write(lun,sfmt)"IsmHotFrac:", trim(adjustl(str))
  write(str,lfmt) CV%IsmAddH2
  write(lun,sfmt)"IsmAddH2:", trim(adjustl(str))
  write(str,dfmt) CV%IsmH2Param1
  write(lun,sfmt) "IsmH2Param1:", trim(adjustl(str))
  write(str,dfmt) CV%IsmH2Param2
  write(lun,sfmt) "IsmH2Param2:", trim(adjustl(str))
  write(lun,*)

  write(str,lfmt) CV%AdjustShldPars
  write(lun,sfmt)"AdjustShldPars:", trim(adjustl(str))
  write(str,dfmt) CV%ShldTemp
  write(lun,sfmt)"ShldTemp:", trim(adjustl(str))
  write(str,dfmt) CV%ShldTauThresh
  write(lun,sfmt)"ShldTauThresh:", trim(adjustl(str))
  write(lun,*) 

  write(lun,sfmt)"InputType:", trim(CV%InputType)
  write(lun,sfmt)"EOSvarName:", trim(CV%EOSvarName)
  write(lun,sfmt)"SnapPath:", trim(CV%SnapPath)
  write(lun,sfmt)"ParFileBase:", trim(CV%ParFileBase)
  write(str,i1fmt) CV%ParFilesPerSnap
  write(lun,sfmt) "ParFilesPerSnap:", trim(adjustl(str))
  write(lun,*)

  write(lun,sfmt)"ItabDir:", trim(CV%ItabDir)
  write(lun,sfmt)"AtomicRatesFile:", trim(CV%AtomicRatesFile)
  write(lun,sfmt)"SlabLookupFile:", trim(CV%SlabLookupFile)
  write(lun,*)

  write(str,lfmt) CV%Monochromatic
  write(lun,sfmt)"Monochromatic:", trim(adjustl(str))
  write(str,dfmt) CV%MonoRydbergs
  write(lun,sfmt) "MonoRydbergs:", trim(adjustl(str))
  write(lun,sfmt)"MonoNormMethod:", trim(CV%MonoNormMethod)
  write(lun,*) 

  write(str,lfmt) CV%DoSelection
  write(lun,sfmt)"DoSelection:", trim(adjustl(str))
  write(str,dfmt) CV%SelectRadius
  write(lun,sfmt) "SelectRadius:", trim(adjustl(str))
  write(str,dfmt) CV%SelectXcoord
  write(lun,sfmt) "SelectXcoord:", trim(adjustl(str))
  write(str,dfmt) CV%SelectYcoord
  write(lun,sfmt) "SelectYcoord:", trim(adjustl(str))
  write(str,dfmt) CV%SelectZcoord
  write(lun,sfmt) "SelectZcoord:", trim(adjustl(str))
  write(lun,*)

  write(str,i1fmt) CV%XcoordBndryCond
  write(lun,sfmt) "XcoordBndryCond:", trim(adjustl(str))
  write(str,i1fmt) CV%YcoordBndryCond
  write(lun,sfmt) "YcoordBndryCond:", trim(adjustl(str))
  write(str,i1fmt) CV%ZcoordBndryCond
  write(lun,sfmt) "ZcoordBndryCond:", trim(adjustl(str))
  write(lun,*) 

  write(str,lfmt) CV%RayStats
  write(lun,sfmt) "RayStats:", trim(adjustl(str))
  write(str,dfmt) CV%MaxRayDist
  write(lun,sfmt) "MaxRayDist:", trim(adjustl(str))
  write(lun,*)

  write(str,i1fmt) CV%HealpixNside
  write(lun,sfmt) "HealpixNside:", trim(adjustl(str))
  write(str,i1fmt) CV%NumIterations
  write(lun,sfmt) "NumIterations:", trim(adjustl(str))
  write(str,dfmt) CV%ConvTol
  write(lun,sfmt) "ConvTol:", trim(adjustl(str))
  write(lun,*)

  write(str,dfmt) CV%SpecLoRyd
  write(lun,sfmt) "SpecLoRyd:", trim(adjustl(str))
  write(str,dfmt) CV%SpecHiRyd
  write(lun,sfmt) "SpecHiRyd:", trim(adjustl(str))
  write(str,dfmt) CV%GammaMultiplier
  write(lun,sfmt) "GammaMultiplier:", trim(adjustl(str))
  write(str,i1fmt) CV%FluxDistribution
  write(lun,sfmt) "FluxDistribution:", trim(adjustl(str))
  write(str,dfmt) CV%RhoMultiplier
  write(lun,sfmt) "RhoMultiplier:", trim(adjustl(str))
  write(lun,*)

  write(str,dfmt) CV%NeBackground
  write(lun,sfmt) "NeBackground:", trim(adjustl(str))
  write(str,dfmt) CV%H_mf
  write(lun,sfmt) "H_mf:", trim(adjustl(str))
  write(str,dfmt) CV%He_mf
  write(lun,sfmt) "He_mf:", trim(adjustl(str))
  write(lun,*)



  write(lun,sfmt) "OutputDir:", trim(CV%OutputDir)
  write(lun,sfmt) "OutputFilBase:", trim(CV%OutputFileBase)
  write(lun,sfmt)  "OutputType:", trim(CV%OutputType)
  write(lun,*) 

  write(str,i1fmt) CV%PartPerCell
  write(lun,sfmt)  "PartPerCell: ", trim(adjustl(str))

  write(lun,*) 
  write(lun,*) 

  if (CV%UseIsoTemp) then
     write(lun,*) "***********************************************************"
     write(lun,*) "you are running a constant temperature simulation."
     write(lun,*) "the temperature is fixed at T (K) = ", CV%IsoTemp
     write(lun,*) "***********************************************************"
  end if

  write(lun,*) 
  write(lun,*) 
  write(lun,*)"====================================="
  write(lun,*)"   End URCHIN Configuration Output   "
  write(lun,*)"====================================="
  write(lun,*)
  write(lun,*)

  close(lun)

end subroutine config_info_to_file



!> writes the configuration file information to the snapshot output
!! hdf5 file
!====================================================================
subroutine write_config_hdf5_lun(lun)
  integer(i4b) :: lun
  
  call hdf5_write_attribute( lun, 'Config/Verbosity',         CV%Verbosity )

  if (CV%JustInit) then
     call hdf5_write_attribute( lun, 'Config/JustInit',       1)
  else
     call hdf5_write_attribute( lun, 'Config/JustInit',       0)
  endif
  
  if (CV%Comoving) then
     call hdf5_write_attribute( lun, 'Config/Comoving',       1)
  else
     call hdf5_write_attribute( lun, 'Config/Comoving',       0)
  endif

  call hdf5_write_attribute( lun, 'Config/IntSeed',           CV%IntSeed)
  
  if (CV%UseIsoTemp) then
     call hdf5_write_attribute( lun, 'Config/UseIsoTemp',     1)
  else
     call hdf5_write_attribute( lun, 'Config/UseIsoTemp',     0)
  endif
  call hdf5_write_attribute( lun, 'Config/IsoTemp',           CV%IsoTemp)

  if (CV%AdjustIsmPars) then
     call hdf5_write_attribute( lun, 'Config/AdjustIsmPars',  1)
  else
     call hdf5_write_attribute( lun, 'Config/AdjustIsmPars',  0)
  endif
  call hdf5_write_attribute( lun, 'Config/IsmTemp',           CV%IsmTemp)
  call hdf5_write_attribute( lun, 'Config/IsmHotFrac',        CV%IsmHotFrac)  
  if (CV%IsmAddH2) then
     call hdf5_write_attribute( lun, 'Config/IsmAddH2',  1)
  else
     call hdf5_write_attribute( lun, 'Config/IsmAddH2',  0)
  endif
  call hdf5_write_attribute( lun, 'Config/IsmH2Param1',          CV%IsmH2Param1)
  call hdf5_write_attribute( lun, 'Config/IsmH2Param2',          CV%IsmH2Param2)

  if (CV%AdjustShldPars) then
     call hdf5_write_attribute( lun, 'Config/AdjustShldPars', 1)
  else
     call hdf5_write_attribute( lun, 'Config/AdjustShldPars', 0)
  endif
  call hdf5_write_attribute( lun, 'Config/ShldTemp',          CV%ShldTemp)
  call hdf5_write_attribute( lun, 'Config/ShldTauThresh',     CV%ShldTauThresh)


  call hdf5_write_attribute( lun, 'Config/InputType',         CV%InputType)
  call hdf5_write_attribute( lun, 'Config/EOSvarName',        CV%EOSvarName)
  call hdf5_write_attribute( lun, 'Config/SnapPath',          CV%SnapPath)
  call hdf5_write_attribute( lun, 'Config/ParFileBase',       CV%ParFileBase)  
  call hdf5_write_attribute( lun, 'Config/ParFilesPerSnap',   CV%ParFilesPerSnap)

  call hdf5_write_attribute( lun, 'Config/ItabDir',           CV%ItabDir)
  call hdf5_write_attribute( lun, 'Config/AtomicRatesFile',   CV%AtomicRatesFile)
  call hdf5_write_attribute( lun, 'Config/SlabLookupFile',    CV%SlabLookupFile)
  
  if (CV%Monochromatic) then
     call hdf5_write_attribute( lun, 'Config/Monochromatic',  1)
  else
     call hdf5_write_attribute( lun, 'Config/Monochromatic',  0)
  endif
  call hdf5_write_attribute( lun, 'Config/MonoRydbergs',      CV%MonoRydbergs)
  call hdf5_write_attribute( lun, 'Config/MonoNormMethod',    CV%MonoNormMethod)


  if (CV%DoSelection) then
     call hdf5_write_attribute( lun, 'Config/DoSelection',    1)
  else
     call hdf5_write_attribute( lun, 'Config/DoSelection',    0)
  endif
  call hdf5_write_attribute( lun, 'Config/SelectRadius',      CV%SelectRadius)
  call hdf5_write_attribute( lun, 'Config/SelectXcoord',      CV%SelectXcoord)
  call hdf5_write_attribute( lun, 'Config/SelectYcoord',      CV%SelectYcoord)
  call hdf5_write_attribute( lun, 'Config/SelectZcoord',      CV%SelectZcoord)
  
  call hdf5_write_attribute( lun, 'Config/XcoordBndryCond',   CV%XcoordBndryCond)
  call hdf5_write_attribute( lun, 'Config/YcoordBndryCond',   CV%YcoordBndryCond)
  call hdf5_write_attribute( lun, 'Config/ZcoordBndryCond',   CV%ZcoordBndryCond)

  if (CV%RayStats) then
     call hdf5_write_attribute( lun, 'Config/RayStats',       1)
  else
     call hdf5_write_attribute( lun, 'Config/RayStats',       0)
  endif
  call hdf5_write_attribute( lun, 'Config/MaxRayDist',        CV%MaxRayDist)
  
  call hdf5_write_attribute( lun, 'Config/HealpixNside',      CV%HealpixNside)
  call hdf5_write_attribute( lun, 'Config/NumIterations',     CV%NumIterations)
  call hdf5_write_attribute( lun, 'Config/ConvTol',           CV%ConvTol)

  call hdf5_write_attribute( lun, 'Config/SpecLoRyd      ',   CV%SpecLoRyd)
  call hdf5_write_attribute( lun, 'Config/SpecHiRyd      ',   CV%SpecHiRyd)
  call hdf5_write_attribute( lun, 'Config/GammaMultiplier',   CV%GammaMultiplier)
  call hdf5_write_attribute( lun, 'Config/FluxDistribution',  CV%FluxDistribution)
  call hdf5_write_attribute( lun, 'Config/RhoMultiplier',     CV%RhoMultiplier)  
  call hdf5_write_attribute( lun, 'Config/NeBackground',      CV%NeBackground)
  
  call hdf5_write_attribute( lun, 'Config/H_mf',              CV%H_mf)
  call hdf5_write_attribute( lun, 'Config/He_mf',             CV%He_mf)
    
  call hdf5_write_attribute( lun, 'Config/OutputDir',         CV%OutputDir)
  call hdf5_write_attribute( lun, 'Config/OutputFileBase',    CV%OutputFileBase)
  call hdf5_write_attribute( lun, 'Config/OutputType',        CV%OutputType)

  call hdf5_write_attribute( lun, 'Config/PartPerCell',       CV%PartPerCell)

  
end subroutine write_config_hdf5_lun





end module config_mod
