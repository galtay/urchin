!> \file output_gadget_hdf5.F90

!> \brief The module that handles output in 
!! the style of the hdf5 version of Gadget
!<
module output_gadget_hdf5_mod
! libs
use myf90_mod
use hdf5_wrapper
! types
use particle_system_mod, only: particle_system_type
use gadget_sphray_header_class, only: gadget_sphray_header_type
use gadget_general_class, only: gadget_constants_type
use gadget_general_class, only: gadget_units_type
! routines
use particle_returns_mod, only: return_H_mf
use global_mod, only: recall_saved_ghead
use gadget_sphray_header_class, only: gadget_sphray_header_hdf5_write_lun
use gadget_general_class, only: gadget_units_write_lun
use gadget_general_class, only: gadget_data_attributes_set
use gadget_general_class, only: gadget_data_attributes_write_lun
use gadget_general_class, only: gadget_constants_write_lun
use gadget_general_class, only: gadget_owls_parameters_write_lun
use config_mod, only: write_config_hdf5_lun
! variables
use config_mod, only: CV
use global_mod, only: GV
use global_mod, only: global_owls_parameters
use gadget_general_class, only: gattrs
implicit none
private


public :: output_snap_gadget_hdf5


contains


subroutine output_snap_gadget_hdf5(psys, NumFiles)

  character(clen), parameter :: myname="output_snap_gadget_hdf5"
  logical, parameter :: crash = .true.
  integer, parameter :: verb = 2

  type(particle_system_type), intent(in) :: psys
  integer(i4b) :: NumFiles  !< number of files to put output in

  type(gadget_constants_type) :: gconst

  integer(i4b) :: ifile

  integer(i8b) :: Nwrote   
  integer(i8b) :: ngas1
  integer(i8b) :: ngood1
  integer(i8b) :: ngood
  integer(i8b) :: igood
  integer(i8b) :: indx
  integer(i8b) :: i

  character(3) :: label
  character(4) :: ext
  character(clen) :: filename
  
  type(gadget_sphray_header_type) :: ghead
  type(gadget_units_type) :: gunits

  real(r4b), allocatable :: rblock3(:,:)
  real(r4b), allocatable :: rblock(:)
#ifdef LONGIDs
  integer(i8b), allocatable :: iblock(:)  !< particle id
#else
  integer(i4b), allocatable :: iblock(:)  !< particle id
#endif

  character(200) :: RunLabelDummy
  character(clen) :: tag
  character(clen) :: group_name
  integer(i4b) :: file_id





  RunLabelDummy = 'Urchin Output'

  group_name = '/PartType0/'

  write(label,'(I3.3)') GV%OutputIndx

  Nwrote = 0

  ! count good particles to output (usually all)
  !----------------------------------------------
  ngood = count( psys%par(:)%mask )



  over_files: do ifile = 0, NumFiles-1

     call recall_saved_ghead( ifile, ghead )


     ! adjust if we have a selection
     !-------------------------------
     if (CV%DoSelection) then
        ngas1 = size(psys%par)
        ngood1 = ngood
        ghead%nfiles = 1
        ghead%npar_file    = 0
        ghead%npar_file(0) = ngood1
        ghead%npar_all     = 0
        ghead%npar_all(0)  = ngood
     else
        ngas1 = ghead%npar_file(0)
        ngood1 = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(indx)%mask) then
              ngood1 = ngood1 + 1
           endif
        end do
        ghead%npar_file    = 0
        ghead%npar_file(0) = ngood1
        ghead%npar_all     = 0
        ghead%npar_all(0)  = ngood
     endif
     

     
     
     ! form output file name
     !--------------------------------     
     write(ext,'(I3)') ifile
     
     filename = trim(CV%OutputDir) // "/" // trim(CV%OutputFileBase) // "_" // label 
     if (NumFiles > 1) then
        filename = trim(filename) // "." // trim(adjustl(ext))
     end if
     filename = trim(filename) // ".hdf5"
     
     call mywrite("writing snapshot state to "//trim(filename), verb)
     call mywrite('',verb) 


     call hdf5_create_file( file_id, filename )
     

     !================================================================
     ! write Header
     !================================================================    
     call hdf5_create_group( file_id, '/Header' )
     call gadget_sphray_header_hdf5_write_lun( ghead, file_id )
     if (trim(CV%InputType) == 'owls') then 
        call hdf5_write_attribute(file_id,'/Header/RunLabel', RunLabelDummy)
     endif
     

     !================================================================
     ! write Units
     !================================================================
     call hdf5_create_group( file_id, 'Units/' )
     gunits%cgs_length   = GV%cgs_len
     gunits%cgs_mass     = GV%cgs_mass
     gunits%cgs_velocity = GV%cgs_vel
     gunits%cgs_density  = GV%cgs_rho
     gunits%cgs_energy   = GV%cgs_enrg
     gunits%cgs_pressure = GV%cgs_prs
     gunits%cgs_time     = GV%cgs_time
     call gadget_units_write_lun( gunits, file_id )          

     call gadget_data_attributes_set(gunits)


     !================================================================
     ! write Constants
     !================================================================
     call hdf5_create_group( file_id, 'Constants/' )
     call gadget_constants_write_lun( gconst, file_id )


     !================================================================
     ! write OWLS/GIMIC Parameters 
     !================================================================
     if (trim(CV%InputType) == 'owls') then 
        call hdf5_create_group( file_id, 'Parameters/' )
        call hdf5_create_group( file_id, 'Parameters/ChemicalElements' )
        call hdf5_create_group( file_id, 'Parameters/NumericalParameters' )
        call hdf5_create_group( file_id, 'Parameters/StellarEvolutionParameters' )
        call hdf5_create_group( file_id, 'Parameters/WindParameters' )
        call gadget_owls_parameters_write_lun( global_owls_parameters, file_id )
     endif


     !================================================================
     ! write config parameters
     !================================================================
     call hdf5_create_group( file_id, 'Config/' )
     call write_config_hdf5_lun(file_id)

     if (ngas1 == 0) then
        call hdf5_close_file(file_id)
        cycle
     endif




     !================================================================
     ! write data
     !================================================================



     ! Only do position and velocity for selections
     !----------------------------------------------------------------
     if (CV%DoSelection) then

        allocate( rblock3(3,ngood1) )

        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(indx)%mask) then
              igood = igood + 1
              rblock3(1,igood) = psys%par( indx )%pos(1)
              rblock3(2,igood) = psys%par( indx )%pos(2)
              rblock3(3,igood) = psys%par( indx )%pos(3)
           endif
        end do
        tag = trim(group_name)//'Coordinates'
        call hdf5_write_data( file_id, tag, rblock3 )
        call gadget_data_attributes_write_lun( gattrs%pos, file_id, group_name, 'Coordinates' )


#ifdef incVel
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock3(1,igood) = psys%par( indx )%vel(1)
              rblock3(2,igood) = psys%par( indx )%vel(2)
              rblock3(3,igood) = psys%par( indx )%vel(3)
           endif
        end do

        tag = trim(group_name)//'Velocity'
        call hdf5_write_data( file_id, tag, rblock3 )
        call gadget_data_attributes_write_lun( gattrs%vel, file_id, group_name, 'Velocity' )

#else

#endif

        deallocate( rblock3 )

     endif


     ! Always do particle IDs just as a sanity check
     !----------------------------------------------------------------
     allocate( iblock(ngood1) )
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           iblock(igood) = psys%par( indx )%id
        endif
     end do
     tag = trim(group_name)//'ParticleIDs'
     call hdf5_write_data( file_id, tag, iblock )
     call gadget_data_attributes_write_lun( gattrs%id, file_id, group_name, 'ParticleIDs' )
     deallocate(iblock)




     allocate(rblock(ngood1))


     ! Only do mass and density for selections
     !----------------------------------------------------------------
     if (CV%DoSelection) then
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%mass
           endif
        end do
        tag = trim(group_name)//'Mass'
        call hdf5_write_data( file_id, tag, rblock )
        call gadget_data_attributes_write_lun( gattrs%mass, file_id, group_name, 'Mass' )

        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%rho
           endif
        end do
        tag = trim(group_name)//'Density'
        call hdf5_write_data( file_id, tag, rblock )
        call gadget_data_attributes_write_lun( gattrs%rho, file_id, group_name, 'Density' )
     endif


     ! Always do electron fraction
     !----------------------------------------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%ye
        endif
     end do
     tag = trim(group_name)//'ElectronFraction'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%ye, file_id, group_name, 'ElectronFraction' )
     

     ! Always do xHI
     !----------------------------------------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%xHI
        endif
     end do
     tag = trim(group_name)//'HydrogenOneFraction'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%xHI, file_id, group_name, 'HydrogenOneFraction' )


     ! Only do Smoothing Length for selections
     !----------------------------------------------------------------
     if (CV%DoSelection) then

        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%hsml
           endif
        end do
        tag = trim(group_name)//'SmoothingLength'
        call hdf5_write_data( file_id, tag, rblock )
        call gadget_data_attributes_write_lun( gattrs%hsml, file_id, group_name, 'SmoothingLength' ) 
     
     endif


     ! Always do temperature as we may change it
     !----------------------------------------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%T
        endif
     end do
     tag = trim(group_name)//'Temperature'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%T, file_id, group_name, 'Temperature' )


     ! Only do Element Mass Fractions for selections (we'll put in nH) 
     !----------------------------------------------------------------
     if (CV%DoSelection) then

#ifdef incHmf
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%Hmf
           endif
        end do
        tag = trim(group_name)//'ElementAbundance/Hydrogen'
        call hdf5_write_data( file_id, tag, rblock )
        call gadget_data_attributes_write_lun( gattrs%Hmf, file_id, group_name, 'ElementAbundance/Hydrogen' )
#endif


#ifdef incHemf
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%Hemf
           endif
        end do
        tag = trim(group_name)//'ElementAbundance/Helium'
        call hdf5_write_data( file_id, tag, rblock )
        call gadget_data_attributes_write_lun( gattrs%Hemf, file_id, group_name, 'ElementAbundance/Helium' )
#endif

     endif


     ! Always put Helium ionization states 
     !----------------------------------------------------------------
#ifdef incHe
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%xHeI
        endif
     end do
     tag = trim(group_name)//'HeliumOneFraction'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%xHeI, file_id, group_name, 'HeliumOneFraction' )

     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%xHeII
        endif
     end do
     tag = trim(group_name)//'HeliumTwoFraction'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%xHeII, file_id, group_name, 'HeliumTwoFraction' )
#endif


     ! Always put photoionization rate
     !----------------------------------------------------------------
#ifdef outGammaHI
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%gammaHI
        endif
     end do
     tag = trim(group_name)//'HydrogenOneGamma'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%gammaHI, file_id, group_name, 'HydrogenOneGamma' )
#endif


     ! Always put optically thin xHI
     !----------------------------------------------------------------
#ifdef incCloudy
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%xHI_cloudy
        endif
     end do
     tag = trim(group_name)//'HydrogenOneFraction_Cloudy'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%xHI_cloudy, file_id, group_name, 'HydrogenOneFraction_Cloudy' )
#endif


     ! Only put SF variable for selections
     !----------------------------------------------------------------
     if (CV%DoSelection) then
#ifdef incEOS
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%eos
           endif
        end do
        tag = trim(group_name)//trim(CV%EOSvarName)
        call hdf5_write_data( file_id, tag, rblock )
        call gadget_data_attributes_write_lun( gattrs%eos, file_id, group_name, trim(CV%EOSvarName) ) 
#endif


     ! Only put SFR for selections
     !----------------------------------------------------------------
#ifdef incSFR
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%sfr
        endif
     end do
     tag = trim(group_name)//'StarFormationRate'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%sfr, file_id, group_name, 'StarFormationRate' )
#endif

     endif


     ! Always put molecular Hydrogen Mass fraction
     !----------------------------------------------------------------
#ifdef incH2
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%fH2
        endif
     end do
     tag = trim(group_name)//'MolecularHydrogenMassFraction'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%fh2, file_id, group_name, 'MolecularHydrogenMassFraction' )
#endif


     ! Always put nH
     !----------------------------------------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%nH
        endif
     end do
     tag = trim(group_name)//'HydrogenNumberDensity'
     call hdf5_write_data( file_id, tag, rblock )
     call gadget_data_attributes_write_lun( gattrs%nH, file_id, group_name, 'HydrogenNumberDensity' )



     deallocate( rblock )
     call hdf5_close_file(file_id)
     Nwrote = Nwrote + ngas1


  end do over_files

  



end subroutine output_snap_gadget_hdf5


end module output_gadget_hdf5_mod
