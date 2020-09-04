!> \file gadget_sphray_header_class.F90

!> \brief Handles Gadget style headers with extra SPHRAY info.  
!!
!!
!< 

module gadget_sphray_header_class
use myf90_mod
use hdf5_wrapper
!use hdf5
use gadget_general_class
use gadget_public_header_class
use gadget_owls_header_class





implicit none
private

public :: gadget_sphray_header_type
public :: gadget_sphray_header_copy_public
public :: gadget_sphray_header_copy_owls
public :: gadget_sphray_header_print_lun
public :: gadget_sphray_header_write_lun
public :: gadget_sphray_header_hdf5_write_lun


!!$interface write_hdf5_header_attrib
!!$   module procedure write_hdf5_header_attrib_int32_s, write_hdf5_header_attrib_int32_v, &
!!$        write_hdf5_header_attrib_int64_s, &
!!$        write_hdf5_header_attrib_real64_s, write_hdf5_header_attrib_real64_v
!!$end interface


!!$integer(i4b), parameter :: header_arrs_num_dims = 1
!!$integer(HSIZE_T), parameter :: header_arrs_sizes(1) = 6


!> Gadget header type
!-----------------------------------------
type gadget_sphray_header_type
#ifdef EAGLE
   integer(i8b) :: npar_file(0:5)  !< number of particles in snapshot file
#else
   integer(i4b) :: npar_file(0:5)  !< number of particles in snapshot file
#endif
   real(r8b)    :: mass(0:5)       !< mass of each particle type if constant
   real(r8b)    :: a               !< scale factor or time
   real(r8b)    :: z               !< redshift
   integer(i4b) :: flag_sfr        !< flag for star formation
   integer(i4b) :: flag_feedback   !< flag for feedback
#ifdef EAGLE
   integer(i8b) :: npar_all(0:5)   !< number of particles in whole snapshot
#else
   integer(i4b) :: npar_all(0:5)   !< number of particles in whole snapshot
#endif
   integer(i4b) :: flag_cooling    !< flag for radiative cooling
   integer(i4b) :: nfiles          !< number of files in a this snapshot
   real(r8b)    :: boxlen          !< box length
   real(r8b)    :: OmegaM          !< omega matter
   real(r8b)    :: OmegaL          !< omega lambda
   real(r8b)    :: h               !< little hubble
   integer(i4b) :: flag_age        !< flag for stellar age
   integer(i4b) :: flag_metals     !< flag for metallicity
#ifdef EAGLE
   integer(i8b) :: npar_hw(0:5)    !< 64 bit part of npar 
#else
   integer(i4b) :: npar_hw(0:5)    !< 64 bit part of npar
#endif 
   integer(i4b) :: flag_entr_ics   !< ICs contain entropy instead of energy

   real(r8b)    :: OmegaB          !< omega baryon 
   integer(i8b) :: rays_traced     !< number of rays traced
   integer(i4b) :: flag_Hmf        !< output has Hydorgen mass fraction? 
   integer(i4b) :: flag_Hemf       !< output has Helium mass fraction? 
   integer(i4b) :: flag_helium     !< output has xHeI and xHeII?
   integer(i4b) :: flag_gammaHI    !< output has HI photoionization rate? 
   integer(i4b) :: flag_cloudy     !< output has Cloudy xHIeq? 
   integer(i4b) :: flag_eos        !< output has EOS info?
   integer(i4b) :: flag_incsfr     !< output has SFR info?
   real(r8b)    :: time_gyr        !< time since BB from cosmo variables [Gyr]
   integer(i4b) :: unused(2)       !< spacer
end type gadget_sphray_header_type


contains





!> writes a gadget header to an open file
!--------------------------------------------------------------
subroutine gadget_sphray_header_write_lun(this, lun)
  type(gadget_sphray_header_type), intent(in) :: this
  integer(i4b), intent(in) :: lun

  write(lun) this%npar_file(:), this%mass(:), this%a, this%z, &              
       this%flag_sfr, this%flag_feedback, this%npar_all(:), &    
       this%flag_cooling, this%nfiles, this%boxlen, this%OmegaM, &         
       this%OmegaL, this%h, this%flag_age, this%flag_metals, &    
       this%npar_hw(:), this%flag_entr_ics, this%OmegaB, &
       this%rays_traced, this%flag_Hmf, this%flag_Hemf, this%flag_helium, &
       this%flag_gammaHI, this%flag_cloudy, this%flag_eos, this%flag_incsfr, &
       this%time_gyr, this%unused(:)      

end subroutine gadget_sphray_header_write_lun









!> writes a gadget header to an open HDF5 file
!--------------------------------------------------------------
subroutine gadget_sphray_header_hdf5_write_lun(this, file_id)
  type(gadget_sphray_header_type), intent(in) :: this
  integer, intent(in) :: file_id

!!$  integer(HID_T) :: group_id
!!$  integer(i4b) :: hdferr


!!$  ! Create a group named "Header" in the file.
!!$  !-------------------------------------------------
!!$  call h5gcreate_f(file_id, "/Header", group_id, hdferr)
!!$
!!$  ! write attributes
!!$  !-------------------------------------------------
!!$  call write_hdf5_header_attrib( file_id, group_id, 'NumPart_ThisFile', this%npar_file )
!!$  call write_hdf5_header_attrib( file_id, group_id, 'NumPart_Total', this%npar_all )
!!$  call write_hdf5_header_attrib( file_id, group_id, 'NumPart_HighWord', this%npar_hw )
!!$  call write_hdf5_header_attrib( file_id, group_id, 'MassTable', this%mass )
!!$  call write_hdf5_header_attrib( file_id, group_id, 'ExpansionFactor', this%a )
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Time_GYR', this%time_gyr )
!!$
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Redshift', this%z )
!!$  call write_hdf5_header_attrib( file_id, group_id, 'BoxSize', this%boxlen )
!!$  call write_hdf5_header_attrib( file_id, group_id, 'NumFilesPerSnapshot', this%nfiles)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Omega0', this%OmegaM)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'OmegaBaryon', this%OmegaB)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'OmegaLambda', this%OmegaL)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'HubbleParam', this%h)
!!$
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_Sfr', this%flag_sfr)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_Cooling', this%flag_cooling)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_StellarAge', this%flag_age)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_Metals', this%flag_metals)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_Feedback', this%flag_feedback)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_Entropy_In_ICs', this%flag_entr_ics)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'RaysTraced', this%rays_traced)
!!$
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_Hmf', this%flag_Hmf)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_Hemf', this%flag_Hemf)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_Helium', this%flag_helium)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_GammaHI', this%flag_gammaHI)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_Cloudy', this%flag_cloudy)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_EoS', this%flag_eos)
!!$  call write_hdf5_header_attrib( file_id, group_id, 'Flag_IncSFR', this%flag_incsfr)
!!$
!!$
!!$  CALL h5gclose_f( group_id, hdferr )
  

  call hdf5_write_attribute(file_id,'Header/NumPart_ThisFile',this%npar_file)
  call hdf5_write_attribute(file_id,'Header/NumPart_Total',this%npar_all)
  call hdf5_write_attribute(file_id,'Header/NumPart_Total_HighWord',this%npar_hw)
  call hdf5_write_attribute(file_id,'Header/MassTable',this%mass)
  call hdf5_write_attribute(file_id,'Header/ExpansionFactor',this%a)
  call hdf5_write_attribute(file_id,'Header/Time_GYR',this%time_gyr)

  call hdf5_write_attribute(file_id,'Header/Redshift',this%z)
  call hdf5_write_attribute(file_id,'Header/BoxSize',this%boxlen)
  call hdf5_write_attribute(file_id,'Header/NumFilesPerSnapshot',this%nfiles)
  call hdf5_write_attribute(file_id,'Header/Omega0',this%OmegaM)
  call hdf5_write_attribute(file_id,'Header/OmegaBaryon',this%OmegaB)
  call hdf5_write_attribute(file_id,'Header/OmegaLambda',this%OmegaL)
  call hdf5_write_attribute(file_id,'Header/HubbleParam',this%h)

  call hdf5_write_attribute(file_id,'Header/Flag_Sfr',this%flag_sfr)
  call hdf5_write_attribute(file_id,'Header/Flag_Cooling',this%flag_cooling)
  call hdf5_write_attribute(file_id,'Header/Flag_StellarAge',this%flag_age)
  call hdf5_write_attribute(file_id,'Header/Flag_Metals',this%flag_metals)
  call hdf5_write_attribute(file_id,'Header/Flag_Feedback',this%flag_feedback)
  call hdf5_write_attribute(file_id,'Header/Flag_Entropy_In_ICs',this%flag_entr_ics)
  call hdf5_write_attribute(file_id,'Header/RaysTraced',this%rays_traced)

  call hdf5_write_attribute(file_id,'Header/Flag_Hmf',this%flag_Hmf)
  call hdf5_write_attribute(file_id,'Header/Flag_Hemf',this%flag_Hemf)
  call hdf5_write_attribute(file_id,'Header/Flag_Helium',this%flag_helium)
  call hdf5_write_attribute(file_id,'Header/Flag_GammaHI',this%flag_gammaHI)
  call hdf5_write_attribute(file_id,'Header/Flag_Cloudy',this%flag_cloudy)
  call hdf5_write_attribute(file_id,'Header/Flag_EoS',this%flag_eos)
  call hdf5_write_attribute(file_id,'Header/Flag_IncSFR',this%flag_incsfr)


end subroutine gadget_sphray_header_hdf5_write_lun


!> copies an OWLS/GIMIC header into the sphray style header
!--------------------------------------------------------------
subroutine gadget_sphray_header_copy_owls( this, owlshead )
  type(gadget_sphray_header_type) :: this
  type(gadget_owls_header_type) :: owlshead

  this%npar_file(0:5)   = owlshead%npar_file(0:5)   
  this%mass(0:5)        = owlshead%mass(0:5)        
  this%a                = owlshead%a                
  this%z                = owlshead%z                
  this%flag_sfr         = owlshead%flag_sfr         
  this%flag_feedback    = owlshead%flag_feedback    
  this%npar_all(0:5)    = owlshead%npar_all(0:5)    
  this%flag_cooling     = owlshead%flag_cooling     
  this%nfiles           = owlshead%nfiles           
  this%boxlen           = owlshead%boxlen           
  this%OmegaM           = owlshead%OmegaM           
  this%OmegaL           = owlshead%OmegaL           
  this%h                = owlshead%h                
  this%flag_age         = owlshead%flag_age         
  this%flag_metals      = owlshead%flag_metals      
  this%npar_hw(0:5)     = owlshead%npar_hw(0:5)     

  this%flag_entr_ics    = 0

  this%OmegaB           = owlshead%OmegaB
  this%time_gyr         = owlshead%time_gyr



end subroutine gadget_sphray_header_copy_owls


!> copies a Public header into the sphray style header
!--------------------------------------------------------------
subroutine gadget_sphray_header_copy_public( this, pubhead )
  type(gadget_sphray_header_type) :: this
  type(gadget_public_header_type) :: pubhead

  this%npar_file(0:5)   = pubhead%npar_file(0:5)   
  this%mass(0:5)        = pubhead%mass(0:5)        
  this%a                = pubhead%a                
  this%z                = pubhead%z                
  this%flag_sfr         = pubhead%flag_sfr         
  this%flag_feedback    = pubhead%flag_feedback    
  this%npar_all(0:5)    = pubhead%npar_all(0:5)    
  this%flag_cooling     = pubhead%flag_cooling     
  this%nfiles           = pubhead%nfiles           
  this%boxlen           = pubhead%boxlen           
  this%OmegaM           = pubhead%OmegaM           
  this%OmegaL           = pubhead%OmegaL           
  this%h                = pubhead%h                
  this%flag_age         = pubhead%flag_age         
  this%flag_metals      = pubhead%flag_metals      
  this%npar_hw(0:5)     = pubhead%npar_hw(0:5)     
  this%flag_entr_ics    = pubhead%flag_entr_ics    

end subroutine gadget_sphray_header_copy_public




!> formatted print of header to lun (including standard out)
!---------------------------------------------------------------
subroutine gadget_sphray_header_print_lun(this, lun)
  type(gadget_sphray_header_type), intent(in) :: this 
  integer(i4b), intent(in) :: lun
  integer(i8b) :: i

  character(clen) :: n1,n2,n3
  character(clen) :: clmn_fmt
  character(clen) :: type_fmt
  character(clen) :: totl_fmt
  character(clen) :: scal_fmt
  character(clen) :: flag_fmt
  character(clen) :: star_fmt
  character(clen) :: line_fmt

  clmn_fmt = "(T17,A,T32,A,T47,A)"
  type_fmt = "(A,I1,A,A,A,T17,A,T32,A,T47,A)"
  totl_fmt = "(A,T17,A,T47,A)"
  scal_fmt = "(A,A,T25,A,A,T45,A,A)"
  flag_fmt = "(A,6(A,I1))"
  star_fmt = "(78('='))"
  line_fmt = "(78('-'))"

  write(lun,*)
  write(lun,star_fmt)
  write(lun,clmn_fmt) "n_all", "mass", "n_file"
  write(lun,line_fmt)
  do i = 0,5
     write(n1,'(I20)') this%npar_all(i)
     write(n2,'(E12.5)') this%mass(i)
     write(n3,'(I20)') this%npar_file(i)
     write(lun,type_fmt) "type",i,"(",gadget_ptype_names(i),")",&
          trim(adjustl(n1)), trim(adjustl(n2)), trim(adjustl(n3))
  end do
  write(lun,line_fmt)
  write(n1,'(I20)') sum(this%npar_all)
  write(n2,'(I20)') sum(this%npar_file)
  write(lun,totl_fmt) "total:", trim(adjustl(n1)), trim(adjustl(n2)) 
  write(lun,*)



  write(n1,'(F12.5)') this%a
  write(n2,'(F12.5)') this%z
  write(n3,'(F12.5)') this%h

  write(lun,scal_fmt) "a = ", trim(adjustl(n1)), &
       "z = ", trim(adjustl(n2)), &
       "h = ", trim(adjustl(n3))

  write(n1,'(F12.5)') this%OmegaM
  write(n2,'(F12.5)') this%OmegaL
  write(n3,'(F12.5)') this%boxlen

  write(lun,scal_fmt) "OmegaM = ", trim(adjustl(n1)), &
       "OmegaL = ", trim(adjustl(n2)), &
       "BoxSize = ", trim(adjustl(n3))
  write(lun,*)

  write(n1,'(I20)') this%nfiles
  write(lun,'(A,A)') "nfiles = ", trim(adjustl(n1))

  write(lun,flag_fmt) "/flags/ ",&
       "  sfr=",this%flag_sfr, &
       ", feedback=",this%flag_feedback, &
       ", cooling=",this%flag_cooling, &
       ", age=",this%flag_age, &
       ", metals=",this%flag_metals, &
       ", entr_ics=", this%flag_entr_ics
  write(lun,*) 
  write(n1,'(F12.5)') this%time_gyr
  write(lun,'(A,A)') "  time[Gyr]: ", trim(n1)

  write(lun,star_fmt)

end subroutine gadget_sphray_header_print_lun





!===================================================================================
!   HDF5 output routines
!===================================================================================



!!$! HDF5 write attribute integer 4B - scalar
!!$!-------------------------------------------------------
!!$subroutine write_hdf5_header_attrib_int32_s( file_id, group_id, attr_name, attr_val )
!!$  integer(HID_T), intent(in) :: file_id
!!$  integer(HID_T), intent(in) :: group_id
!!$  character(*), intent(in) :: attr_name
!!$  integer(i4b), intent(in) :: attr_val
!!$
!!$  integer(HID_T) :: dspace_id, attrib_id   
!!$  integer(i4b) :: hdferr
!!$
!!$  integer(HID_T) :: I4_DTYPE
!!$
!!$  I4_DTYPE = H5T_NATIVE_INTEGER
!!$
!!$  call h5screate_f(H5S_SCALAR_F, dspace_id, hdferr)
!!$  call h5acreate_f( group_id, trim(attr_name), I4_DTYPE, dspace_id, attrib_id, hdferr )
!!$  call h5awrite_f( attrib_id, I4_DTYPE, attr_val, header_arrs_sizes, hdferr )
!!$  call h5aclose_f( attrib_id, hdferr )
!!$  
!!$end subroutine write_hdf5_header_attrib_int32_s



!!$! HDF5 write attribute integer 8B - scalar
!!$!-------------------------------------------------------
!!$subroutine write_hdf5_header_attrib_int64_s( file_id, group_id, attr_name, attr_val )
!!$  integer(HID_T), intent(in) :: file_id
!!$  integer(HID_T), intent(in) :: group_id
!!$  character(*), intent(in) :: attr_name
!!$  integer(i8b), intent(in) :: attr_val
!!$
!!$  integer(HID_T) :: dspace_id, attrib_id   
!!$  integer(i4b) :: hdferr
!!$
!!$  integer(HID_T) :: I8_DTYPE
!!$
!!$  I8_DTYPE = H5T_STD_I64LE
!!$
!!$  call h5screate_f(H5S_SCALAR_F, dspace_id, hdferr)
!!$  call h5acreate_f( group_id, trim(attr_name), I8_DTYPE, dspace_id, attrib_id, hdferr )
!!$  call h5awrite_f( attrib_id, I8_DTYPE, attr_val, header_arrs_sizes, hdferr )
!!$  call h5aclose_f( attrib_id, hdferr )
!!$  
!!$end subroutine write_hdf5_header_attrib_int64_s








!!$! HDF5 write attribute float 8B - scalar
!!$!-------------------------------------------------------
!!$subroutine write_hdf5_header_attrib_real64_s( file_id, group_id, attr_name, attr_val )
!!$  integer(HID_T), intent(in) :: file_id
!!$  integer(HID_T), intent(in) :: group_id
!!$  character(*), intent(in) :: attr_name
!!$  real(r8b), intent(in) :: attr_val
!!$
!!$  integer(HID_T) :: dspace_id, attrib_id   
!!$  integer(i4b) :: hdferr
!!$
!!$  integer(HID_T) :: R8_DTYPE
!!$
!!$  R8_DTYPE = H5T_NATIVE_DOUBLE
!!$
!!$  call h5screate_f( H5S_SCALAR_F, dspace_id, hdferr )
!!$  call h5acreate_f( group_id, trim(attr_name), R8_DTYPE, dspace_id, attrib_id, hdferr )
!!$  call h5awrite_f( attrib_id, R8_DTYPE, attr_val, header_arrs_sizes, hdferr )
!!$  call h5aclose_f( attrib_id, hdferr )
!!$  
!!$end subroutine write_hdf5_header_attrib_real64_s
!!$
!!$
!!$! HDF5 write attribute integer 4B - vector
!!$!-------------------------------------------------------
!!$subroutine write_hdf5_header_attrib_int32_v( file_id, group_id, attr_name, attr_val )
!!$  integer(HID_T), intent(in) :: file_id
!!$  integer(HID_T), intent(in) :: group_id
!!$  character(*), intent(in) :: attr_name
!!$  integer(i4b), intent(in) :: attr_val(header_arrs_num_dims)
!!$
!!$  integer(HID_T) :: dspace_id, attrib_id   
!!$  integer(i4b) :: hdferr
!!$
!!$  integer(HID_T) :: I4_DTYPE
!!$
!!$  I4_DTYPE = H5T_NATIVE_INTEGER
!!$
!!$  call h5screate_simple_f( header_arrs_num_dims, header_arrs_sizes, dspace_id, hdferr )
!!$  call h5acreate_f( group_id, trim(attr_name), I4_DTYPE, dspace_id, attrib_id, hdferr )
!!$  call h5awrite_f( attrib_id, I4_DTYPE, attr_val, header_arrs_sizes, hdferr )
!!$  call h5aclose_f( attrib_id, hdferr )
!!$  
!!$end subroutine write_hdf5_header_attrib_int32_v
!!$
!!$
!!$! HDF5 write attribute float 8B - vector
!!$!-------------------------------------------------------
!!$subroutine write_hdf5_header_attrib_real64_v( file_id, group_id, attr_name, attr_val )
!!$  integer(HID_T), intent(in) :: file_id
!!$  integer(HID_T), intent(in) :: group_id
!!$  character(*), intent(in) :: attr_name
!!$  real(r8b), intent(in) :: attr_val(header_arrs_num_dims)
!!$
!!$  integer(HID_T) :: dspace_id, attrib_id   
!!$  integer(i4b) :: hdferr
!!$
!!$  integer(HID_T) :: R8_DTYPE
!!$
!!$  R8_DTYPE = H5T_NATIVE_DOUBLE
!!$
!!$  call h5screate_simple_f( header_arrs_num_dims, header_arrs_sizes, dspace_id, hdferr )
!!$  call h5acreate_f( group_id, trim(attr_name), R8_DTYPE, dspace_id, attrib_id, hdferr )
!!$  call h5awrite_f( attrib_id, R8_DTYPE, attr_val, header_arrs_sizes, hdferr )
!!$  call h5aclose_f( attrib_id, hdferr )
!!$  
!!$end subroutine write_hdf5_header_attrib_real64_v















end module gadget_sphray_header_class
