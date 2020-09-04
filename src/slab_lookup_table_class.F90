!> \file slab_lookup_table_class.F90

!> \brief For constructing lookup tables of analytic solutions.  In this 
!! case, plane parallel radiation with an HM01 spectrum incident on 
!! a constant density slab. 
!! 
!<


module slab_lookup_table_class
use timers
use omp_lib
use myf90_mod
use hdf5_wrapper
use ion_table_class, only: read_ion_table_file
use ion_table_class, only: ion_table_type
use ion_table_class, only: mini_spec_type
use ion_table_class, only: set_mini_spectrum
use ion_table_class, only: gammaHI_from_mini_spec_shield_plane
use ion_table_class, only: gammaHI_from_mini_spec_thin
use ion_table_class, only: HI_photo_cs_verner
use analytic_ionization_solutions_mod, only: analytic_xHIeq
use hui_gnedin_atomic_rates_mod, only: Hui_HI_col_ion
use hui_gnedin_atomic_rates_mod, only: Hui_HII_recombA
use config_mod, only: CV

implicit none
private



public :: create_slab_lookup_table
public :: load_slab_lookup_table
public :: interpolate_slab_table
public :: create_and_solve_slab

public :: return_min_table_vals
public :: return_max_table_vals


type table_dimension_type
   real(r8b) :: logmin
   real(r8b) :: logmax
   real(r8b) :: min
   real(r8b) :: max
   integer(i8b) :: nbins
   real(r4b), allocatable :: values(:)
end type table_dimension_type


type slab_table_type
   type(table_dimension_type) :: dim(0:3)    ! tauH, K=GHI/nH, tauHI, T
   real(r4b), allocatable :: table(:,:,:,:)  ! the table [tauH, K, tauHI, T]
end type slab_table_type

integer(i4b), parameter :: ndims = 4

character(10), parameter :: table_name = 'log_GHI_eq'

character(10), parameter :: dim_names(0:ndims-1) = (/&
     'log_tauH  ', &
     'log_K     ', &
     'log_tauHI ', &
     'log_T     ' /)


type(slab_table_type) :: stab
type(ion_table_type) :: itab
real(r8b), allocatable :: xHI_eq_ci(:)

real(r8b) :: sigma_th 

real(r8b), parameter :: kpc_in_cm = 3.0856780d21







contains
  
  
! reads and loads slab lookup table from file
!======================================================
subroutine load_slab_lookup_table(file)
  character(clen), parameter :: myname = "load_slab_lookup_table"
  integer, parameter :: verb=1
  logical, parameter :: crash = .true.
  character(clen) :: str
 
  character(clen) :: file
  integer(i4b) :: fh
  integer(i4b) :: rank, dims(3),i
  logical :: fthere

  call mywrite('', verb)
  write(str,'(A,A)') 'using slab lookup file: ', trim(file)
  call mywrite(str, verb) 
  call mywrite('', verb)
  
  inquire(file=file, exist=fthere)
  if (.not. fthere) call myerr('cant find slab lookup file',myname,crash)
  
  
  call hdf5_open_file(fh, file, readonly=.true.)

  
  do i = 0, ndims-1
     call hdf5_get_dimensions(fh, trim(dim_names(i)), rank, dims)
     allocate( stab%dim(i)%values( 0:dims(1)-1 ) )
     call hdf5_read_data(fh, trim(dim_names(i)), stab%dim(i)%values )
     stab%dim(i)%nbins = size(stab%dim(i)%values)
     
     stab%dim(i)%logmin = stab%dim(i)%values(0)
     stab%dim(i)%logmax = stab%dim(i)%values(stab%dim(i)%nbins-1)
     stab%dim(i)%min    = 10.0d0**stab%dim(i)%logmin
     stab%dim(i)%max    = 10.0d0**stab%dim(i)%logmax     
  end do
     
  allocate( stab%table( 0:stab%dim(0)%nbins-1, 0:stab%dim(1)%nbins-1, &
                        0:stab%dim(2)%nbins-1, 0:stab%dim(3)%nbins-1 ) )

  call hdf5_read_data(fh, table_name, stab%table )

  write(*,'(A)') '  slab table dims, log min/max: '
  write(*,'(T10,A,I5,2X,2F9.3)') ' tau_H:    ', &
       stab%dim(0)%nbins, stab%dim(0)%logmin, stab%dim(0)%logmax
  write(*,'(T10,A,I5,2X,2F9.3)') ' G_HI/n_H: ', &
       stab%dim(1)%nbins, stab%dim(1)%logmin, stab%dim(1)%logmax
  write(*,'(T10,A,I5,2X,2F9.3)') ' tau_HI:   ', &
       stab%dim(2)%nbins, stab%dim(2)%logmin, stab%dim(2)%logmax
  write(*,'(T10,A,I5,2X,2F9.3)') ' T:        ', &
       stab%dim(3)%nbins, stab%dim(3)%logmin, stab%dim(3)%logmax
  write(*,*) 

end subroutine load_slab_lookup_table







! solves for table values as a function of H optical 
! depth = NH * sigma_th given a set of input parameters
!======================================================
function create_and_solve_slab( nH, GHI_0, Lslab, Tslab, Nlayers, &
     mini_spec, tauH_bins, tauH_arr ) result(log_table_eq)

  real(r8b) :: nH
  real(r8b) :: GHI_0
  real(r8b) :: Lslab
  real(r8b) :: Tslab
  integer(i8b) :: Nlayers 
  type(mini_spec_type) :: mini_spec
  integer(i8b) :: tauH_bins
  real(r8b) :: tauH_arr(tauH_bins)
  real(r8b) :: log_table_eq(tauH_bins)

  real(r8b) :: GHI_thin
  real(r8b) :: tauHI_eff
  real(r8b) :: dz

  real(r8b) :: dtauH
  real(r8b) :: tauH_prior, tauH, tauH_table
  real(r8b) :: tauH_lo, tauH_hi
  real(r8b) :: frac

  real(r8b) :: dtauHI
  real(r8b) :: tauHI

  real(r8b) :: xHI_prior, xHI
  real(r8b) :: GHI_prior, GHI


  integer(i8b) :: i, j, itauH

  real(r8b), parameter :: y = 0.0d0
  


  ! check input GHI is consistent with mini_spec
  !-------------------------------------------------
  sigma_th = HI_photo_cs_verner( 1.0d0 )        
  GHI_thin = gammaHI_from_mini_spec_shield_plane( mini_spec, tauHI_th=0.0d0) 
  if ( (GHI_thin-GHI_0)/GHI_0 > 1.0d-5) then
     write(*,*) 'GHI_thin = ', GHI_thin
     write(*,*) 'GHI_0    = ', GHI_0
     stop 'GHI mismatch'
  endif

  dz    = Lslab / Nlayers
  dtauH = dz * kpc_in_cm * nH * sigma_th
 
  tauH   = 0.0d0
  tauHI  = 0.0d0
  GHI    = gammaHI_from_mini_spec_shield_plane( mini_spec, tauHI_th=tauHI ) 
  tauHI_eff = -log( GHI/GHI_thin )
  xHI    = analytic_xHIeq( Tslab, GHI, nH, y, tauHI_eff )
  dtauHI = dz * kpc_in_cm * nH * xHI * sigma_th

  itauH = 1
  tauH_table = tauH_arr(iTauH)        


  do 

     tauH_prior = tauH
     xHI_prior  = xHI
     GHI_prior  = GHI

     tauH   = tauH + dtauH
     tauHI  = tauHI + dtauHI
     GHI    = gammaHI_from_mini_spec_shield_plane( mini_spec, tauHI_th=tauHI ) 
     tauHI_eff = -log( GHI/GHI_thin )
     xHI    = analytic_xHIeq( Tslab, GHI, nH, y, tauHI_eff )
     dtauHI = dz * kpc_in_cm * nH * xHI * sigma_th
    
     if ( tauH >= tauH_table ) then

        tauH_lo = tauH_prior
        tauH_hi = tauH

        frac = (tauH_table - tauH_lo) / (tauH_hi - tauH_lo)        

        if (GHI > 1.0d-99) then
           log_table_eq(itauH) = &
                log10(GHI_prior) + frac * ( log10(GHI) - log10(GHI_prior) )
        else
           log_table_eq(itauH) = -99.0d0
        endif

        if (itauH == tauH_bins) exit
        itauH = itauH + 1
        tauH_table = tauH_arr(itauH)
        
     endif
     

  end do


 
end function create_and_solve_slab






! loop through all the input parameters calculating
! slabs for them and then interpolates into the table
!======================================================
subroutine create_slab_lookup_table(z, qmax, log_tauH_max, Nlayers, &
                                    N_tauH, N_Kin, N_tauHI, N_T)

  real(r8b), intent(in) :: z    ! redshift we want UV background
  real(r8b), intent(in) :: qmax ! cut-off frequency in Rydbergs
  real(r8b), intent(in) :: log_tauH_max ! max tauH input value
  integer(i8b), intent(in) :: Nlayers  ! number of layers / slab
  integer(i8b), intent(in) :: N_tauH
  integer(i8b), intent(in) :: N_Kin
  integer(i8b), intent(in) :: N_tauHI
  integer(i8b), intent(in) :: N_T

  real(r8b) :: nH_fix, Lslab, Tslab
  real(r8b) :: GHI_goal, GHI_trial, GHI_thin
  real(r8b) :: tauHI_th
  real(r8b) :: tauH, tauH_hi, tauH_lo, frac, xHIeq, logxHIeq
  real(r8b) :: T, CI, RC
  integer(i8b) :: i, j, iK, iHI, iT, iS, tid, nthreads
  character(clen) :: msg

  real(r8b) :: log_qmin, log_qmax

  type(mini_spec_type) :: mini_spec

  integer(i8b) :: Ns, iS_count

  real(r8b), allocatable :: tauH_arr(:)

  real(r8b) :: log_min_vals(0:ndims-1) 
  real(r8b) :: log_max_vals(0:ndims-1) 
                                                                                
  real(r8b) :: min_vals(0:ndims-1) 
  real(r8b) :: max_vals(0:ndims-1) 

  integer(i8b) :: n_bins(0:ndims-1) 
  
  logical :: caseA

  character(clen) :: out_file_name
  character(clen) :: tau_str
  character(clen) :: z_str
  character(clen) :: qmax_str
  character(clen) :: N_str(4)
  character(clen) :: N_name
  integer :: indx

  logical :: Monochromatic
  real(r8b) :: MonoRydbergs
  character(clen) :: MonoNormMethod

  ! set config variables that are used in slab calc
  !----------------------------------------------
  Monochromatic = .false.
  MonoRydbergs = 0.0
  MonoNormMethod = 'not_used'

  ! set default values for table min and max values
  ! [ tauH, K, tauHI, T]
  !----------------------------------------------
  log_min_vals(0:ndims-1) = (/ 0.0d0,  -17.0d0,  0.0d0,  4.0d0 /) 
  log_max_vals(0:ndims-1) = (/ 0.0d0,  -8.0d0,   4.0d0,  5.0d0 /) 
  
  min_vals(0:ndims-1) = 10.0d0**log_min_vals                       
  max_vals(0:ndims-1) = 10.0d0**log_max_vals                       

  n_bins = (/ N_tauH, N_Kin, N_tauHI, N_T /)

  allocate( tauH_arr(n_bins(0) ) )

  log_max_vals(0) = log_tauH_max
  max_vals(0) = 10.0d0**log_tauH_max


  ! create output file name
  !----------------------------------------------
  write(qmax_str,'(I20)') int( qmax )

  write(N_str(1),'(I20)') N_tauH
  write(N_str(2),'(I20)') N_Kin
  write(N_str(3),'(I20)') N_tauHI
  write(N_str(4),'(I20)') N_T

  write(N_name,*) 'N'  // trim(adjustl(N_str(1))) // &
                  'xN' // trim(adjustl(N_str(2))) // &
                  'xN' // trim(adjustl(N_str(3))) // &
                  'xN' // trim(adjustl(N_str(4))) 

  write(*,*) trim(adjustl(N_name))


  write(tau_str,'(ES12.4)') log_tauH_max
  write(z_str, '(A,F6.2)') "z", z 
  indx = scan( z_str, "." )
  z_str(indx:indx) = "p"

  write(out_file_name,'(A,I10.10,A)') &
       '../data/slab_lookup_tables/hm01qg_qmax' // &
       trim(adjustl(qmax_str)) // '/GHI/' // trim(adjustl(N_name)) // &
       '/' // trim(adjustl(z_str)) // '/Nl_', Nlayers, '_hmax_' // &
       adjustl(trim(tau_str)) // '.hdf5'

  ! report some input
  !-----------------------------------------------

  write(*,*) 
  write(*,*) 'log tauH max: ', log_max_vals(0)
  write(*,*) 'Nlayers:      ', Nlayers
  write(*,*) 'z = ', z
  write(*,*) 'qmax = ', qmax
  write(*,*) 'tauH max:     ', max_vals(0)
  write(*,*) 
  write(*,*) 'outfile: ', trim(out_file_name)


  write(*,*) 
  call time_message(stdout,' initializing') 
  write(*,*) 

  if (Nlayers < n_bins(0)) then
     write(*,*) 'Nlayers needs to be greater than tauH_bins!'
     stop
  endif



  ! read spectra file 
  !--------------------------------------------------------------
  call read_ion_table_file( "../data/ionization_tables/h1.hdf5", itab )
  write(*,*) 'read h1 ionization table'

  ! initialize table
  !---------------------------
  do j = 0, ndims-1

     stab%dim(j)%min    = min_vals(j)
     stab%dim(j)%max    = max_vals(j)
     stab%dim(j)%logmin = log_min_vals(j)
     stab%dim(j)%logmax = log_max_vals(j)
     stab%dim(j)%nbins  = n_bins(j)

     allocate( stab%dim(j)%values(0:n_bins(j)-1) )
     do i = 0, n_bins(j)-1 
        stab%dim(j)%values(i) = stab%dim(j)%logmin + &
             i * (stab%dim(j)%logmax - stab%dim(j)%logmin) / &
             (stab%dim(j)%nbins-1)
     end do

  end do


  ! do tauHI individually.  want the first bin to be 0.01, 
  ! second to be 0.1, and the rest log spaced between 1.0 and 4.0
  !--------------------------------------------------------
  stab%dim(2)%values(0) = -2.0d0
  stab%dim(2)%values(1) = -1.0d0
  do i = 2, n_bins(2)-1
     stab%dim(2)%values(i) = stab%dim(2)%logmin + &
          (i-2) * (stab%dim(2)%logmax - stab%dim(2)%logmin) / &
          (stab%dim(2)%nbins-3)
  end do

  write(*,*) 
  write(*,*) 'log tauH bins:  ', n_bins(0)
  write(*,*) 'log K bins:     ', n_bins(1)
  write(*,*) 'log tauHI bins: ', n_bins(2)
  write(*,*) 'log T bins:     ', n_bins(3)
  write(*,*) 

  allocate( stab%table( 0:n_bins(0)-1, &
                        0:n_bins(1)-1, &
                        0:n_bins(2)-1, &
                        0:n_bins(3)-1 ) )

  allocate( xHI_eq_ci( 0:n_bins(3) ) )


  ! calculate collisional equilibrium values
  !--------------------------------------------------------------
  write(*,*) 
  write(*,*) 'Collisional Ionization Equilibrium'
  write(*,*) '    T              xHI            log xHI        '
  write(*,*) '-------------------------------------------------'
  do iT = 0, n_bins(3)-1 
     T = 10.0d0**stab%dim(3)%values(iT)
     RC = Hui_HII_recombA( T )
     CI = Hui_HI_col_ion( T )
     xHI_eq_ci(iT) = RC / (RC+CI)
     write(*,'(3ES15.4)') T, xHI_eq_ci(iT), log10(xHI_eq_ci(iT))
  end do
  write(*,*) 



  
  ! report limits
  !------------------------------------------------------
  nH_fix = 1.0d0
  sigma_th = HI_photo_cs_verner( 1.0d0 )        
  Lslab = max_vals(0) / (nH_fix * sigma_th) / kpc_in_cm * (1.01d0)


  write(*,*) 'calculating slab lookup table'
  write(*,*) ' entries = ', product( n_bins )
  write(*,*) 
  write(*,*) 'min/max log tauH:  ', log_min_vals(0), log_max_vals(0)
  write(*,*) 'min/max log K:     ', log_min_vals(1), log_max_vals(1)
  write(*,*) 'min/max log tauHI: ', log_min_vals(2), log_max_vals(2)
  write(*,*) 'min/max log T:     ', log_min_vals(3), log_max_vals(3)

  write(*,*) 
  write(*,*) 'sig_th      = ', sigma_th 
  write(*,*) 'nH          = ', nH_fix
  write(*,*) 'tauH thru   = ', Lslab * kpc_in_cm * nH_fix * sigma_th
  write(*,*) 'tauH/layer  = ', Lslab * kpc_in_cm * nH_fix * sigma_th / Nlayers 
  write(*,*) 'Lslab [kpc] = ', Lslab 
  write(*,*) 
  call time_message(stdout,' beginning parallel section') 
  write(*,*)   

  log_qmin = log10( 1.0d0 )
  log_qmax = log10( qmax  )

  iS_count = 0

!$OMP parallel shared(nthreads,Ns,itab) private(tid, iS, iK, iHI, iT, GHI_goal, tauHI_th, Tslab, GHI_trial, mini_spec)



  tid = omp_get_thread_num()
  write(*,*) 'Im thread ', tid
  if (tid==0) then 
     nthreads = omp_get_num_threads()
     write(stdout,*) '  nthreads = ', nthreads
     write(stdout,*) ''
     call flush(stdout)
  endif
  

!$OMP BARRIER       


  !==============================
  ! make all the slabs
  !==============================

  ! for every choice of GammaHI/nH, tauHI and Temp. 
  ! solve a new slab for multiple tauH's
  !--------------------------------------------------------------

  Ns = n_bins(1) * n_bins(2) * n_bins(3)
  write(*,*) 'Nslabs = ', Ns


!$OMP do
  slab_loop: do iS = 0, Ns - 1

     write(*,*) 'tid, iS: ', tid, iS

     iS_count = iS_count + 1
     if ( mod(iS_count,Ns/100) == 0 ) then
        write(*,'(A,I10,A,I10)') 'done ', iS_count, ' of ', Ns
     endif
     
     iK  = mod( iS, n_bins(1) ) 
     iHI = mod( iS/n_bins(1), n_bins(2))
     iT  = iS / (n_bins(1)*n_bins(2))

     write(*,*) 'iK,iHI,iT: ', iK,iHI,iT


     ! find the correct GHI for this K
     !----------------------------------
     GHI_goal = 10.0d0**stab%dim(1)%values(iK) * nH_fix
     
     ! set tauHI @ nu_th
     !----------------------------------
     tauHI_th = 10.0d0**stab%dim(2)%values(iHI)
!     if (tauHI_th > 1.0d0) then 
!        caseA = .false.
!     else
        caseA = .true.
!     endif

     write(*,*) 'tauHI_th = ', tauHI_th


     ! set temperature
     !-------------------------
     Tslab = 10.0**stab%dim(3)%values(iT)


     ! attenuate mini spectrum flux from tauHI_th
     !---------------------------------------------------------------
     write(*,*) 'log_qmin: ', log_qmin
     write(*,*) 'log_qmax: ', log_qmax
     write(*,*) 'z: ', z
     write(*,*) 'Monochromatic: ', Monochromatic
     write(*,*) 'MonoRydbergs: ', MonoRydbergs
     write(*,*) 'MonoNormMethod: ', trim(MonoNormMethod)

     call set_mini_spectrum(log_qmin, log_qmax, z, 1.0d0, Monochromatic, &
          MonoRydbergs, MonoNormMethod, itab, mini_spec)
     write(*,*) 'set minispec'


     mini_spec%flux = mini_spec%flux * exp( -tauHI_th * mini_spec%sigma_ratio )



!     write(*,*) minval(mini_spec%flux)
!     mini_spec%logflux = log10( mini_spec%flux )
     
     ! normalize mini spectrum to produce GHI in optically thin limit
     !---------------------------------------------------------------
     GHI_trial = gammaHI_from_mini_spec_shield_plane( mini_spec,tauHI_th=0.0d0) 
     mini_spec%flux = ( GHI_goal / GHI_trial ) * mini_spec%flux
     GHI_trial = gammaHI_from_mini_spec_shield_plane( mini_spec,tauHI_th=0.0d0) 

         
     ! solve a new slab each time we get here.
     ! we want a slab that covers the required 
     ! range in tauH
     !------------------------------------------
     
     tauH_arr(:) = 10.0d0**stab%dim(0)%values(:)
          
     stab%table(:, iK, iHI, iT) = &
          create_and_solve_slab( nH_fix, GHI_goal, Lslab, Tslab, &
          Nlayers, mini_spec, stab%dim(0)%nbins, tauH_arr ) 

     
  end do slab_loop

     
!$OMP end do
!$OMP end parallel


  call time_message(stdout,' finished parallel part') 

  call time_message(stdout,' writing lookup table') 

  call write_slab_lookup_table(out_file_name)


end subroutine create_slab_lookup_table



! writes the table to an HDF5 file
!======================================================
subroutine write_slab_lookup_table( filename )
  character(clen), intent(in) :: filename

  integer(i4b) :: lun
  character(clen) :: data_tag
  character(clen) :: attr_tag
  character(clen) :: VarDescription


  call hdf5_create_file( lun, filename )

  ! write table values
  !---------------------------------------------------------------
  data_tag = table_name
  call hdf5_write_data( lun, data_tag, stab%table )

  VarDescription = 'Log10 of Equilibrium Photo-ionization rate / nH [cm^3/s] dimensions=[tauH, K_in=GHI/nH, tauHI, T]'
  attr_tag = trim(data_tag) // '/VarDescription'
  call hdf5_write_attribute( lun, attr_tag, VarDescription )


  ! write tauH = H optical depth at nu_th
  !---------------------------------------------------------------
  data_tag = 'log_tauH'
  call hdf5_write_data( lun, data_tag, stab%dim(0)%values )

  VarDescription = 'Log10(tauH) = Log10( depth THRU slab * n_H * sigma_HI_th )'
  attr_tag = trim(data_tag) // '/VarDescription'
  call hdf5_write_attribute( lun, attr_tag, VarDescription )


  ! write K = GHI/nH
  !---------------------------------------------------------------
  data_tag = 'log_K'
  call hdf5_write_data( lun, data_tag, stab%dim(1)%values )

  VarDescription = 'Log10(K) = Log10(GHI/nH) = ' // &
       'Log10(incident H photo-ionization rate / Hydrogen number density) [cm^3/s]'
  attr_tag = trim(data_tag) // '/VarDescription'
  call hdf5_write_attribute( lun, attr_tag, VarDescription )


  ! write tauHI = HI optical depth at nu_th
  !---------------------------------------------------------------
  data_tag = 'log_tauHI'
  call hdf5_write_data( lun, data_tag, stab%dim(2)%values )

  VarDescription = 'Log10(tauHI) = Log10( depth TO slab * n_HI * sigma_HI_th )'
  attr_tag = trim(data_tag) // '/VarDescription'
  call hdf5_write_attribute( lun, attr_tag, VarDescription )


  ! write T = Temperature [K]
  !---------------------------------------------------------------
  data_tag = 'log_T'
  call hdf5_write_data( lun, data_tag, stab%dim(3)%values )

  VarDescription = 'Log10(T) = Log10(temperature) [Kelvin]'
  attr_tag = trim(data_tag) // '/VarDescription'
  call hdf5_write_attribute( lun, attr_tag, VarDescription )



  call hdf5_close_file(lun)


end subroutine write_slab_lookup_table



! returns min look up values
!-------------------------------
subroutine return_min_table_vals(min_vals) 
  real(r8b) :: min_vals(0:3)
  integer(i4b) :: i
  do i = 0,3
     min_vals(i) = stab%dim(i)%min
  end do
end subroutine return_min_table_vals

! returns max look up values
!-------------------------------
subroutine return_max_table_vals(max_vals)
  real(r8b) :: max_vals(0:3)
  integer(i4b) :: i
  do i = 0,3
     max_vals(i) = stab%dim(i)%max
  end do
end subroutine return_max_table_vals



! given K = GHI/nH, tauHI, tauH, and T, gets an 
! xHI_eq value from the table (HI=tauHI, H=tauH)
!======================================================
function interpolate_slab_table(H, K, HI, T) result(intrp)
  real(r8b) :: H, K, HI, T, intrp

  real(r8b) :: w_llll, w_hlll, w_lhll, w_hhll
  real(r8b) :: w_llhl, w_hlhl, w_lhhl, w_hhhl
  real(r8b) :: w_lllh, w_hllh, w_lhlh, w_hhlh
  real(r8b) :: w_llhh, w_hlhh, w_lhhh, w_hhhh

  real(r8b) :: vals(0:ndims-1), log_vals(0:ndims-1)
  real(r8b) :: wlo(0:ndims-1), whi(0:ndims-1), delta(0:ndims-1)
  integer(i8b) :: ilo(0:ndims-1), ihi(0:ndims-1)
  real(r8b) :: log_intrp

  integer(i8b) :: i, j


  ! special case for temperature
  !--------------------------------
!  if ( T > stab%dim(3)%max ) T = stab%dim(3)%max

  ! special case for K
  !--------------------------------
!  if ( K < stab%dim(1)%min ) K = stab%dim(1)%min


  vals = (/ H, K, HI, T /)

  do j = 0, ndims-1


     ! make sure we're in bounds
     !--------------------------
     if ( vals(j) < stab%dim(j)%min ) then
        vals(j) = stab%dim(j)%min
!        write(*,*) ' interpolate slab error: '
!        write(*,*) 'j, val, min: ', j, vals(j), stab%dim(j)%min
!        stop 
     endif

     if ( vals(j) > stab%dim(j)%max ) then
        vals(j) = stab%dim(j)%max
!        write(*,*) ' interpolate slab error: '
!        write(*,*) 'j, val, max: ', j, vals(j), stab%dim(j)%max
!        stop 
     endif

     log_vals(j) = log10( vals(j) )


     ! find bracketing indices
     !--------------------------
     over_bins: do i = 0, stab%dim(j)%nbins-1
        if ( stab%dim(j)%values(i) > log_vals(j) ) then
           ihi(j) = i
           ilo(j) = ihi(j)-1
           exit over_bins
        end if
        ihi(j) = stab%dim(j)%nbins-1
        ilo(j) = ihi(j)-1
     end do over_bins


     ! calculate weights
     !--------------------------
     delta(j) = stab%dim(j)%values( ihi(j) ) - stab%dim(j)%values( ilo(j) )
     wlo(j)    = ( log_vals(j) - stab%dim(j)%values( ilo(j) ) ) / delta(j)
     whi(j)    = ( stab%dim(j)%values( ihi(j) ) - log_vals(j) ) / delta(j)
     
     
  end do


  w_llll = whi(0) * whi(1) * whi(2) * whi(3)
  w_hlll = wlo(0) * whi(1) * whi(2) * whi(3)
  w_lhll = whi(0) * wlo(1) * whi(2) * whi(3)
  w_hhll = wlo(0) * wlo(1) * whi(2) * whi(3)
  w_llhl = whi(0) * whi(1) * wlo(2) * whi(3)
  w_hlhl = wlo(0) * whi(1) * wlo(2) * whi(3)
  w_lhhl = whi(0) * wlo(1) * wlo(2) * whi(3)
  w_hhhl = wlo(0) * wlo(1) * wlo(2) * whi(3)
  w_lllh = whi(0) * whi(1) * whi(2) * wlo(3)
  w_hllh = wlo(0) * whi(1) * whi(2) * wlo(3)
  w_lhlh = whi(0) * wlo(1) * whi(2) * wlo(3)
  w_hhlh = wlo(0) * wlo(1) * whi(2) * wlo(3)
  w_llhh = whi(0) * whi(1) * wlo(2) * wlo(3)
  w_hlhh = wlo(0) * whi(1) * wlo(2) * wlo(3)
  w_lhhh = whi(0) * wlo(1) * wlo(2) * wlo(3)
  w_hhhh = wlo(0) * wlo(1) * wlo(2) * wlo(3)



  log_intrp = &
       w_llll * stab%table( ilo(0), ilo(1), ilo(2), ilo(3) ) + &
       w_hlll * stab%table( ihi(0), ilo(1), ilo(2), ilo(3) ) + &
       w_lhll * stab%table( ilo(0), ihi(1), ilo(2), ilo(3) ) + &
       w_hhll * stab%table( ihi(0), ihi(1), ilo(2), ilo(3) ) + &
       w_llhl * stab%table( ilo(0), ilo(1), ihi(2), ilo(3) ) + &
       w_hlhl * stab%table( ihi(0), ilo(1), ihi(2), ilo(3) ) + &
       w_lhhl * stab%table( ilo(0), ihi(1), ihi(2), ilo(3) ) + &
       w_hhhl * stab%table( ihi(0), ihi(1), ihi(2), ilo(3) ) + &
       w_lllh * stab%table( ilo(0), ilo(1), ilo(2), ihi(3) ) + &
       w_hllh * stab%table( ihi(0), ilo(1), ilo(2), ihi(3) ) + &
       w_lhlh * stab%table( ilo(0), ihi(1), ilo(2), ihi(3) ) + &
       w_hhlh * stab%table( ihi(0), ihi(1), ilo(2), ihi(3) ) + &
       w_llhh * stab%table( ilo(0), ilo(1), ihi(2), ihi(3) ) + &
       w_hlhh * stab%table( ihi(0), ilo(1), ihi(2), ihi(3) ) + &
       w_lhhh * stab%table( ilo(0), ihi(1), ihi(2), ihi(3) ) + &
       w_hhhh * stab%table( ihi(0), ihi(1), ihi(2), ihi(3) ) 

  
  intrp = 10.0d0**log_intrp




end function interpolate_slab_table




end module slab_lookup_table_class






