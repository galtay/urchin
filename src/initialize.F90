!> \file initialize.F90 

!> \brief Initialization module
!! 
!! This module contains the subroutines to initialize all of the non  
!! particle system variables.     
!<
module initialize_mod
use myf90_mod
use config_mod, only: read_config_file
use mt19937_mod, only: init_mersenne_twister
use atomic_rates_mod, only: read_atomic_rates_file
use atomic_rates_mod, only: write_atomic_rates_to_log_file
use ion_table_class, only: read_ion_table_file
use slab_lookup_table_class, only: load_slab_lookup_table
use main_input_mod, only: get_planning_data

use global_mod, only: write_units_to_log_file
use global_mod, only: write_particle_headers_to_log_file

use atomic_rates_mod, only: get_atomic_rates 
use config_mod, only: CV
use global_mod, only: GV
use global_mod, only: gconst
use global_mod, only: rtable
use global_mod, only: xHII_k, cmbT_k, isoT_k
use global_mod, only: h1_itab, he1_itab, he2_itab
use particle_system_mod, only: calc_bytes_per_particle
implicit none
private

public :: initialize

  real, parameter :: zero = 0.0d0

  contains

!> Main initialization routine
!---------------------------------
subroutine initialize(config_file)

  character(clen), intent(in) :: config_file !< configuration file
  
  ! these subroutines are in other files 
  ! indicated in the use statements above

  call read_config_file(config_file)

  call init_mersenne_twister(CV%IntSeed)

  call read_atomic_rates_file(rtable, CV%AtomicRatesFile)

  call write_atomic_rates_to_log_file(rtable, CV%OutputDir)

  call mywrite('reading CLOUDY ionization tables:',verb=1)
  call mywrite('',verb=2)
  call read_ion_table_file( trim(CV%ItabDir)//'/h1.hdf5',  h1_itab )
  call read_ion_table_file( trim(CV%ItabDir)//'/he1.hdf5', he1_itab )
  call read_ion_table_file( trim(CV%ItabDir)//'/he2.hdf5', he2_itab )

#ifndef oldRT
  call load_slab_lookup_table( CV%SlabLookupFile )
#endif

  call get_planning_data()

  call write_units_to_log_file()

  call write_particle_headers_to_log_file()

  call get_atomic_rates(1.0d4, rtable, xHII_k)
  call get_atomic_rates( GV%Tcmb_cur, rtable, cmbT_k )
  if (CV%IsoTemp > 0.0) then
     call get_atomic_rates(CV%IsoTemp, rtable, isoT_k)
  end if  

  call calc_bytes_per_particle(GV%bytesperpar)

  GV%rayn = 0
  GV%MB = zero
  GV%OutputIndx = 0

  ! open log files
  !-----------------
  GV%pardata_file = trim(CV%OutputDir) // "/particle_data.log"
  call open_formatted_file_w(GV%pardata_file,GV%pardatalun)


          
end subroutine initialize







  





end module initialize_mod
