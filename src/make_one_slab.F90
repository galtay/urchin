! short program to write a single slab
!------------------------------------

program make_one_slab
use timers
use myf90_mod
use slab_lookup_table_class

use ion_table_class, only: ion_table_type
use ion_table_class, only: read_ion_table_file
use ion_table_class, only: mini_spec_type
use ion_table_class, only: set_mini_spectrum
use ion_table_class, only: gammaHI_from_mini_spec_thin
use ion_table_class, only: HI_photo_cs_verner
implicit none


type(ion_table_type) :: itab
type(mini_spec_type) :: mini_spec
type(slab_type) :: slab

real(r8b) :: nH
real(r8b) :: GHI_0
real(r8b) :: Lslab
real(r8b) :: Tslab
integer(i8b) :: Nlayers

real(r8b) :: qmax = 4.0d0
real(r8b) :: z = 3.0d0

real(r8b) :: GHI_thin

integer, parameter :: N = 512
real(r8b) :: tauH_min
real(r8b) :: tauH_max
real(r8b) :: dtauH

integer :: i
real(r8b) :: tauH(0:N-1)
real(r8b) :: xHI(0:N-1)

real(r8b) :: sigma_th

real(r8b), parameter :: kpc_in_cm = 3.0856780d21

sigma_th = HI_photo_cs_verner( 1.0d0 )        

! read spectra file 
!--------------------------------------------------------------
call read_ion_table_file( "../data/ionization_tables/h1.hdf5", itab )

! make mini spectrum
!--------------------------------------------------------------
call set_mini_spectrum(min_logryd=0.0d0, max_logryd=log10(qmax), z=z, itab=itab, mini_spec=mini_spec)
GHI_thin = gammaHI_from_mini_spec_thin( mini_spec ) 

! solve for slab
!--------------------------------------------------------------
nH = 1.0d0
GHI_0 = 1.0d-12
Nlayers = 1000000

tauH_min = 1.0d0
tauH_max = 1.0d2
dtauH = (tauH_max - tauH_min) / (N-1)

do i = 0, N-1
   tauH(i) = tauH_min + i * dtauH
end do

Lslab = tauH_max / (sigma_th * nH) / kpc_in_cm * 1.01d0

write(*,*) 
write(*,*) 'Nlayers: ', Nlayers
write(*,*) 'Lslab:   ', Lslab

call time_message(stdout,' start') 
xHI = create_and_solve_slab( nH, GHI_0, Lslab, Tslab, Nlayers, mini_spec, tauH ) 
call time_message(stdout,' finish') 
write(*,*) 

stop 'finished one slab'


end program make_one_slab
