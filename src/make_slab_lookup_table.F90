! short program to write a table
!------------------------------------

program make_slab_lookup_table
use myf90_mod
use slab_lookup_table_class
implicit none


type(command_line_type) :: cmnd
real(r8b) :: log_tauH_max
integer(i8b) :: Nlayers
real(r8b) :: z
real(r8b) :: qmax
integer(i8b) :: N_tauH
integer(i8b) :: N_Kin
integer(i8b) :: N_tauHI
integer(i8b) :: N_T

! read command line and reset tauH MAX val.
!----------------------------------------------
cmnd = myf90_initialize_command_line(verb=3)

if (cmnd%nargs /= 8) then 
   write(*,*) 'usage: '
   write(*,*) 'make_slab_lookup_table log_tauH_max Nlayers z qmax'
   write(*,*) '                       N_tauH, N_Kin, N_tauHI, N_T'
   stop
endif



read(cmnd%args(1),*) log_tauH_max

read(cmnd%args(2),*) Nlayers

read(cmnd%args(3),*) z

read(cmnd%args(4),*) qmax

read(cmnd%args(5),*) N_tauH

read(cmnd%args(6),*) N_Kin

read(cmnd%args(7),*) N_tauHI

read(cmnd%args(8),*) N_T



call create_slab_lookup_table(z, qmax, log_tauH_max, Nlayers, &
                              N_tauH, N_Kin, N_tauHI, N_T)



end program make_slab_lookup_table
