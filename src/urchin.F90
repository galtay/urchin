!> \file urchin.F90

!> \brief URCHIN - calculates equilibrium background shielding
!! 
!! This program calls initialize and main_loop 
!<
program urchin
use myf90_mod
use initialize_mod, only: initialize
use main_loop_mod, only: main_loop
implicit none


  
  character(clen) :: config_file    !< configuration file
  type(command_line_type) :: cmnd   !< command line variable
  
  integer, parameter :: verb = 0    !< verbosity before config file is read
  logical, parameter :: crash = .true.             !< stop on error
  character(clen), parameter :: myname = "urchin"  !< routine name
  


  ! set initial myf90_verbosity (will be changed once config file is read)
  ! 3 lets any statement with verbosity 3 or less through
  !========================================================================
  myf90_verbosity = verb
    


  ! initialize command line variable
  !========================================================================
  cmnd = myf90_initialize_command_line(verb=3)
  


  ! check command line arguments
  !========================================================================
  if (cmnd%nargs /= 1 ) then 
     call myerr( 'usage: ./urchin config_file', myname, crash )
  end if
  


  ! clear terminal and print splash screen
  !========================================================================
  write(*,*) 
  write(*,*) "**************************************************************"
  write(*,*) "**************************************************************"
  write(*,*) "                        URCHIN 1.0                            "
  write(*,*)
  write(*,*) "  developed by,                                               "
  write(*,*)
  write(*,*) "  Gabriel Altay (gabriel.altay@gmail.com)                     "
  write(*,*) "  Tom Theuns    (tom.theuns@googlemail.com)                   "
  write(*,*) "**************************************************************"
  write(*,*) "**************************************************************"
  write(*,*) 
  write(*,*) 



  ! initialize and call mainloop()
  !========================================================================
  config_file = cmnd%args(1)

  call initialize(config_file)
  
  call main_loop()
  
  
  
    
end program urchin
