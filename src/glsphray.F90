!> \file glsphray.F90

!> \brief SPHRAY basic operation program + OpenGL viewer
!! 
!! this program does the same thing as sphray.f90 but contains an 
!! extra subroutine to launch a thread that calls the OpenGL viewer.  
!<
program glsphray
use myf90_mod
use initialize_mod, only: initialize
use mainloop_mod, only: mainloop
use viewermod
implicit none

   character(clen) :: config_file    !< configuration file
   type(command_line_type) :: cmnd   !< command line variable

   integer, parameter :: verb = 3    !< verbosity before config file is read
   logical, parameter :: crash = .true.               !< stop on error
   character(clen), parameter :: myname = "glsphray"  !< routine name

   ! set initial myf90_verbosity (will be changed once config file is read)
   ! 3 lets any statement with verbosity 3 or less through
   !========================================================================
   myf90_verbosity = verb

   ! test user defined variables sizes
   !========================================================================
   if (myf90_test_var_sizes() /= 0) then
      call mywrite('  ** warning selected____kind sizes are off', verb )
      call system( 'sleep(2)' )     
   endif

   ! initialize command line variable
   !========================================================================
   cmnd = myf90_initialize_command_line(verb)

   ! check command line arguments
   !========================================================================
   if (cmnd%nargs /= 1 ) then 
      call myerr( 'usage: ./glsphray config_file', myname, crash )
   end if

   ! clear terminal and print splash screen
   !========================================================================
   call system("clear")


   write(*,*) "**************************************************************"
   write(*,*) "**************************************************************"
   write(*,*) "         Smoothed Particle Hydrodynamics Ray Tracer           "
   write(*,*) "                        SPHRAY 0.9                            "
   write(*,*) "                  (with OpenGL Viewer)                        "
   write(*,*) 
   write(*,*) "  developed by,                                               "
   write(*,*)
   write(*,*) "  Gabriel Altay (gabriel.altay@gmail.com) and                 "
   write(*,*) "  Inti Pelupessy                                              "
   write(*,*) "**************************************************************"
   write(*,*) "**************************************************************"
   write(*,*) 
   write(*,*) 
   write(*,*) 

   ! initialize and call viewer and mainloop()
   !========================================================================
   config_file = cmnd%args(1)
   call initialize(config_file)
      
   call viewer()

   call mainloop()

   write(*,*) " *** GAME OVER ***"
   write(*,*) "you can close the viewer safely now.."

   call viewer()


end program glsphray
