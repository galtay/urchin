!> \file output_gadget_public.F90

!> \brief The module that handles output in 
!! the style of the public version of Gadget
!<
module output_gadget_public_mod
use myf90_mod
! types
use particle_system_mod, only: particle_system_type
use gadget_sphray_header_class, only: gadget_sphray_header_type
use gadget_general_class, only: gadget_constants_type
! routines
use particle_returns_mod, only: return_H_mf
use global_mod, only: recall_saved_ghead
use gadget_sphray_header_class, only: gadget_sphray_header_write_lun
! variables
use config_mod, only: CV
use global_mod, only: GV
implicit none
private


public :: output_snap_gadget_public


contains


subroutine output_snap_gadget_public(psys, NumFiles)

  character(clen), parameter :: myname="output_snap_gadget_public"
  logical, parameter :: crash = .true.
  integer, parameter :: verb = 2

  type(particle_system_type), intent(in) :: psys
  integer(i4b) :: NumFiles  !< number of files to put output in

  type(gadget_constants_type) :: gconst

  integer(i4b) :: Nwrote   
  integer(i4b) :: ifile
  integer(i4b) :: ngas1
  integer(i4b) :: ngood1
  integer(i4b) :: igood
  integer(i4b) :: indx
  integer(i4b) :: lun
  integer(i4b) :: i
  real(r4b) :: Hmf
  real(r4b) :: mu

  character(3) :: label
  character(4) :: ext
  character(clen) :: filename
  
  type(gadget_sphray_header_type) :: ghead

  real(r4b), allocatable :: rblock3(:,:)
  real(r4b), allocatable :: rblock(:)
  integer(i4b), allocatable :: iblock(:)

  
  write(label,'(I3.3)') GV%OutputIndx

  Nwrote = 0

  do ifile = 0, NumFiles-1

     call recall_saved_ghead( ifile, ghead )


     ! adjust if we have a selection
     !-------------------------------
     if (CV%SelectRadius > 0) then
        ghead%nfiles = 1
        ghead%npar_file    = 0
        ghead%npar_file(0) = size(psys%par)
        ghead%npar_all     = 0
        ghead%npar_all(0)  = size(psys%par)
     endif

     ngood1 = count( psys%par(:)%mask )
     ngas1 = ghead%npar_file(0)


     ! form output file name
     !--------------------------------
     write(ext,'(I3)') ifile

     filename = trim(CV%OutputDir) // "/" // trim(CV%OutputFileBase) // "_" // label 
     if (NumFiles > 1) then
        filename = trim(filename) // "." // trim(adjustl(ext))
     end if

     call mywrite("writing snapshot state to "//trim(filename), verb)
     call mywrite('',verb) 


    

     ! open and write to file
     !--------------------------------
     call open_unformatted_file_w(filename,lun)
     


     !================================================================
     ! write Header
     !================================================================
     call gadget_sphray_header_write_lun(ghead,lun)



     !================================================================
     ! write data
     !
     ! out_list makes sure the particles are output in the 
     ! order they were read in.  mask is default true, but 
     ! can be set to false during the calculation
     ! for default output we only write URCHIN specific data
     ! for selections, we write all the data
     !================================================================
     if (CV%DoSelection) then

        allocate( rblock3(3,ngood1) )
        
        ! write positions (selection only)
        !-----------------------------------
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
        write(lun) rblock3


        ! write velocities (selection only)
        !-----------------------------------
#ifdef incVel
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(indx)%mask) then
              igood = igood + 1
              rblock3(1,igood) = psys%par( indx )%vel(1)
              rblock3(2,igood) = psys%par( indx )%vel(2)
              rblock3(3,igood) = psys%par( indx )%vel(3)
           endif
        end do
#else
        rblock3 = 0.0e0
#endif
        write(lun) rblock3
        deallocate( rblock3 )

     endif


     ! write IDs
     !---------------------------------
     allocate( iblock(ngood1) )
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           iblock(igood) = psys%par( indx )%id
        endif
     end do
     write(lun) iblock
     deallocate( iblock )
     

     allocate( rblock(ngood1) )


     if (CV%DoSelection) then

        ! write mass (selection only)
        !-------------------------------------
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%mass
           endif
        end do
        write(lun) rblock

     endif
        
     ! write internal energy 
     !------------------------------------------------
     igood = 0     
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           Hmf = return_H_mf( psys%par( indx ) )
           mu  = 4.0d0 / ( 3.0d0 * Hmf + 1.0d0 + 4.0d0 * Hmf * psys%par(indx)%ye )
           rblock(i) = ( gconst%BOLTZMANN * psys%par(indx)%T ) / &
                ( (gconst%GAMMA - 1.0d0) * mu * gconst%PROTONMASS  )
           rblock(i) = rblock(i) * GV%cgs_mass / GV%cgs_enrg
        endif
     end do
     write(lun) rblock


     if (CV%DoSelection) then

        ! write density (selection only)
        !---------------------------------
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%rho
           endif
        end do
        write(lun) rblock
        
     endif


     ! write electon fraction ye = ne/nH
     !--------------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%ye
        endif
     end do    
     write(lun) rblock
     

     ! write neutral H fraction xHI = nHI/nH
     !----------------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%xHI
        endif
     end do
     write(lun) rblock
     

     if (CV%DoSelection) then

        ! write smoothing lengths (selection only)
        !--------------------------------------------
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%hsml
           endif
        end do
        write(lun) rblock

     endif


     ! write temperatures
     !---------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%T
        endif
     end do    
     write(lun) rblock
     


     if (CV%DoSelection) then

#ifdef incHmf
        ! write Hydrogen mass fractions (selection only)
        !---------------------------------
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%Hmf
           endif
        end do
        write(lun) rblock
#endif
        

#ifdef incHemf
        ! write Helium mass fractions (selection only)
        !---------------------------------
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%Hemf
           endif
        end do
        write(lun) rblock
#endif
        
     endif



#ifdef incHe
     ! write Helium neutral fractions
     !---------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%xHeI
        endif
     end do
     write(lun) rblock

     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%xHeII
        endif
     end do
     write(lun) rblock
#endif
     
     
#ifdef outGammaHI
     ! write photo-ionization rates
     !---------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%gammaHI
        endif
     end do
     write(lun) rblock
#endif

 
#ifdef incCloudy
     ! write optically thin xHI
     !---------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%xHI_cloudy
        endif
     end do
     write(lun) rblock
#endif


     if (CV%DoSelection) then
#ifdef incEOS
        ! write EoS variable
        !---------------------------------
        igood = 0
        do i = 1, ngas1
           indx = psys%out_list(Nwrote+i)
           if (psys%par(i)%mask) then
              igood = igood + 1
              rblock(igood) = psys%par( indx )%eos
           endif
        end do
        write(lun) rblock
#endif

                
#ifdef incSFR
     ! write SFR 
     !---------------------------------
     igood = 0
     do i = 1, ngas1
        indx = psys%out_list(Nwrote+i)
        if (psys%par(i)%mask) then
           igood = igood + 1
           rblock(igood) = psys%par( indx )%sfr
        endif
     end do
     write(lun) rblock
#endif

     endif

     close(lun)
     Nwrote = Nwrote + ngas1
     deallocate( rblock )
     

  end do


end subroutine output_snap_gadget_public


end module output_gadget_public_mod
