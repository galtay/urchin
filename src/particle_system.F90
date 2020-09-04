!> \file particle_system.F90

!> \brief Particles and box, types and subroutines.
!! 
!<

!> particle, source, and box descriptions
module particle_system_mod
use myf90_mod
use atomic_rates_mod, only: calc_colion_eq_fits
use mt19937_mod, only: genrand_real1
use m_mrgrnk, only: mrgrnk
use config_mod, only: CV
implicit none
private

public :: particle_type
public :: box_type
public :: particle_system_type
public :: transformation_type
public :: make_output_list
public :: orderparticles
public :: orderpsys
public :: calc_bytes_per_particle
public :: scale_comoving_to_physical
public :: scale_physical_to_comoving
public :: set_ye
public :: set_collisional_ionization_equilibrium
public :: enforce_x_and_T_minmax
public :: particle_info_to_screen
public :: adjustbox
public :: create_particle_random_access_list
public :: create_particle_density_access_list
public :: create_particle_radial_access_list
public :: create_particle_slab_access_list
public :: constrain_psys_to_selection


!> Particle type. 
!=========================
type particle_type
   real(r4b)    :: pos(3)     !< x,y,z coordinates
#ifdef incVel
   real(r4b)    :: vel(3)     !< x,y,z velocities
#endif
#ifdef LONGIDs
   integer(i8b) :: id         !< particle id
   integer(i8b) :: orig_order !< readin order
#else
   integer(i4b) :: id         !< particle id
   integer(i4b) :: orig_order !< readin order
#endif
   real(r4b)    :: mass       !< particle mass
   real(r4b)    :: T          !< temperature in K       
   real(r4b)    :: rho        !< density 
   real(r4b)    :: ye         !< electron fraction
   real(r4b)    :: xHI        !< nHI  / nH 
   real(r4b)    :: xHII       !< nHII / nH
   real(r4b)    :: hsml       !< smoothing length

   real(r4b)    :: xHI_tmp    !< to check convergence

#ifdef incH2
   real(r4b)    :: fH2        !< H2 mass fraction = mass of H2 / mass of H
#endif

#ifdef incCloudy
   real(r4b)    :: xHI_cloudy !< cloudy eq solutions
#endif

#ifdef Tinit
   real(r4b)    :: Ti         !< initial temperature
#endif

#ifdef incHmf   
   real(r4b)    :: Hmf        !< Hydrogen mass fraction
#endif

#ifdef incHe
   real(r4b)    :: xHeI       !< nHeI   / nHe 
   real(r4b)    :: xHeII      !< nHeII  / nHe
   real(r4b)    :: xHeIII     !< nHeIII / nHe
#endif

#ifdef incHemf
   real(r4b)    :: Hemf       !< Helium mass fraction
#endif

   real(r4b)    :: gammaHI    !< HI photoionization rate

#ifdef incEOS
   real(r4b)    :: eos        !< equation of state variable
#endif


   logical      :: mask       !< just flag for test output

   real(r4b)    :: nH         !< physical atomic H number density [cm^-3]

#ifdef EAGLE
   integer(i4b) :: file_index
   integer(i4b) :: file_offset
#endif  

end type particle_type



!> simulation box and boundary conditions
!==========================================
type box_type
   real(r8b)    :: top(1:3)     !< upper x,y,z coordinates
   real(r8b)    :: bot(1:3)     !< lower x,y,z coordinates
   integer(i8b) :: bbound(1:3)  !< BCs for upper faces (0:vac 1:per -1: ref) 
   integer(i8b) :: tbound(1:3)  !< BCs for lower faces (0:vac 1:per -1: ref) 
end type box_type



!> particles and box
!----------------------------
  type particle_system_type
     type(box_type) :: box                          !< the simulation box     
     type(particle_type), allocatable :: par(:)     !< all particles
#ifdef LONGIDs
     integer(i8b), allocatable :: acc_list(:)       !< access list
     integer(i8b), allocatable :: out_list(:)       !< output indexx
#else
     integer(i4b), allocatable :: acc_list(:)       !< access list
     integer(i4b), allocatable :: out_list(:)       !< output indexx
#endif
     integer(i8b)              :: ngas_max          !< maximum size of gas particle array
  end type particle_system_type                     !< (ngas_max only used for selections)


!> transformation type
!=======================
type transformation_type
   integer(i8b) :: fac(1:3)  !< newpos = fac * (oldpos - shift)
   real(r8b) :: shift(1:3)   !< newpos = fac * (oldpos - shift)
end type transformation_type


real(r8b), parameter, private :: zero = 0.0d0


contains


!> ranks the orig order array to make the output list
!---------------------------------------------------------------------------
 subroutine make_output_list(psys)
   type(particle_system_type) :: psys
   integer(i8b) :: N 

   N = size(psys%par(:))

   if ( allocated( psys%out_list) ) then
      deallocate( psys%out_list )
   endif

   allocate( psys%out_list(N) )
      
   call mrgrnk(psys%par(:)%orig_order, psys%out_list(:))

 end subroutine make_output_list


!> figures out how many bytes  are needed per particle
!=================================================================
subroutine calc_bytes_per_particle(bpp)
  character(clen), parameter :: myname="calc_bytes_per_particle"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2

  integer(i4b), intent(out) :: bpp !< bytes per particle
  character(clen) :: str

  bpp = 12       ! positions

#ifdef incVel
  bpp = bpp + 12 ! velocities  
#endif

#ifdef LONGIDs
  bpp = bpp + 8  ! ID
  bpp = bpp + 8  ! orig order
#else
  bpp = bpp + 4  ! ID
  bpp = bpp + 4  ! orig order
#endif

  bpp = bpp + 4  ! mass
  bpp = bpp + 4  ! temperature
  bpp = bpp + 4  ! rho
  bpp = bpp + 4  ! ye
  bpp = bpp + 8  ! H ionization fractions
  bpp = bpp + 4  ! hsml

  bpp = bpp + 4  ! xHI tmp

#ifdef incH2
  bpp = bpp + 4  ! H2 mass fraction = mass of H2 / mass of H
#endif

#ifdef incCloudy
  bpp = bpp + 4  ! cloudy table xHI
#endif

#ifdef Tinit
  bpp = bpp + 4  ! initial temperature
#endif

#ifdef incHmf
  bpp = bpp + 4  ! H mass fraction
#endif

#ifdef incHe
  bpp = bpp + 12 ! He ionization fractions
#endif

#ifdef incHemf
  bpp = bpp + 4  ! He mass fraction
#endif

  bpp = bpp + 4  ! GammaHI tracking

#ifdef incEOS
  bpp = bpp + 4  ! Equation of State variable
#endif

  write(str,'(A,I4)') "   bytes per particle = ", bpp
  call mywrite(str, verb)
  call mywrite('', verb)
  
end subroutine calc_bytes_per_particle




!> resize psys with selected region only
!------------------------------------------------------
subroutine constrain_psys_to_selection( psys , mode)
  character(clen), parameter :: myname="constrain_psys_to_selection" 
  logical, parameter :: crash=.true.

  type(particle_system_type) :: psys
  integer(i8b), optional     :: mode
  type(particle_system_type) :: psys_temp

  integer(i8b) :: ngas
  logical, allocatable :: keep(:)
  real(r8b) :: center(3), boxlens(3), dr(3), rad2 
  integer(i8b) :: ic, ngaskeep, i, err

  ngas = size(psys%par)
  allocate( keep(ngas) )
  keep = .false. 

  rad2   = CV%SelectRadius**2
  center = (/ CV%SelectXcoord, CV%SelectYCoord, CV%SelectZcoord /)          
  boxlens = psys%box%top - psys%box%bot

  write(*,*) 
  write(*,*) '  Only keeping particles in selection '
  write(*,'(A,3ES12.5)') '     Center: ', center
  write(*,'(A,ES12.5)') '     Radius: ', CV%SelectRadius 
  write(*,'(A,3I4)') '     Top BCs: ', psys%box%tbound
  write(*,'(A,3I4)') '     Bot BCs: ', psys%box%bbound
  write(*,*) 



  ! decide which particles are naughty and which are nice
  !-----------------------------------------------------------!
  do i = 1, ngas

     ! calc seperations in each dimension w/ possible periodic BCs
     !----------------------------------------------------------------
     do ic=1,3
        dr(ic) = abs( psys%par(i)%pos(ic) - center(ic) )
        if ( psys%box%bbound(ic) == 1 .and. psys%box%tbound(ic) == 1 ) then
           if ( dr(ic) > 0.5d0 * boxlens(ic) ) then
              dr(ic) = dr(ic) - 0.5 * boxlens(ic)
           endif
        endif
     enddo
     
     ! decide to keep or not
     !-----------------------------------------
     if ( sum(dr*dr) < rad2 ) then
        keep(i) = .true. 
     endif
        
  end do

  
  ! move good particles to front of par array.
  !------------------------------------
  ngaskeep = 0
  do i = 1, ngas
     if (keep(i)) then        
        ngaskeep = ngaskeep + 1
        psys%par(ngaskeep) = psys%par(i)
     endif
  end do
  if (present(mode)) then
     mode = ngaskeep
     write(*,*) 
     write(*,*) '  ... constrain: ngas particles kept = ', ngaskeep
     write(*,*)
     
     deallocate( keep )
     return
  endif
  
     
  ! re-allocate par array if we need to
  !-----------------------------------------------------------!
  if (ngaskeep /= ngas) then

     allocate(psys_temp%par(ngaskeep), stat=err)
     if (err /= 0) call myerr("failed to allocate psys_temp%par",myname,crash=.true.)

     do i = 1, ngaskeep
        psys_temp%par(i) = psys%par(i)
     end do

     deallocate(psys%par)
     allocate(psys%par(ngaskeep), stat=err)

     do i = 1, ngaskeep
        psys%par(i) = psys_temp%par(i)
     end do

     deallocate(psys_temp%par)

  endif


  write(*,*) 
  write(*,*) '  ngas particles kept = ', ngaskeep
  write(*,*)
  
  deallocate( keep )

  

end subroutine constrain_psys_to_selection





!> allows for accessing the particles in a random order
!------------------------------------------------------
subroutine create_particle_random_access_list( psys )
  type(particle_system_type) :: psys

  integer(i8b) :: i
  integer(i8b) :: n
  real(r4b), allocatable :: randoms(:)

  n = size( psys%par )
  if (.not. allocated(psys%acc_list) ) allocate( psys%acc_list(n) )
  allocate( randoms(n) )

  do i = 1, size(randoms)
     randoms(i) = genrand_real1()
  end do
  
  call mrgrnk( randoms, psys%acc_list )

  deallocate( randoms )


end subroutine create_particle_random_access_list



!> allows for accessing the particles from least to most dense
!--------------------------------------------------------------
subroutine create_particle_density_access_list( psys )
  type(particle_system_type) :: psys

  integer(i8b) :: i
  integer(i8b) :: n
  real(r4b), allocatable :: rhos(:)

  n = size( psys%par )
  if (.not. allocated(psys%acc_list) ) allocate( psys%acc_list(n) )
  allocate( rhos(n) )

  do i = 1,n     
     rhos(i) = psys%par(i)%rho
  end do
  call mrgrnk( rhos, psys%acc_list )
  
  deallocate( rhos )

end subroutine create_particle_density_access_list



!> allows for accessing particles from furthers to nearest the center
!--------------------------------------------------------------
subroutine create_particle_radial_access_list( psys )
  type(particle_system_type) :: psys

  integer(i8b) :: i
  integer(i8b) :: n
  real(r8b), allocatable :: r(:)
  real(r8b) :: center(3)
  real(r8b) :: radius

  center = (/500.0d0, 500.0d0, 500.0d0/)
  radius = 500.0d0

  n = size( psys%par )
  if (.not. allocated(psys%acc_list) ) allocate( psys%acc_list(n) )
  allocate( r(n) )

  do i = 1,n     
     r(i) = radius - sqrt( (center(1)-psys%par(i)%pos(1))**2 + &
                           (center(2)-psys%par(i)%pos(2))**2 + &
                           (center(3)-psys%par(i)%pos(3))**2   )
  end do



  call mrgrnk( r, psys%acc_list )
  
  deallocate( r )

end subroutine create_particle_radial_access_list


!> allows for accessing particles from lowest to highest
!! z-coordinate
!--------------------------------------------------------------
subroutine create_particle_slab_access_list( psys )
  type(particle_system_type) :: psys

  integer(i8b) :: i
  integer(i8b) :: n
  real(r8b), allocatable :: r(:)
!  real(r4b) :: zmax

  n = size( psys%par )
  if (.not. allocated(psys%acc_list) ) allocate( psys%acc_list(n) )
  allocate( r(n) )

  do i = 1,n     
     r(i) = psys%par(i)%pos(3) 
  end do

  call mrgrnk( r, psys%acc_list )
  
  deallocate( r )

end subroutine create_particle_slab_access_list



! this routine rearranges the particles in pars so that they are stored
! in the sequence given by the array order.  ie order = [3,1,2] takes the 
! third particle to the first position, the first particle to the second
! position, and the second particle to the third position.  the array 
! order is not preserved during the routine
!> reorders the particles according to the array order
!===========================================================================
subroutine orderparticles(pars,order)
  type(particle_type), intent(inout) :: pars(:) !< input particles
#ifdef LONGIDs
  integer(i8b), intent(inout) :: order(:)       !< desired order
#else
  integer(i4b), intent(inout) :: order(:)       !< desired order
#endif

  type(particle_type) :: par
  integer(i8b) :: i
  integer(i8b) :: goal
  integer(i8b) :: npar
  
  if (size(pars) /= size(order)) stop "size(pars) /= size(order)"
  npar = size(pars)
  
  do i=1,npar 
     par=pars(i)
     goal=order(i)
     do while(goal < i)
        goal=order(goal)
        order(i)=goal
     enddo
     pars(i)=pars(goal)
     pars(goal)=par 
  enddo
  do i=1,npar
     order(i)=i
  enddo

end subroutine orderparticles

!-------------------------------------------------------------------------- 
! this routine rearranges the particles in a particle system  according to 
! the array order as in the routine above
!> reorders the particles in a particle system
subroutine orderpsys(psys,order)
  type(particle_system_type), intent(inout) :: psys  !< inout particle sys
#ifdef LONGIDs
  integer(i8b), intent(inout) :: order(1:size(psys%par))  !< desired order
#else
  integer(i4b), intent(inout) :: order(1:size(psys%par))  !< desired order
#endif  
  type(particle_type) :: par
  integer(i8b) :: i
  integer(i8b) :: goal
  
  do i=1,size(psys%par)
     par=psys%par(i)
     goal=order(i)
     do while(goal < i)
        goal=order(goal)
        order(i)=goal
     enddo
     psys%par(i)=psys%par(goal)
     psys%par(goal)=par 
  enddo
  do i=1,size(psys%par)
     order(i)=i
  enddo
  
end subroutine orderpsys

!> creates a transformed particle from an input particle
!============================================================
subroutine copypart(parclone,par,transform)
  type(particle_type), intent(out) :: parclone  !< particle copy
  type(particle_type), intent(in)  :: par       !< input particle
  type(transformation_type), intent(in), optional :: transform !< optional transformation
  parclone = par
  if(present(transform)) parclone%pos = transform%fac * (par%pos - transform%shift)
end subroutine copypart


!> resets the box limits where the BCs are vacuum
!==================================================================
subroutine adjustbox(box,bot,top)
  type(box_type), intent(inout) :: box !< input box
  real(r4b), intent(in) :: bot(3)      !< new bottoms
  real(r4b), intent(in) :: top(3)      !< new tops

  where (box%bbound==0) box%bot = bot
  where (box%tbound==0) box%top = top

end subroutine adjustbox



!> scales particles, sources, and the box from comoving to physical values.
! velocity is taken from Gadget code value to peculiar. 
!==========================================================================
subroutine scale_comoving_to_physical(this, a, h)

!  character(clen), parameter :: myname="scale_comoving_to_physical"
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  type(particle_system_type) :: this
  real(r8b), intent(in) :: a    !< scale factor
  real(r8b), intent(in) :: h    !< hubble parameter (little h)

  call mywrite("   scaling comoving to physical coordinates", verb)
  fmt = "(A,F12.5,T22,A,T25,F12.5)"
  write(str,fmt) "   a = ", a, "h = ", h
  call mywrite(str,verb)

  ! particles
  !------------------------------------------------
  this%par%pos(1) = this%par%pos(1) * a / h
  this%par%pos(2) = this%par%pos(2) * a / h
  this%par%pos(3) = this%par%pos(3) * a / h

#ifdef incVel
  this%par%vel(1) = this%par%vel(1) * sqrt(a) 
  this%par%vel(2) = this%par%vel(2) * sqrt(a) 
  this%par%vel(3) = this%par%vel(3) * sqrt(a) 
#endif

  this%par%mass = this%par%mass / h
  this%par%hsml = this%par%hsml * a / h
  this%par%rho  = ( this%par%rho / (a*a*a) ) * (h*h)


  ! box
  !------------------------------------------------
  this%box%top = this%box%top * a / h
  this%box%bot = this%box%bot * a / h
  

  

end subroutine scale_comoving_to_physical


!> scales particles, sources, and the box from physical to comoving values.
! velocity is taken from peculiar to Gadget code value. 
!==========================================================================
subroutine scale_physical_to_comoving(this, a, h)

  character(clen), parameter :: myname="scale_physical_to_comoving"
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  type(particle_system_type) :: this
  real(r8b), intent(in) :: a   !< scale factor 
  real(r8b), intent(in) :: h   !< hubble parameter (little h)

  call mywrite("   scaling physical to comoving coordinates", verb)
  fmt = "(A,F12.5,T22,A,T25,F12.5)"
  write(str,fmt) "   a = ", a, "h = ", h
  call mywrite(str,verb)

  ! particles
  !------------------------------------------------
  this%par%pos(1) = this%par%pos(1) / a * h
  this%par%pos(2) = this%par%pos(2) / a * h
  this%par%pos(3) = this%par%pos(3) / a * h

#ifdef incVel
  this%par%vel(1) = this%par%vel(1) / sqrt(a)
  this%par%vel(2) = this%par%vel(2) / sqrt(a)
  this%par%vel(3) = this%par%vel(3) / sqrt(a)
#endif

  this%par%mass = this%par%mass * h
  this%par%hsml = this%par%hsml / a * h
  this%par%rho = this%par%rho * (a*a*a) / (h*h)



  ! box
  !------------------------------------------------
  this%box%top = this%box%top / a * h
  this%box%bot = this%box%bot / a * h
  

end subroutine scale_physical_to_comoving


!> set electron fraction from ionization fractions
!=======================================================================================
subroutine set_ye(psys)

  type(particle_system_type) :: psys
  integer(i8b) :: i
#ifdef incHe
  real(r8b) :: Hmf
  real(r8b) :: Hemf
  real(r8b) :: nHe_over_nH
#endif

  do i = 1,size(psys%par)
     psys%par(i)%ye = psys%par(i)%xHII + CV%NeBackground
  end do
  
#ifdef incHe

  do i = 1,size(psys%par)
     
     Hmf  = return_H_mf( psys%par(i)%Hmf )
     Hemf = return_He_mf( psys%par(i)%Hemf ) 
     
     nHe_over_nH = 0.25d0 * Hemf / Hmf
     
     psys%par(i)%ye = psys%par(i)%ye + &
          ( psys%par(i)%xHeII + 2.0d0 * psys%par(i)%xHeIII ) * nHe_over_nH
     
  end do

#endif


end subroutine set_ye










!> sets ionization fractions to their collisional equilibrium values
!=======================================================================================
subroutine set_collisional_ionization_equilibrium(psys, caseA, IsoTemp, DoHydrogen, fit)

  type(particle_system_type) :: psys
  logical, intent(in) :: caseA(2)   
  real(r8b), intent(in) :: IsoTemp
  logical, intent(in) :: DoHydrogen
  character(*), intent(in) :: fit
  real(r8b) :: xvec(5)
  real(r8b) :: Tdum
  integer(i8b) :: i


  ! if we have a single temperature
  !------------------------------------
  if (IsoTemp > 0.0) then
     call calc_colion_eq_fits(fit, IsoTemp, caseA, xvec)
     if (DoHydrogen) then
        psys%par(:)%xHI    = xvec(1)
        psys%par(:)%xHII   = xvec(2)
     endif
#ifdef incHe
     psys%par(:)%xHeI   = xvec(3)
     psys%par(:)%xHeII  = xvec(4)
     psys%par(:)%xHeIII = xvec(5)
#endif

  ! if we have individual temperatures
  !------------------------------------
  else
     do i = 1,size(psys%par)
        Tdum = psys%par(i)%T
        call calc_colion_eq_fits(fit, Tdum, caseA, xvec)
        if (DoHydrogen) then
           psys%par(i)%xHI    = xvec(1)
           psys%par(i)%xHII   = xvec(2)
        endif
#ifdef incHe
        psys%par(i)%xHeI   = xvec(3)
        psys%par(i)%xHeII  = xvec(4)
        psys%par(i)%xHeIII = xvec(5)
#endif
     end do

  end if

end subroutine set_collisional_ionization_equilibrium


!>   enforces a minimum and maximum value on the ionization fractions and temperatures
!=======================================================================================
subroutine enforce_x_and_T_minmax(par,xmin,xmax,tmin,tmax)

  type(particle_type), intent(inout) :: par(:)  !< inout particle system
  real(r8b), intent(in) :: xmin, xmax, tmin, tmax
  integer(i8b) :: i

  do i = 1,size(par)
     if (par(i)%xHI < xmin) par(i)%xHI = xmin
     if (par(i)%xHI > xmax) par(i)%xHI = xmax
#ifdef incHe
     if (par(i)%xHeI  < xmin) par(i)%xHeI  = xmin
     if (par(i)%xHeII < xmin) par(i)%xHeII = xmin
     if (par(i)%xHeI  > xmax) par(i)%xHeI  = xmax
     if (par(i)%xHeII > xmax) par(i)%xHeII = xmax
#endif

     if (par(i)%T < tmin) par(i)%T = tmin
     if (par(i)%T > tmax) par(i)%T = tmax

  end do

end subroutine enforce_x_and_T_minmax





!> outputs currently loaded particle data to the screen
!=================================================================
subroutine particle_info_to_screen(psys,str,lun)

  type(particle_system_type), intent(in) :: psys     !< particle system
  character(*), optional, intent(in) :: str          !< arbitrary string
  integer(i4b), optional, intent(in) :: lun          !< if present goes to file
  integer(i4b) :: outlun
  
  
  outlun=stdout
  if (present(lun)) outlun=lun


  99  format(72("-"))  
  100 format(A,T10,3ES15.5)
!  101 format(A,T10,2I15,ES15.5)
  102 format(A,T10,2I15)
  103 format(A,T10,3I15)


  write(outlun,99) 
  if (present(str)) write(outlun,"(A)") trim(str)
  write(outlun,"(A,I15,A)") "particle data for ", size(psys%par), "  particles"
!  write(outlun,99) 
  write(outlun,*) 


  write(outlun,100) "xpos", minval(psys%par%pos(1)), maxval(psys%par%pos(1)), &
       meanval_real(psys%par%pos(1))

  write(outlun,100) "ypos", minval(psys%par%pos(2)), maxval(psys%par%pos(2)), &
       meanval_real(psys%par%pos(2))
  
  write(outlun,100) "zpos", minval(psys%par%pos(3)), maxval(psys%par%pos(3)), &
       meanval_real(psys%par%pos(3))   
  
#ifdef incVel
  write(outlun,100) "xvel", minval(psys%par%vel(1)), maxval(psys%par%vel(1)), &
       meanval_real(psys%par%vel(1))
  
  write(outlun,100) "yvel", minval(psys%par%vel(2)), maxval(psys%par%vel(2)), &
       meanval_real(psys%par%vel(2))
  
  write(outlun,100) "zvel", minval(psys%par%vel(3)), maxval(psys%par%vel(3)), &
       meanval_real(psys%par%vel(3))
#endif
  
  write(outlun,102) "id",   minval(psys%par%id), maxval(psys%par%id)
  
  write(outlun,100) "mass", minval(psys%par%mass), maxval(psys%par%mass), &
       meanval_real(psys%par%mass)
  
  write(outlun,100) "T",    minval(psys%par(:)%T), maxval(psys%par(:)%T), &
       meanval_real(psys%par%T)
  
  write(outlun,100) "rho",  minval(psys%par%rho), maxval(psys%par%rho), &
       meanval_real(psys%par%rho)
  
  write(outlun,100) "nH",   minval(psys%par%nH), maxval(psys%par%nH), &
       meanval_real(psys%par%nH)

  write(outlun,100) "ye",    minval(psys%par%ye), maxval(psys%par%ye), &
       meanval_real(psys%par%ye)
  
  write(outlun,100) "xHI",    minval(psys%par%xHI), maxval(psys%par%xHI), &
       meanval_real(psys%par%xHI)


#ifdef incCloudy
  write(outlun,100) "xHI_cld",    minval(psys%par%xHI_cloudy), maxval(psys%par%xHI_cloudy), &
       meanval_real(psys%par%xHI_cloudy)
#endif
  


  write(outlun,100) "xHII",   minval(psys%par%xHII), maxval(psys%par%xHII), &
       meanval_real(psys%par%xHII)

#ifdef incHmf  
  write(outlun,100) "Hmf",    minval(psys%par%Hmf), maxval(psys%par%Hmf), &
       meanval_real(psys%par%Hmf)
#endif

  write(outlun,100) "hsml", minval(psys%par%hsml), maxval(psys%par%hsml), &
       meanval_real(psys%par%hsml)

  
#ifdef incH2
  write(outlun,100) "fH2", minval(psys%par%fH2), maxval(psys%par%fH2), &
       meanval_real(psys%par%fH2)
#endif


#ifdef incHe
  write(outlun,100) "xHeI",   minval(psys%par%xHeI), maxval(psys%par%xHeI), &
       meanval_real(psys%par%xHeI)
  
  write(outlun,100) "xHeII",  minval(psys%par%xHeII), maxval(psys%par%xHeII), &
       meanval_real(psys%par%xHeII)
  
  write(outlun,100) "xHeIII", minval(psys%par%xHeIII), maxval(psys%par%xHeIII), &
       meanval_real(psys%par%xHeIII)
#endif


#ifdef incHemf  
  write(outlun,100) "Hemf",    minval(psys%par%Hemf), maxval(psys%par%Hemf), &
       meanval_real(psys%par%Hemf)
#endif


#ifdef outGammaHI
  write(outlun,100) "gammaHI",   minval(psys%par%gammaHI), maxval(psys%par%gammaHI), &
       meanval_real(psys%par%gammaHI)
#endif


#ifdef incEOS
  write(outlun,100) "eos",   minval(psys%par%eos), maxval(psys%par%eos), &
       meanval_real(psys%par%eos)  
#endif


  write(outlun,*)
  
  write(outlun,100) "Box Uppers = ", psys%box%top
  write(outlun,100) "Box Lowers = ", psys%box%bot
  write(outlun,103) "Upr BCs    = ", psys%box%tbound
  write(outlun,103) "Lwr BCs    = ", psys%box%bbound
  
  write(outlun,99) 
  
end subroutine particle_info_to_screen

!> calculates the mean value w/o using the intrinsics
!=================================================================
function meanval_real(arr) result (mean)
  real(r4b), dimension(:), intent(in) :: arr  !< array to average
  real(r8b) :: mean                           !< mean value to return
  integer(i8b) :: i
  
  mean = 0.0d0
  do i = 1,size(arr)
     mean = mean + arr(i)
  end do
  mean = mean / size(arr)

end function meanval_real





end module particle_system_mod
