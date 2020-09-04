!> \file ray.F90

!> \brief the ray module 
!!
!<

module ray_mod
use myf90_mod
use particle_system_mod
use oct_tree_mod
use mt19937_mod, only: genrand_real1
use pix_tools, only: pix2vec_ring
implicit none
private

public :: ray_type
public :: cell_intersection
public :: part_intersection
public :: transform_ray
public :: dist2ray
public :: make_healpix_ray
public :: make_probe_ray

real(r8b), parameter :: one = 1.0d0
real(r8b), parameter :: zero = 0.0d0
  

!> a ray to measure optical depth along path
!-----------------------------------------------------------------------
  type ray_type
     integer(i4b) :: class  !< based on direction signs (MMM, PMM, ...)
     real(r8b) :: start(3)  !< starting position of active region
     real(r8b) :: dir(3)    !< UNIT vector for direction
     real(r8b) :: length    !< length of ray
  end type ray_type



contains


!> creates a healpix ray
!-----------------------------------------------------------------------  
  subroutine make_healpix_ray(start, nside, ipix, ray, length)

    real(r8b) :: start(3)                     !< the ray starting position
    integer(i4b), intent(in) :: nside         !< the healpix nsides
    integer(i4b), intent(in) :: ipix          !< the healpix pixel number
    type(ray_type), intent(out) :: ray        !< output ray
    real(r8b), intent(in), optional :: length !< optional active length, default=huge

    call pix2vec_ring( nside, ipix, ray%dir )
    ray%start = start 

    if (present(length)) then
       ray%length = length
    else
       ray%length = huge(1.0d0) * 0.1d0
    endif

    ray%class = ray_class_from_dir( ray%dir )

  end subroutine make_healpix_ray





!> creates probe ray
!-----------------------------------------------------------------------  
  subroutine make_probe_ray(pos, dir, ray, length)

    real(r8b), intent(in) :: pos(3)           !< ray start
    real(r8b), intent(in) :: dir(3)           !< ray start
    type(ray_type), intent(out) :: ray        !< output ray
    real(r8b), intent(in), optional :: length !< optional active length, default=huge

    ray%start = pos 
    ray%dir = dir

    if (present(length)) then
       ray%length = length
    else
       ray%length = huge(1.0d0) * 0.1d0
    endif

    ray%class = ray_class_from_dir( ray%dir )

  end subroutine make_probe_ray




  ! pre computes the class of the ray for the Pluecker test
  ! ray label    class
  !   MMM          0
  !   PMM          1
  !   MPM          2
  !   PPM          3
  !   MMP          4
  !   PMP          5
  !   MPP          6
  !   PPP          7
  !-----------------------------------------------------------
  function ray_class_from_dir( dir ) result( class ) 
    real(r8b) :: dir(3)
    real(i4b) :: class
    integer(i4b) :: i

    class = 0
    do i = 1, 3
       if ( dir(i) >= zero ) class = class + 2**(i-1)
    end do

  end function ray_class_from_dir





!> properly sets the starting point, direction, class, and length of a ray
!--------------------------------------------------------------------------
  subroutine set_ray(ray, start, dir, len)
    type(ray_type) :: ray     !< inout ray
    real(r8b) :: start(3)     !< starting position
    real(r8b) :: dir(3)       !< direction 
    real(r8b),optional :: len !< length    
    real(r8b) :: unit(3)

    unit = dir / sqrt( dot_product(dir,dir) )

    ray%start = start
    ray%dir   = unit
    ray%class = ray_class_from_dir( ray%dir )

    if (present(len)) then
       ray%length = len
    else
       ray%length = huge(1.0d0) * 0.1d0
    endif


  end subroutine set_ray

  
!> returns a transformed ray while preserving the initial ray
!--------------------------------------------------------------
  subroutine transform_ray(ray, inray, pm)
    type(ray_type) :: ray           !< output transformed ray
    type(ray_type) :: inray         !< input ray
    type(transformation_type) :: pm !< transformation

    ray%start  = inray%start * pm%fac + pm%shift
    ray%dir    = inray%dir   * pm%fac
    ray%length = inray%length  
    ray%class  = ray_class_from_dir( ray%dir )

  end subroutine transform_ray


!> returns the minimum (dperp) distance between a point and a ray,
!! the projected distance (dproj), and the distance from the ray 
!! start to the point (dstrt).
!!
!!  ray start
!!  |
!!  >---------------------------------->
!!                               |
!!  |----------------------------| 
!!             dproj             | dperp
!!                               |
!!                        pos -  *
!!
!--------------------------------------------------------------  
  subroutine dist2ray(ray, pos, dperp, dproj, dstrt) 
    type(ray_type), intent(in) :: ray    !< input ray
    real(r4b), intent(in) :: pos(3)      !< input position
    real(r8b), intent(out) :: dperp      !< perpendicular distance
    real(r8b), intent(out) :: dproj      !< projected distance 
    real(r8b), intent(out) :: dstrt      !< distance to ray start

    real(r8b) :: vec1(3)
    real(r8b) :: vec2(3)
    real(r8b) :: vperp(3)

    ! vec1 goes from the ray start to pos
    vec1 = pos - ray%start 

    ! vec2 is unit vector in ray direction
    vec2 = ray%dir

    ! dproj is the projected distance from ray start
    ! to where dperp intersects the ray
    dproj = dot_product( vec1, vec2 )

    ! dstart is the distace from ray start to pos
    dstrt = sqrt( dot_product( vec1, vec1 ) )

    ! vperp is the minimum length vector from 
    ! pos to the ray
    vperp  = vec1 - dproj * ray%dir 

    ! dperp is the norm of vperp 
    dperp = sqrt( dot_product( vperp, vperp ) )

  end subroutine dist2ray


!> returns the minimum distance between a particle and a ray
!--------------------------------------------------------------  
  subroutine pdist2ray(ray, part, dperp, dproj, dstrt) 
    type(ray_type), intent(in) :: ray        !< input ray
    type(particle_type), intent(in) :: part  !< input particle
    real(r8b), intent(out) :: dperp          !< perpendicular distance
    real(r8b), intent(out) :: dproj          !< projected distance 
    real(r8b), intent(out) :: dstrt          !< distance to ray start
    real(r4b) :: pos(3)

    pos = part%pos
    call dist2ray( ray, pos, dperp, dproj, dstrt )

  end subroutine pdist2ray



  
!> tests for ray / particle intersection. 
!----------------------------------------------  
  subroutine part_intersection(ray, part, dperp, dproj, dstrt, hit) 
  type(ray_type), intent(in) :: ray
  type(particle_type), intent(in) :: part   !< particle
  real(r8b), intent(out) :: dperp           !< perpendicular distance to particle 
  real(r8b), intent(out) :: dproj           !< projected distance along ray
  real(r8b), intent(out) :: dstrt           !< distance to ray start
  logical, intent(out) :: hit               !< true or false result
  
  real(r8b) :: end2cen         !< distance^2 from ray end to part position
  real(r8b) :: diff(3)         !< vector from ray end to part center

  
  ! calculate distances from particle to ray
  !----------------------------------------------------
  call dist2ray( ray, part%pos, dperp, dproj, dstrt )

  
  ! if the perpendicular distance to the point is larger than hsml exit
  !---------------------------------------------------------------------
  if (dperp > part%hsml) then
     hit = .false.
     return
  endif
  
  ! now our only concern is the position of the particle
  ! along the ray.  
  !--------------------------------------------------------


  ! accept particles that must be intersections
  !----------------------------------------------------------------
  if (dproj > part%hsml .or. dproj < ray%length - part%hsml) then
     hit = .true.
     return
  endif


  ! reject particles that cannot possibly be intersections
  !----------------------------------------------------------------
  if ( dproj < -part%hsml ) then
     hit = .false. 
     return
  endif

  if ( dproj > ray%length + part%hsml ) then
     hit = .false. 
     return
  endif



  ! if the particle smoothing volume contains the origin of the ray
  !---------------------------------------------------------------------
  if ( dstrt <= part%hsml ) then
     hit = .true.
     return
  endif

  ! if the particle smoothing volume could contain the terminus of 
  ! the ray, test for that.
  !---------------------------------------------------------------------
  if ( abs(dproj - ray%length) < part%hsml ) then
     diff = part%pos - (ray%start + ray%length * ray%dir)
     end2cen = sum( diff*diff )
     if (end2cen <= part%hsml*part%hsml) then
        hit = .true.
        return
     endif
  endif

  ! at this point we cant  have a hit
  !-------------------------------------
  hit = .false.
  return


  end subroutine part_intersection
  


!> tests for ray / cell intersection. 
!----------------------------------------------  
  function cell_intersection(ray,cell) result(hit)
    logical :: hit           !< true or false result
    type(ray_type) :: ray    !< ray
    type(cell_type) :: cell  !< cell
    real(r8b) bot(3),top(3)
    bot = cell%botrange - ray%start
    top = cell%toprange - ray%start
    hit = pluecker(ray, bot, top)
  end function cell_intersection
  

!> pluecker test for line segment / cell intersection
!-----------------------------------------------------    
  function pluecker(ray, s2b, s2t) result( hit )

    logical :: hit            !< true or false result
    type(ray_type) :: ray     !< input ray
    real(r8b) :: s2b(3)       !< vector from ray start to lower cell corner
    real(r8b) :: s2t(3)       !< vector from ray start to upper cell corner

    real(r8b) :: dir(3)
    real(r8b) :: dist

    real(r8b) :: e2b(3)       !< vector from ray end to lower cell corner
    real(r8b) :: e2t(3)       !< vector from ray end to upper cell corner

    dir  = ray%dir  
    dist = ray%length

    e2b = s2b - dir * dist
    e2t = s2t - dir * dist

    hit = .false.

    ! branch on ray direction
    !---------------------------
    select case( ray%class )

    ! MMM
    !-----------
    case(0)

       if(s2b(1) > zero .or. s2b(2) > zero .or. s2b(3) > zero) return ! on negative part of ray 
       if(e2t(1) < zero .or. e2t(2) < zero .or. e2t(3) < zero) return ! past length of ray      

       if ( dir(1)*s2b(2) - dir(2)*s2t(1) < zero .or.  &
            dir(1)*s2t(2) - dir(2)*s2b(1) > zero .or.  &
            dir(1)*s2t(3) - dir(3)*s2b(1) > zero .or.  &
            dir(1)*s2b(3) - dir(3)*s2t(1) < zero .or.  &
            dir(2)*s2b(3) - dir(3)*s2t(2) < zero .or.  &
            dir(2)*s2t(3) - dir(3)*s2b(2) > zero       ) return

    ! PMM
    !-----------
    case(1)

       if(s2t(1) < zero .or. s2b(2) > zero .or. s2b(3) > zero) return ! on negative part of ray 
       if(e2b(1) > zero .or. e2t(2) < zero .or. e2t(3) < zero) return ! past length of ray      

       if ( dir(1)*s2t(2) - dir(2)*s2t(1) < zero .or.  &
            dir(1)*s2b(2) - dir(2)*s2b(1) > zero .or.  &
            dir(1)*s2b(3) - dir(3)*s2b(1) > zero .or.  &
            dir(1)*s2t(3) - dir(3)*s2t(1) < zero .or.  &
            dir(2)*s2b(3) - dir(3)*s2t(2) < zero .or.  &
            dir(2)*s2t(3) - dir(3)*s2b(2) > zero       ) return

    ! MPM
    !-----------
    case(2)

       if(s2b(1) > zero .or. s2t(2) < zero .or. s2b(3) > zero) return ! on negative part of ray 
       if(e2t(1) < zero .or. e2b(2) > zero .or. e2t(3) < zero) return ! past length of ray      

       if ( dir(1)*s2b(2) - dir(2)*s2b(1) < zero .or.  &
            dir(1)*s2t(2) - dir(2)*s2t(1) > zero .or.  &
            dir(1)*s2t(3) - dir(3)*s2b(1) > zero .or.  &
            dir(1)*s2b(3) - dir(3)*s2t(1) < zero .or.  &
            dir(2)*s2t(3) - dir(3)*s2t(2) < zero .or.  &
            dir(2)*s2b(3) - dir(3)*s2b(2) > zero       ) return

    ! PPM
    !-----------
    case(3)

       if(s2t(1) < zero .or. s2t(2) < zero .or. s2b(3) > zero) return ! on negative part of ray 
       if(e2b(1) > zero .or. e2b(2) > zero .or. e2t(3) < zero) return ! past length of ray      

       if ( dir(1)*s2t(2) - dir(2)*s2b(1) < zero .or.  &
            dir(1)*s2b(2) - dir(2)*s2t(1) > zero .or.  &
            dir(1)*s2b(3) - dir(3)*s2b(1) > zero .or.  &
            dir(1)*s2t(3) - dir(3)*s2t(1) < zero .or.  &
            dir(2)*s2t(3) - dir(3)*s2t(2) < zero .or.  &
            dir(2)*s2b(3) - dir(3)*s2b(2) > zero       ) return

    ! MMP
    !-----------
    case(4)

       if(s2b(1) > zero .or. s2b(2) > zero .or. s2t(3) < zero) return ! on negative part of ray 
       if(e2t(1) < zero .or. e2t(2) < zero .or. e2b(3) > zero) return ! past length of ray      

       if ( dir(1)*s2b(2) - dir(2)*s2t(1) < zero .or.  &
            dir(1)*s2t(2) - dir(2)*s2b(1) > zero .or.  &
            dir(1)*s2t(3) - dir(3)*s2t(1) > zero .or.  &
            dir(1)*s2b(3) - dir(3)*s2b(1) < zero .or.  &
            dir(2)*s2b(3) - dir(3)*s2b(2) < zero .or.  &
            dir(2)*s2t(3) - dir(3)*s2t(2) > zero       ) return


    ! PMP
    !-----------
    case(5)

       if(s2t(1) < zero .or. s2b(2) > zero .or. s2t(3) < zero) return ! on negative part of ray 
       if(e2b(1) > zero .or. e2t(2) < zero .or. e2b(3) > zero) return ! past length of ray      

       if ( dir(1)*s2t(2) - dir(2)*s2t(1) < zero .or.  &
            dir(1)*s2b(2) - dir(2)*s2b(1) > zero .or.  &
            dir(1)*s2b(3) - dir(3)*s2t(1) > zero .or.  &
            dir(1)*s2t(3) - dir(3)*s2b(1) < zero .or.  &
            dir(2)*s2b(3) - dir(3)*s2b(2) < zero .or.  &
            dir(2)*s2t(3) - dir(3)*s2t(2) > zero       ) return


    ! MPP
    !-----------
    case(6)

       if(s2b(1) > zero .or. s2t(2) < zero .or. s2t(3) < zero) return ! on negative part of ray 
       if(e2t(1) < zero .or. e2b(2) > zero .or. e2b(3) > zero) return ! past length of ray      

       if ( dir(1)*s2b(2) - dir(2)*s2b(1) < zero .or.  &
            dir(1)*s2t(2) - dir(2)*s2t(1) > zero .or.  &
            dir(1)*s2t(3) - dir(3)*s2t(1) > zero .or.  &
            dir(1)*s2b(3) - dir(3)*s2b(1) < zero .or.  &
            dir(2)*s2t(3) - dir(3)*s2b(2) < zero .or.  &
            dir(2)*s2b(3) - dir(3)*s2t(2) > zero       ) return

    ! PPP
    !-----------
    case(7)

       if(s2t(1) < zero .or. s2t(2) < zero .or. s2t(3) < zero) return ! on negative part of ray 
       if(e2b(1) > zero .or. e2b(2) > zero .or. e2b(3) > zero) return ! past length of ray      

       if ( dir(1)*s2t(2) - dir(2)*s2b(1) < zero .or.  &
            dir(1)*s2b(2) - dir(2)*s2t(1) > zero .or.  &
            dir(1)*s2b(3) - dir(3)*s2t(1) > zero .or.  &
            dir(1)*s2t(3) - dir(3)*s2b(1) < zero .or.  &
            dir(2)*s2t(3) - dir(3)*s2b(2) < zero .or.  &
            dir(2)*s2b(3) - dir(3)*s2t(2) > zero       ) return

    case default
       call rayError('ray class.')

    end select

    hit=.true.

  end function pluecker







!> error handling
!-----------------------------      
  subroutine rayError(string,i)
    character(*) :: string  !< error string
    integer, optional :: i  !< error number
    
    print*,' Error detected:'
    
    if(present(i)) then
       print*,string,i
    else
       print*,string
    endif
    
    stop
  end subroutine rayError
  
end module ray_mod
