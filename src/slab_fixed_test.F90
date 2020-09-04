! runs slab fixed test and quits
!----------------------------------


module slab_fixed_test
use myf90_mod
use config_mod, only: CV
use global_mod, only: psys, tree, GV
use mt19937_mod, only: init_mersenne_twister
use particle_system_mod, only: calc_bytes_per_particle
use gadget_public_input_hdf5_mod, only: get_planning_data_gadget_public_hdf5
use gadget_public_input_hdf5_mod, only: read_Gpubhdf5_particles
use oct_tree_mod, only: buildtree
use oct_tree_mod, only: setparticleorder
use raylist_mod, only: prepare_raysearch
use raylist_mod, only: trace_ray
use ray_mod, only: make_probe_ray
use ionpar_mod, only: par_to_tau_par
use ionpar_mod, only: tauHI_from_erf
use w4_gadget_spline_kernel_class, only: kernel
use particle_system_mod, only: create_particle_slab_access_list

use ray_mod, only: ray_type
use raylist_mod, only: raylist_type
use ionpar_mod, only: tau_particle_type
use gadget_general_class, only: gadget_constants_type
use particle_system_mod, only: particle_type
implicit none
private


public :: run_slab_fixed_test


real(r8b) :: one = 1.0d0


contains





subroutine run_slab_fixed_test()
  real(r8b) :: MB
  integer(i8b) :: i, j, ihit, ni, nj, hit_indx, Nsegments, iseg, nnb
  
  
  type(ray_type) :: ray
  type(raylist_type) :: raylist
  type(particle_type) :: hitpar
  type(tau_particle_type) :: tpar
  type(gadget_constants_type) :: gconst

  real(r8b), parameter :: sigmaHIth = 6.3462964d-18 ! verner  
  real(r8b), parameter :: xHI_slab_fixed_test = 5.0d-2
  real(r8b), parameter :: F0 = 1.0d4

  real(r8b) :: total_mass, mean_rho, mean_nH, min_nH, max_nH
  real(r8b) :: wall_padding
  
  character(clen) :: log_file
  integer(i4b) :: log_lun

  real(r8b) :: path_len, ray_len, h2, h, r2, r

  real(r8b) :: nH, nHI, tausum_ana, tausum_erf
  
  real(r8b) :: di, dj, ds

  real(r8b) :: start(3), dir(3), pos(3)

  real(r8b) :: z_ceiling

  integer(i8b) :: list_indx, ipar, ndone

  real(r8b) :: GHI_ana, GHI_erf

  character(clen) :: fmt

  write(*,*) 
  write(*,*) 'running slab_fixed test ...'
  write(*,*) 


  ! check Makefile options
  !------------------------------
#ifdef incHmf
  write(*,*) 'incHmf should not be defined in Makefile for "slab_fixed_test"'
  stop
#endif
  
#ifdef incH2
  write(*,*) 'incH2 should not be defined in Makefile for "slab_fixed_test"'
  stop
#endif

  ! set global variables
  !------------------------------
  CV%OutputDir = '../data/test_output/slab_fixed_test'
  CV%SnapNum = 7
  CV%ParFilesPerSnap = 1
  CV%SnapPath = '../data/density_fields/slab'
  CV%ParFileBase = 'glass-h-N32'
  CV%PartPerCell = 8
  CV%H_mf = 1.0d0

  CV%Verbosity = 3
  myf90_verbosity = CV%Verbosity

  ! do input
  !------------------------------
  call calc_bytes_per_particle(GV%bytesperpar)
  call get_planning_data_gadget_public_hdf5()   
  call read_Gpubhdf5_particles()
  GV%cgs_ilen2 = one / (GV%cgs_len * GV%cgs_len)

  !  make the z boundaries vacuum
  !----------------------------------
  psys%box%tbound(3) = 0
  psys%box%bbound(3) = 0

  !  make the x and y boundaries periodic
  !----------------------------------
  psys%box%tbound(1) = 1
  psys%box%bbound(1) = 1
  psys%box%tbound(2) = 1
  psys%box%bbound(2) = 1

  ! set the box dimensions
  !-----------------------------------------
  psys%box%bot = 0.0d0
  psys%box%top = GV%BoxUprsComoh(1)

  ! make tree
  !-----------------------------------------
  call buildtree(psys, tree, MB, CV%PartPerCell)
  call setparticleorder(psys, tree)             
  call prepare_raysearch(psys, raylist)


  ! set test ionization states
  !-----------------------------
  psys%par(:)%xHI = xHI_slab_fixed_test
  psys%par(:)%xHII = 1.0d0 - psys%par(:)%xHI
  psys%par(:)%ye = psys%par(:)%xHII  



  ! do density calculations
  ! this sort of mean rho is biased towards high values. 
  !-----------------------------
  total_mass = 0.0d0
  do i = 1, size(psys%par)
     total_mass = total_mass + psys%par(i)%mass
  end do
  mean_rho = total_mass / GV%BoxUprsComoh(1)**3
  mean_nH  = mean_rho * GV%cgs_rho / gconst%ProtonMass

  min_nH = minval( psys%par(:)%rho ) * GV%cgs_rho / gconst%ProtonMass
  max_nH = maxval( psys%par(:)%rho ) * GV%cgs_rho / gconst%ProtonMass

  write(*,*) 'min/max nH: ', min_nH, max_nH
  write(*,*) 'mean nH:    ', mean_nH


  ! set ray wall_padding (just to stay away from edges)
  !-------------------------------------------------
  wall_padding = maxval( psys%par(:)%hsml )
  z_ceiling = GV%BoxUprsComoh(1) - wall_padding
  write(*,*) 
  write(*,*) 'wall_padding = ', wall_padding

  ! set log file
  !-------------------------------------------------
  log_file = '../data/test_output/slab_fixed_test/tau_GHI_plane_N32.txt'
  write(*,*) 'log file: ', trim(log_file)
  call open_formatted_file_w( log_file, log_lun )
  write(log_lun,'(A,T20,A,T35,A,T50,A,T65,A)') '#   path length', 'tau analytic', 'GHI analytic', &
       'tau erf', 'GHI erf' 



  ! set particle access
  !-------------------------------------------------
  call create_particle_slab_access_list( psys )




  ndone = 0

  do list_indx = 1, size(psys%par)
        
     ipar = psys%acc_list(list_indx)  

     if ( psys%par(ipar)%pos(3) >= z_ceiling    ) cycle
     if ( psys%par(ipar)%pos(3) <= wall_padding ) cycle

     ndone = ndone + 1

     path_len = z_ceiling - psys%par(ipar)%pos(3) 
     ray_len = path_len 
          
     start = psys%par(ipar)%pos        
     dir = (/ 0.0, 0.0, 1.0 /)
     
     ! Do error function solution
     !---------------------------------------------------
     call make_probe_ray( start, dir, ray, length=ray_len )
     raylist%ray = ray
     call trace_ray(ray, raylist, psys, tree, dosort=.false.) 

     ! Calculate tau
     !----------------------------------------------------------     
     tausum_erf = 0.0d0


     do ihit = 1, raylist%nnb    
        hit_indx = raylist%intersection(ihit)%pindx
        hitpar = psys%par(hit_indx)
        tpar = par_to_tau_par( hitpar, raylist%intersection(ihit), raylen=ray%length )
!        if (raylist%intersection(ihit)%s .lt. hitpar%hsml) cycle
        call tauHI_from_erf( tpar, GV%cgs_ilen2 )                           
        tausum_erf = tausum_erf + tpar%tau_erf                        
     enddo



     ! Calculate analytic solution and write to file
     !---------------------------------------------------
     tausum_ana = path_len * GV%cgs_len * mean_nH * xHI_slab_fixed_test * sigmaHIth     
     
     GHI_ana = F0 * sigmaHIth * exp(-tausum_ana)
     GHI_erf = F0 * sigmaHIth * exp(-tausum_erf)
 
     fmt = "(5ES15.5)"
     write(log_lun,fmt) path_len, tausum_ana, GHI_ana, tausum_erf, GHI_erf
     
  enddo
  
  close(log_lun)
       
  write(*,*) 'finished slab_fixed test'   
  stop
  


  
  
  


end subroutine run_slab_fixed_test



end module slab_fixed_test
