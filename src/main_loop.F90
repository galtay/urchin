!> \file main_loop.F90

!> \brief The main driver of URCHIN
!! 
!<

module main_loop_mod

  ! utilities
  use myf90_mod
  use omp_lib
  use timers

  ! routines
  use main_input_mod, only: readin_snapshot
  use oct_tree_mod, only: buildtree
  use oct_tree_mod, only: setparticleorder
  use raylist_mod, only: prepare_raysearch
  use output_mod, only: output_total_snap
  use ray_mod, only: make_healpix_ray
  use ray_mod, only: make_probe_ray
  use particle_system_mod, only: particle_info_to_screen
  use particle_system_mod, only: particle_type
  use particle_system_mod, only: create_particle_density_access_list
  use particle_system_mod, only: create_particle_slab_access_list
  use particle_system_mod, only: create_particle_random_access_list

  use mt19937_mod, only: genrand_real1
  use particle_returns_mod, only: return_nH

  use ionpar_mod, only: par_to_tau_par
  use ionpar_mod, only: tauH_from_erf
  use ionpar_mod, only: tauHI_from_erf
  use ionpar_mod, only: tauHIxHI_from_erf
  use ionpar_mod, only: initialize_bckgnd_particle
  use ionpar_mod, only: solve_bckgnd_particle_xH_eq

  use ion_table_class, only: set_mini_spectrum
  use ion_table_class, only: gammaHI_from_mini_spec_thin
  use ion_table_class, only: gammaHI_from_mini_spec_shield

  use raylist_mod, only: trace_ray
  use w4_gadget_spline_kernel_class, only: test_kernels
  use w4_gadget_spline_kernel_class, only: kernel
  use slab_lookup_table_class, only: interpolate_slab_table
  use slab_lookup_table_class, only: return_min_table_vals
  use slab_lookup_table_class, only: return_max_table_vals
  use analytic_ionization_solutions_mod, only: bracket_GHI_plane_parallel
  use hui_gnedin_atomic_rates_mod

  ! variables and types
  use config_mod, only: CV
  use global_mod, only: psys
  use global_mod, only: tree
  use global_mod, only: GV
  use global_mod, only: gconst
  use global_mod, only: saved_gheads
  use global_mod, only: h1_itab

  use ray_mod, only: ray_type
  use ionpar_mod, only: tau_particle_type
  use ionpar_mod, only: bckgnd_particle_type
  use raylist_mod, only: raylist_type
  use ion_table_class, only: mini_spec_type

  implicit none
  
  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: one = 1.0d0
  real(r8b), parameter :: TAU_HI_TH_MAX = 3000.0d0
  real(r8b), parameter :: xHI_hard_pass = 1.0d-6


  real(r8b), parameter :: MIN_GHI_OUTPUT = 1.0d-35
  real(r8b) :: r, ctheta, stheta, phi


  logical :: first_good_par
  integer(i8b) :: first_good_indx
  
  
contains



  
!> this is the main driver of URCHIN
!======================================
subroutine main_loop()    
  character(clen), parameter :: myname="main_loop"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  ! ray tracing
  !---------------
  integer(i4b) :: npix      !< number of Healpix pixels
  integer(i4b) :: nside     !< Healpix nside
  real(r8b) :: sky_frac     !< fraction of sky for each pixel

  !  local counters 
  !-----------------  
  integer(i8b) :: list_indx !< index in particle access list
  integer(i8b) :: psys_indx !< index in particle system 
  integer(i8b) :: ndone     !< particles done this iteration
        
  real(r8b) :: MB
 
  ! work variables
  !-----------------
  integer(i4b) :: ipix

  real(r8b), allocatable :: tauHI_ray_pix(:)
  real(r8b), allocatable :: tauHI_nbr_pix(:)
  real(r8b), allocatable :: tauHI_far_pix(:)
  real(r8b) :: tauHI_far_sum

  real(r8b), allocatable :: tauH_nbr_pix(:)
  real(r8b) :: tauH_nbr_sum

  real(r8b), allocatable :: GHI_ray_pix(:)
  real(r8b), allocatable :: GHI_far_pix(:)
  real(r8b), allocatable :: GHI_use_pix(:)

  real(r8b) :: GHI_far_sum
  real(r8b) :: GHI_use_sum
  real(r8b) :: GHI_int_thin
  real(r8b) :: GHI_total
  



  type(raylist_type), allocatable :: raylist_pix(:)   

  real(r8b) :: slab_table_max_vals(0:3)
  real(r8b) :: slab_table_min_vals(0:3)

  type(mini_spec_type) :: minispec
  type(bckgnd_particle_type) :: bck_par

  real(r8b) :: redshift

  real(r8b) :: tauHI_eff
  real(r8b) :: tauH_nbr_eff

  real(r8b) :: nH, Temp
  real(r8b) :: min_log_tauHI, max_log_tauHI
  real(r8b) :: min_log_K, max_log_K
  real(r8b) :: min_log_tauH, max_log_tauH
  real(r8b) :: min_GHI_pix, max_GHI_pix
  real(r8b) :: min_xHI, max_xHI
  real(r8b) :: max_tauHI_eff

  real(r8b) :: min_logryd
  real(r8b) :: max_logryd
  
  integer(i8b) :: iter
  integer(i8b) :: err





  integer(i8b) :: nthreads
  integer(i8b) :: tid
  character(len=clen) :: msg

  real(r8b) :: nshld

  integer(i8b) :: nskipped_cut
  integer(i8b) :: n_converged
  real(r8b) :: xHInow
  real(r8b) :: xHIold


  logical :: flag_cycle
  logical :: flag_converged_output






  ! initialize variables to track when the first particle is traced,
  ! if we have reached converged output, and min/max values for 
  ! ray traced quantities
  !---------------------------------------------------------------
  first_good_par = .false.
  first_good_indx = 0
  flag_converged_output = .false.

  min_xHI = huge(1.0d0)
  max_xHI = tiny(1.0d0)

  min_log_K = huge(1.0d0)
  max_log_K = tiny(1.0d0)
  
  min_log_tauH = huge(1.0d0)
  max_log_tauH = tiny(1.0d0)       

  min_log_tauHI = huge(1.0d0)
  max_log_tauHI = tiny(1.0d0)       

  GV%cgs_ilen2 = one / (GV%cgs_len * GV%cgs_len)

  ! first timing message
  !---------------------------------------------------------------
  call time_message(stdout," annnnnd they're off ")


  !  set number of pixels and sky fraction from config file values
  !---------------------------------------------------------------- 
  nside = CV%HealpixNside
  select case( CV%FluxDistribution )

  ! isotropic on HealPixels
  !-------------------------------
  case(0,8)
     npix = 12 * nside**2

  ! in + or - z direction
  !-------------------------------
  case(-3,3)
     npix = 1

  ! in + and - z direction
  !-------------------------------
  case(4)
     npix = 2

  ! in +/- x,y,z directions
  !-------------------------------
  case(6)
     npix = 6

  case default
     write(*,*) 'CV%FluxDistribution not recognized'
     write(*,*) CV%FluxDistribution
     stop
  end select

  sky_frac = 1.0d0 / npix



  !  allocate pixel arrays
  !---------------------------------------------------------------- 
  allocate( tauHI_ray_pix (0:npix-1) )
  allocate( tauHI_nbr_pix (0:npix-1) )
  allocate( tauHI_far_pix (0:npix-1) )
  allocate( tauH_nbr_pix  (0:npix-1) )

  allocate( GHI_ray_pix   (0:npix-1) )
  allocate( GHI_far_pix   (0:npix-1) )
  allocate( GHI_use_pix   (0:npix-1) )
  allocate( raylist_pix   (0:npix-1) )



  !  warn if using old RT method
  !---------------------------------------------------------------- 
#ifdef oldRT
  write(*,*) 
  write(*,*) ' *************************** '
  write(*,*) ' *** USING OLD RT METHOD *** '
  write(*,*) ' *************************** '
  write(*,*) 
#endif


  !  read in particle and source snapshot
  !---------------------------------------------------------------- 
  call readin_snapshot()


  !  write timing message to std out
  !---------------------------------------------------------------- 
  call time_message(stdout,' snapshot data read') 


  !  initialize temporary xHI holder
  !---------------------------------------------------------------- 
  psys%par(:)%xHI_tmp = zero

         
  !  build oct tree.  only need to do this once per snap (for now)
  !----------------------------------------------------------------
  call buildtree(psys,tree,MB,CV%PartPerCell)
  GV%MB = GV%MB + MB


  ! re-order particles in memory for efficient access
  !----------------------------------------------------------------
  call setparticleorder(psys, tree)             

    
  ! get values from slab table
  !----------------------------------------------------------------
  call return_max_table_vals( slab_table_max_vals )  
  call return_min_table_vals( slab_table_min_vals )


  ! initialize the npix raylists
  !--------------------------------
  do ipix = 0, npix-1
     call prepare_raysearch( psys, raylist_pix(ipix) )
  end do


  !  set up spectrum 
  !---------------------------------------------------------------- 
  redshift = saved_gheads(0)%z

  min_logryd = log10(CV%SpecLoRyd)
  max_logryd = log10(CV%SpecHiRyd)

  call set_mini_spectrum( min_logryd, max_logryd, redshift, &
       CV%GammaMultiplier, CV%Monochromatic, CV%MonoRydbergs, &
       CV%MonoNormMethod, h1_itab , minispec)

  GHI_int_thin = gammaHI_from_mini_spec_thin( minispec )

  call mywrite('H photoionization rate from HM01 Q+G',verb)
  call mywrite('',verb)
  write(str,'(A,ES15.5)') '   integrating mini spectrum: ', &
       GHI_int_thin
  call mywrite(str,verb)
  write(str,'(A,ES15.5)') '   direct from table:  ', &
       minispec%gammaHI_table
  call mywrite(str,verb)
  write(str,'(A,ES15.5)') '   integrated / table: ', &
       GHI_int_thin/minispec%gammaHI_table
  call mywrite(str,verb)
  write(str,'(A,ES15.5)') '   Gamma HI multiplier: ', &
       CV%GammaMultiplier
  call mywrite('',verb)


  write(*,*) 'Flux Distribution  = ', CV%FluxDistribution
  write(*,*) 'Sky Fraction = ', sky_frac
  write(*,*) 'Photons / cm^3 = ', minispec%ngamma_int

  
  call mywrite('',verb)
  write(str,'(T4,A,ES15.5)') 'maximum ray distance: ', CV%MaxRayDist
  call mywrite(str,verb)
  write(str,'(T4,A,I2)') 'X boundary conditions:  ', CV%XcoordBndryCond
  call mywrite(str,verb)
  write(str,'(T4,A,I2)') 'Y boundary conditions:  ', CV%YcoordBndryCond
  call mywrite(str,verb)
  write(str,'(T4,A,I2)') 'Z boundary conditions:  ', CV%ZcoordBndryCond
  call mywrite(str,verb)
  call mywrite('',verb)

  write(str,'(T4,A,ES12.5)') 'xHI Hard Pass = ', xHI_hard_pass
  call mywrite(str,verb)


  ! stop here and report if just initializing
  !----------------------------------------------------------------
  if(CV%JustInit) then
     write(str,"(A,F10.2)") "total memory allocated [MB] = ", GV%MB
     call mywrite(str,verb)
     call mywrite("just initializing", verb)
     call mywrite("",verb)
     stop
  end if


  call mywrite('beginning iterations ...', verb-1)


 
  ! report number of threads 
  !---------------------------------------------------------------------
!$OMP parallel default(shared) &
!$OMP private(tid)

  tid = omp_get_thread_num()
  if (tid==0) then 
     nthreads = omp_get_num_threads()
     write(stdout,*) '  nthreads = ', nthreads
     write(stdout,*) ''
     call flush(stdout)
  endif

!$OMP end parallel


  

  ! NOW begin real URCHIN run
  !=====================================================================
  over_iterations: do iter = 1, CV%NumIterations

     nskipped_cut = 0
     n_converged = 0
     ndone = 0
     nshld = zero
     max_tauHI_eff = zero



     ! create list to access particles from least to most dense
     !-----------------------------------------------------------
     call create_particle_density_access_list( psys )
!     call create_particle_slab_access_list( psys )

       
     ! begin ray tracing 
     !------------------------- 
     over_pars: do list_indx = 1, size(psys%par)


        ! do some screen reporting
        !-------------------------------
        if (mod(list_indx,size(psys%par)/50)==0 .and. first_good_par) then
           write(*,*) 
           fmt = '(" iter: ", I4.4, ", done ", I12, " of ", I12)'
           write (msg,fmt) iter, list_indx, size(psys%par)
           call time_message(stdout, msg)
           fmt = '(A,I12.12,"   ",I12.12,"   ",I12.12)'
           write(*,fmt) 'nskipped_cut, n_converged, sum: ', &
                nskipped_cut, n_converged, nskipped_cut+n_converged
           fmt = '(A,2ES20.6)'
           write(*,fmt) 'nshld/npar, max_tauHI_eff: ', &
                nshld / size(psys%par), max_tauHI_eff

           fmt = '(A,2ES20.6))'
           if (min_log_K > zero) then
              write(*,fmt) 'min/max log K:    ', log10(min_log_K), log10(max_log_K) 
           else
              write(*,fmt) 'min/max log K:    ', -99.0d0, log10(max_log_K) 
           endif

           if (min_xHI > zero) then
              write(*,fmt) 'min/max log xHI:    ', log10(min_xHI), log10(max_xHI) 
           else
              write(*,fmt) 'min/max log xHI:    ', -99.0d0, log10(max_xHI) 
           endif

           if (min_log_tauH > zero) then
              write(*,fmt) 'min/max log tauH: ', log10(min_log_tauH), log10(max_log_tauH) 
           else
              write(*,fmt) 'min/max log tauH: ', -99.0d0, log10(max_log_tauH) 
           endif

           if (min_log_tauHI > zero) then
              write(*,fmt) 'min/max log tauHI:', log10(min_log_tauHI), log10(max_log_tauHI) 
           else
              write(*,fmt) 'min/max log tauHI:', -99.0d0, log10(max_log_tauHI) 
           endif

        endif


        ! get particle system index
        !-----------------------------------------------------
        psys_indx = psys%acc_list(list_indx)        

        ! skip over highly ionized, optically thin particles
        !-----------------------------------------------------
        flag_cycle = .false.

        if (iter == 1 ) then 
           if ( psys%par(psys_indx)%xHI < xHI_hard_pass * 1.0d-1 ) then
              flag_cycle = .true.
           endif           
        else           
           if ( psys%par(psys_indx)%xHI < xHI_hard_pass ) then
              flag_cycle = .true.
           endif           
        endif


        ! set to equilibrium and skip
        !-----------------------------------------------------
        if (flag_cycle) then

           nskipped_cut = nskipped_cut + 1

           bck_par = initialize_bckgnd_particle( psys%par(psys_indx), &
                                                 GHI_int_thin, &
                                                 tauHI_eff = zero )

           err = solve_bckgnd_particle_xH_eq(bck_par)
           psys%par(psys_indx)%xHI_tmp = psys%par(psys_indx)%xHI 
           psys%par(psys_indx)%xHI     = bck_par%xHI
           psys%par(psys_indx)%xHII    = bck_par%xHII
           psys%par(psys_indx)%gammaHI = GHI_int_thin           

           if ( bck_par%xHI > max_xHI ) max_xHI = bck_par%xHI
           if ( bck_par%xHI < min_xHI ) min_xHI = bck_par%xHI

           cycle

        endif



        ! do some stats on the first good particle
        !---------------------------------------------
        if (.not. first_good_par) first_good_indx = list_indx
        first_good_par = .true.
        
        if (CV%RayStats) then
           if (list_indx == first_good_indx) then
              if (iter==1) then           
                 call mywrite('',verb)
           
                 fmt = '(T3,A)'
                 write(str,fmt) 'Stats for first good particle'
                 call mywrite(str,verb)
                 call mywrite('',verb)
                 
                 fmt = '(T3,A,T15,A,T25,A,T40,A)'
                 write(str,fmt) 'pixel', 'hits', 'tau_ana', 'tau_erf'
                 call mywrite(str,verb)
                 
                 fmt = '(A)'        
                 write(str,fmt) '---------------------------------------------------------------------'
                 call mywrite(str,verb)
              endif
           endif
        endif

        ! set some particle values
        ! at this point all factors of h & a 
        ! should have been accounted for
        !---------------------------------------------
        nH = psys%par(psys_indx)%nH
        Temp = psys%par(psys_indx)%T
        ndone = ndone + 1


!------------------------------------------------------------------------------
!$OMP parallel default(shared) &
!$OMP private(ipix,tid)

        tid = omp_get_thread_num()

!$OMP do
        over_pixels: do ipix = 0, npix-1

           call calc_tau_for_healpixel( psys_indx, ipix, nside, &
                                        raylist_pix(ipix),      &
                                        tauHI_ray_pix(ipix),    &
                                        tauHI_nbr_pix(ipix),    &
                                        tauHI_far_pix(ipix),    &
                                        tauH_nbr_pix(ipix)      )


#ifdef oldRT

           ! If we are using the old simple RT method, use tauHI ray
           !------------------------------------------------------------
           GHI_use_pix(ipix) = &
                gammaHI_from_mini_spec_shield( minispec, tauHI_ray_pix(ipix) ) * &
                sky_frac



#else
           
           ! Alway calculate GHI from far particles
           !------------------------------------------------------------
           GHI_far_pix(ipix) = gammaHI_from_mini_spec_shield( minispec, tauHI_far_pix(ipix) ) 

           ! If we are using monochromatic radiation, use analytic 
           ! slab solution to set the usable GHI
           !------------------------------------------------------------
           if (CV%Monochromatic) then

              ! do angular attenuation before analytic slab solution
              !------------------------------------------------------------
              select case( CV%FluxDistribution )                 
              case(0)
                 ! leave alone
              case(-3,3,4,6)
                 GHI_far_pix(ipix) = GHI_far_pix(ipix) * sky_frac                 
              end select

              ! do analytic slab solution
              !------------------------------------------------------------
              GHI_use_pix(ipix) = &
                   bracket_GHI_plane_parallel( Temp, &
                                               nH, &
                                               GHI_far_pix(ipix), &
                                               tauH_nbr_pix(ipix) * minispec%sigma_ratio(1), &
                                               GHI_int_thin )



              
              ! do angular attenuation after analytic slab solution
              !------------------------------------------------------------
              select case( CV%FluxDistribution )                 
              case(0)
                 GHI_use_pix(ipix) = GHI_use_pix(ipix) * sky_frac                 
              case(-3,3,4,6)
                 ! leave alone
              end select
              

           ! otherwise use lookup table for HM01 to set usable GHI
           !------------------------------------------------------------
           else

              ! If GHI/nH out of range of look up table, dont use it.
              !------------------------------------------------------------
              if ( GHI_far_pix(ipix)/nH > slab_table_max_vals(1) .or. &
                   GHI_far_pix(ipix)/nH < slab_table_min_vals(1)      ) then
                 
                   GHI_use_pix(ipix) = gammaHI_from_mini_spec_shield( &
                        minispec, tauHI_ray_pix(ipix) + tauHI_nbr_pix(ipix) )

                   GHI_use_pix(ipix) = GHI_use_pix(ipix) * sky_frac                 


              ! Otherwise use look up table.
              !------------------------------------------------------------
              else

                 ! do angular attenuation before slab solution lookup
                 !------------------------------------------------------------
                 select case( CV%FluxDistribution )                 
                 case(0)
                    ! leave alone
                 case(-3,3,4,6)
                    GHI_far_pix(ipix) = GHI_far_pix(ipix) * sky_frac                 
                 end select

                 ! do slab solution lookup
                 !------------------------------------------------------------
                 GHI_use_pix(ipix) = &
                      interpolate_slab_table( tauH_nbr_pix(ipix), &
                      GHI_far_pix(ipix)/nH, tauHI_far_pix(ipix), Temp  ) * nH                    
                 
                 ! do angular attenuation after slab solution lookup
                 !------------------------------------------------------------
                 select case( CV%FluxDistribution )                 
                 case(0)
                    GHI_use_pix(ipix) = GHI_use_pix(ipix) * sky_frac                 
                 case(-3,3,4,6)
                    ! leave alone
                 end select

                 
              endif


           endif

#endif           


           if (CV%RayStats) then 
              if (iter == 1) then 
                 if (list_indx == first_good_indx) then
                    call calc_ray_stats(psys_indx, ipix, raylist_pix(ipix), nside )
                 endif
              endif
           endif


        end do over_pixels
!$OMP end do
!$OMP end parallel



        

        ! combine tauHI and GHI before using table
        if (CV%FluxDistribution == 8) then

           GHI_far_sum   = sum( GHI_far_pix )   * sky_frac
           tauH_nbr_sum  = sum( tauH_nbr_pix )  * sky_frac

           tauHI_far_sum = sum( tauHI_far_pix ) * sky_frac

           GHI_use_sum = interpolate_slab_table( tauH_nbr_sum, GHI_far_sum/nH, &
                tauHI_far_sum, Temp  ) * nH                    

           GHI_total = GHI_use_sum

        else

           GHI_total = sum( GHI_use_pix ) 

        endif



        
        ! calculate Hydrogen photoionization rate
        !--------------------------------------------               

        if (GHI_total > GHI_int_thin) then
           write(*,*) 'GHI total > GHI total thin'
           write(*,*) 'GHI_total      = ', GHI_total
           write(*,*) 'GHI_int_thin = ', GHI_int_thin 
           stop
        endif
        
        ! calculate effective taus
        !----------------------------------
        tauH_nbr_eff = sum( tauH_nbr_pix ) * sky_frac
       
        if (GHI_total / GHI_int_thin > zero) then
           tauHI_eff = -log(GHI_total / GHI_int_thin)
        else
           tauHI_eff = TAU_HI_TH_MAX
        endif

        ! calculate mins/maxs
        !-----------------------------------------
        if (tauH_nbr_eff < min_log_tauH) min_log_tauH = tauH_nbr_eff
        if (tauH_nbr_eff > max_log_tauH) max_log_tauH = tauH_nbr_eff

        min_GHI_pix = minval(GHI_use_pix)
        max_GHI_pix = maxval(GHI_use_pix)

        if (min_GHI_pix/nH < min_log_K) min_log_K = min_GHI_pix/nH
        if (max_GHI_pix/nH > max_log_K) max_log_K = max_GHI_pix/nH

        if (tauHI_eff > max_tauHI_eff) max_tauHI_eff = tauHI_eff

        if ( minval(tauHI_ray_pix) < min_log_tauHI ) min_log_tauHI = minval(tauHI_ray_pix)
        if ( maxval(tauHI_ray_pix) > max_log_tauHI ) max_log_tauHI = maxval(tauHI_ray_pix)


        ! initialize background particle. 
        ! here we change the temperature if the particle's effective tauHI_th
        ! is above the value in the config file.
        !-------------------------------------------------------------------
        if (tauHI_eff >= CV%ShldTauThresh) then 
           if (CV%ShldTemp > zero .and. psys%par(psys_indx)%T < 1.0d5) then
              psys%par(psys_indx)%T = CV%ShldTemp
              nshld = nshld + 1.0d0
           endif
           bck_par = initialize_bckgnd_particle( psys%par(psys_indx), &
                                                 GHI_total, &
                                                 tauHI_eff  )
        else
           bck_par = initialize_bckgnd_particle( psys%par(psys_indx), &
                                                 GHI_total, &
                                                 tauHI_eff  )
        endif


        ! update neutral fraction
        !------------------------------------------  
        !------------------------------------------
        err = solve_bckgnd_particle_xH_eq( bck_par )

        !------------------------------------------
        !------------------------------------------

        GV%rayn = GV%rayn + npix 


        ! update and check for convergence
        !------------------------------------------
        if ( bck_par%xHI > max_xHI ) max_xHI = bck_par%xHI
        if ( bck_par%xHI < min_xHI ) min_xHI = bck_par%xHI


        psys%par(psys_indx)%xHI     = bck_par%xHI
        psys%par(psys_indx)%xHII    = bck_par%xHII
        psys%par(psys_indx)%gammaHI = max( GHI_total, MIN_GHI_OUTPUT )

        xHInow = psys%par(psys_indx)%xHI
        xHIold = psys%par(psys_indx)%xHI_tmp

        if ( xHInow - xHIold < CV%ConvTol * xHIold ) then
           n_converged = n_converged + 1
        endif

        psys%par(psys_indx)%xHI_tmp = psys%par(psys_indx)%xHI 


!        if (ndone == 100000) stop 'profiling'




        
     end do over_pars
     
     
     call particle_info_to_screen(psys)
     call flush(stdout)

     ! if iteration condition met or converged 
     !----------------------------------------------------
     if ( iter == 1    .or. iter == 3    .or.          &
          iter == 10   .or. iter == 20   .or.          & 
          iter == 40   .or. iter == 80   .or.          &
          iter == 160  .or. iter == 320  .or.          &
          nskipped_cut + n_converged == size(psys%par) ) then

        GV%OutputIndx = iter

        ! if converged relabel output 999 
        !----------------------------------------------------
        if( nskipped_cut + n_converged == size(psys%par) ) then
           GV%OutputIndx = 999 
        endif

        call output_total_snap(psys)      

     endif


     if ( nskipped_cut + n_converged == size(psys%par) ) then

        write(*,*)
        write(*,*) ' *** Convergence Reached on Iteration: ', iter
        write(*,*) 
        stop

     endif

     
  end do over_iterations


  close(GV%pardatalun)
  
  
  
end subroutine main_loop


! seperate out the calculation of tau for each direction
!-----------------------------------------------------------
subroutine calc_tau_for_healpixel(psys_indx, ipix, nside, raylist, &
     tauHI_ray, tauHI_nbr, tauHI_far, tauH_nbr )

  integer(i8b), intent(in) :: psys_indx          !< index of active particle in psys
  integer(i4b), intent(in) :: ipix               !< pixel ID
  integer(i4b), intent(in) :: nside              !< Healpix Nside
  type(raylist_type), intent(inout) :: raylist   !< ray/intersection object

  real(r8b), intent(out) :: tauHI_ray            !< tauHI along whole ray
  real(r8b), intent(out) :: tauHI_nbr            !< tauHI from neighbor particles
  real(r8b), intent(out) :: tauHI_far            !< tauHI from far particles
  real(r8b), intent(out) :: tauH_nbr             !< tauH from neighbor particles
 
  real(r8b) :: start(3)                          !< start position of ray
  real(r8b) :: vec(3)                            !< unit direction of ray

  integer(i8b) :: ihit                           !< intersection counter
  integer(i8b) :: hit_indx                       !< index of hit particle in psys

  type(ray_type) :: ray
  type(particle_type) :: hitpar
  type(tau_particle_type) :: tpar



  ! actions always performed
  !=======================================================================================

     
  ! Error function integrations thru partial particles. 
  !-------------------------------------------------------------------------------
  start = psys%par(psys_indx)%pos


  ! Set rays based on flux distribution
  !-------------------------------------------------------------------------------
  select case( CV%FluxDistribution )

  case(0)
     if (CV%MaxRayDist > zero) then
        call make_healpix_ray( start, nside, ipix, ray, length=CV%MaxRayDist )
     else
        call make_healpix_ray( start, nside, ipix, ray)
     endif

  case(-3)
     vec = (/0.0d0, 0.0d0, -1.0d0/)
     if (CV%MaxRayDist > zero) then
        call make_probe_ray( start, vec, ray, length=CV%MaxRayDist )
     else
        call make_probe_ray( start, vec, ray)
     endif

  case(3)
     vec = (/0.0d0, 0.0d0, 1.0d0/)
     if (CV%MaxRayDist > zero) then
        call make_probe_ray( start, vec, ray, length=CV%MaxRayDist )
     else
        call make_probe_ray( start, vec, ray)
     endif

  case(4)
     select case ( ipix )
     case(0)
        vec = (/0.0d0, 0.0d0, -1.0d0/)
     case(1)
        vec = (/0.0d0, 0.0d0,  1.0d0/)        
     case default
        write(*,*) 'ipix out of range'
        write(*,*) 'ipix = ', ipix
        stop
     end select

     if (CV%MaxRayDist > zero) then
        call make_probe_ray( start, vec, ray, length=CV%MaxRayDist )
     else
        call make_probe_ray( start, vec, ray)
     endif

  case(6)
     select case ( ipix )
     case(0)
        vec = (/-1.0d0, 0.0d0, 0.0d0/)
     case(1)
        vec = (/1.0d0, 0.0d0,  0.0d0/)        
     case(2)
        vec = (/0.0d0, -1.0d0, 0.0d0/)
     case(3)
        vec = (/0.0d0, 1.0d0,  0.0d0/)        
     case(4)
        vec = (/0.0d0, 0.0d0, -1.0d0/)
     case(5)
        vec = (/0.0d0, 0.0d0,  1.0d0/)        
     case default
        write(*,*) 'ipix out of range'
        write(*,*) 'ipix = ', ipix
        stop
     end select

     if (CV%MaxRayDist > zero) then
        call make_probe_ray( start, vec, ray, length=CV%MaxRayDist )
     else
        call make_probe_ray( start, vec, ray)
     endif


  end select
  

  ! trace ray
  !-------------------------------------------------------------------------------  
  raylist%ray = ray
  call trace_ray(ray, raylist, psys, tree, dosort=.false.) 
  
  tauHI_ray = zero
  tauHI_nbr = zero
  tauHI_far = zero
  tauH_nbr = zero  


  ! OLD RT method considers all particles when 
  ! calculating tauHI and hence GHI
  !---------------------------------------------------------------
#ifdef oldRT
  do ihit = 1, raylist%nnb
     hit_indx = raylist%intersection(ihit)%pindx
     hitpar = psys%par(hit_indx)
     tpar = par_to_tau_par( hitpar, raylist%intersection(ihit), raylen=ray%length )

     call tauHI_from_erf( tpar, GV%cgs_ilen2 )
     tauHI_ray = tauHI_ray + tpar%tau_erf
  end do

!  GHI_ray = gammaHI_from_mini_spec_shield( minispec, tauHI_ray ) 

  ! NEW RT method considers only non-neighbor particles 
  ! when calculating tauHI and hence GHI then uses look up table
  !---------------------------------------------------------------
#else
  do ihit = 1, raylist%nnb    
     hit_indx = raylist%intersection(ihit)%pindx
     hitpar   = psys%par(hit_indx)
     tpar     = par_to_tau_par( hitpar, raylist%intersection(ihit), raylen=ray%length )

     if (raylist%intersection(ihit)%s .le. hitpar%hsml .or. hit_indx == psys_indx) then
        call tauH_from_erf( tpar, GV%cgs_ilen2 )   
        tauH_nbr = tauH_nbr + tpar%tau_erf
        call tauHI_from_erf( tpar, GV%cgs_ilen2 ) 
        tauHI_nbr = tauHI_nbr + tpar%tau_erf
     else
        call tauHI_from_erf( tpar, GV%cgs_ilen2 ) 
        tauHI_far = tauHI_far + tpar%tau_erf
     endif
  enddo

  tauHI_ray = tauHI_nbr + tauHI_far


#endif


end subroutine calc_tau_for_healpixel





! calculate some stats on the first iteration
!-----------------------------------------------------------
subroutine calc_ray_stats(psys_indx, ipix, raylist, nside)

  character(clen), parameter :: myname="calc_ray_stats"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  integer(i8b), intent(in) :: psys_indx
  integer(i4b), intent(in) :: ipix
  type(raylist_type), intent(inout) :: raylist
  integer(i4b), intent(in) :: nside

  integer(i8b) :: tid

  character(clen) :: rayfile
  integer(i8b) :: raystatlun
  integer(i8b) :: ihit
  integer(i8b) :: hit_indx

  real(r8b) :: start(3)
  real(r8b) :: ray_len

  type(ray_type) :: ray
  type(particle_type) :: hitpar
  type(tau_particle_type) :: tpar

  real(r8b) :: tauHI_ray_erf
  real(r8b) :: tauHI_ray_ana

  integer(i8b) :: Nsegments
  integer(i8b) :: iseg
  real(r8b) :: ds
  real(r8b) :: pos(3)
  real(r8b) :: rr, hh
  real(r8b) :: nH, nHI
  !real(r8b) :: T, y

  character(clen) :: dum(5)


  ! actions performed on the first iteration and the first good particle if ray stats is on
  !=======================================================================================

        
  tid = omp_get_thread_num()
  raystatlun = 510 + tid
  
  write(rayfile,'(A,I4.4,A)') trim(CV%OutputDir) // '/rayxyz', ipix, '.txt'
  open(unit=raystatlun, file=rayfile, action="write")
  
  do ihit = 1,raylist%nnb
     hit_indx = raylist%intersection(ihit)%pindx
     if (hit_indx == psys_indx) cycle
     write(raystatlun,'(4ES15.5)') psys%par(hit_indx)%pos, raylist%intersection(ihit)%d
  end do
  close(raystatlun)
  
  
  ! do a test against a spec wizard style numerical integration thru the first rays
  !---------------------------------------------------------------------------------
  
  
  
  ! error function method
  !--------------------------
  start = psys%par(psys_indx)%pos
  ray_len = CV%MaxRayDist
  
  if (CV%MaxRayDist > zero) then
     call make_healpix_ray( start, nside, ipix, ray, length=ray_len )
  else
     call make_healpix_ray( start, nside, ipix, ray)
  endif
  raylist%ray = ray
  call trace_ray(ray, raylist, psys, tree, dosort=.false.) 
  
  tauHI_ray_erf = zero
  do ihit = 1, raylist%nnb
     hit_indx = raylist%intersection(ihit)%pindx
     hitpar = psys%par(hit_indx)
     tpar = par_to_tau_par( hitpar, raylist%intersection(ihit), raylen=ray%length )  
     call tauHI_from_erf( tpar, GV%cgs_ilen2 )
     tauHI_ray_erf = tauHI_ray_erf + tpar%tau_erf
  end do
  
  
  ! spec wizard SPH method
  !--------------------------
  Nsegments = 1000
  ds = ray_len / Nsegments
  tauHI_ray_ana = 0.0d0
  do iseg = 0, Nsegments-1
     nH = 0.0d0
     nHI = 0.0d0
     pos = start + ray%dir * (iseg+0.5d0)*ds
     
     do ihit = 1, raylist%nnb
        hit_indx = raylist%intersection(ihit)%pindx
        hitpar = psys%par(hit_indx)
        tpar = par_to_tau_par( hitpar, raylist%intersection(ihit), raylen=ray%length ) 
        rr = sqrt(  sum( (hitpar%pos - pos)**2 ) ) * GV%cgs_len
        hh = hitpar%hsml * GV%cgs_len
        if (rr<hh) then
           nH  = nH  + kernel(rr,hh) * tpar%Hcnt
           nHI = nHI + kernel(rr,hh) * tpar%HIcnt
        endif
     end do
     tauHI_ray_ana = tauHI_ray_ana + nHI * tpar%sigmaHIth * ds * GV%cgs_len
     
  end do
  
  
  write(dum(1),'(I4.4)') ipix
  write(dum(2),'(I8.8)') raylist%nnb
  write(dum(3),'(ES13.4)') tauHI_ray_ana
  write(dum(4),'(ES13.4)') tauHI_ray_erf
  
  
  
  fmt = '(T3,A, T15, A, T25, A, T40, A)'
  write(str,fmt) adjustl(trim(dum(1))), adjustl(trim(dum(2))), adjustl(trim(dum(3))), &
       adjustl(trim(dum(4)))
  call mywrite(str,verb)
  



end subroutine calc_ray_stats










!!$! compute SPH density using provided smoothing lengths
!!$! THIS MAY HAVE A BUG - NEED TO CHECK
!!$!--------------------------------------------------------
!!$subroutine compute_density(ngas, part)
!!$  implicit none
!!$  integer, intent(in) :: ngas
!!$  type(particle_type), intent(inout) :: part(ngas)
!!$  ! local variables
!!$  integer :: ip, ic, jp, iter
!!$  integer :: nlist(1000), nb
!!$  real (r8b) :: dist1, h, h2, hinv, xpart(3), rho
!!$  integer, parameter :: nmin=30,nmax=50
!!$
!!$  write(*,*) 'ngas/100: ', ngas/100
!!$
!!$  do ip=1, ngas
!!$
!!$     if (mod(ip,ngas/100)==0) then
!!$        write(*,*) 'ip,nb:', ip, nb
!!$     endif
!!$
!!$     nb    = 0
!!$     rho   = zero
!!$     !
!!$     h     = part(ip)%hsml
!!$     h2    = h*h
!!$     hinv  = 1.d0/h
!!$     xpart = part(ip)%pos
!!$     iter  = 0
!!$     !
!!$     do while (nb <= nmin .or. nb >= nmax) 
!!$        nb    = 0
!!$        rho   = zero
!!$        do jp=1, ngas
!!$           if (ip /= jp) then
!!$              ic      = 1
!!$              dist1    = abs(part(jp)%pos(ic)-xpart(ic))
!!$              if(dist1 < h) then
!!$                 ic = 2
!!$                 dist1 = dist1 * dist1 + (part(jp)%pos(ic)-xpart(ic)) * (part(jp)%pos(ic)-xpart(ic))
!!$                 if (dist1 < h2) then
!!$                    ic = 3
!!$                    dist1 = dist1 + (part(jp)%pos(ic)-xpart(ic)) * (part(jp)%pos(ic)-xpart(ic))
!!$                    if (dist1 < h2) then
!!$                       nb        = nb + 1
!!$!                       nlist(nb) = jp
!!$                       rho       = rho + part(jp)%mass * kernel(sqrt(dist1), h)
!!$                    endif
!!$                 endif
!!$              endif
!!$           endif
!!$        end do
!!$!        write (*,*) ip,nb,h
!!$        if(nb <= nmin) then
!!$           h     = h * 1.3
!!$        endif
!!$        if(nb >= nmax) then
!!$           h = h / 1.26
!!$        endif
!!$        h2    = h*h
!!$        hinv  = 1.d0/h
!!$        iter  = iter + 1
!!$        !
!!$        if (iter > 300) then
!!$           write (*,*) ' not converging: ',h,nb
!!$           stop
!!$        endif
!!$        
!!$     end do
!!$
!!$     part(ip)%rho = rho
!!$     part(ip)%hsml = h
!!$
!!$  enddo
!!$
!!$end subroutine compute_density




end module main_loop_mod
