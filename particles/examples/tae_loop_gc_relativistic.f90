program tae_loop

use particle_tracer
use mod_particle_diagnostics
use mpi
use mod_atomic_elements
use mod_particle_io
use mod_event
use mod_project_particles
use mod_particle_loop
use mod_jorek_timestepping
use mod_random_seed
use mod_interp, only: mode_moivre, interp_RZ, interp_0
use mod_basisfunctions
use nodes_elements
use phys_module, only: tstep, restart, t_start, restart_particles, nout_projection, nout
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint, amin, F0, R_geo
use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0, T_EP_eV, RHO_EP, Lambda_peak, delta_Lambda, output_pe_E_mu, output_pe_mu_Pphi, use_pe_EmuPphi
use phys_module, only: n_mode_families
use phys_module, only: particle_pusher

use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG

use mod_particle_sputtering, only: particle_sputter, sample_fluid_particle_energy
use mod_projection_functions, only: proj_f_combined_density, &
                                    proj_f_combined_energy, proj_f_combined_par_momentum
use mod_edge_domain
use mod_edge_elements
use data_structure, only: type_bnd_element_list, type_bnd_node_list 
use equil_info
use mod_boundary, only: boundary_from_grid
use mod_phase_space_project
use mod_projection_functions_phase
use mod_math_operators , only: cross_product
use mod_tae_loop_init_callbacks

!$ use omp_lib

implicit none

type(event)                                       :: fieldreader, partreader, partwriter
!type(adf11_all)                                   :: adas
type(pcg32_rng), dimension(:), allocatable        :: rng
type(pcg32_rng)                                   :: rng_part_init
type(count_action)                                :: counter
type(projection), target                          :: jorek_feedback, project_density, project_current
type(jorek_timestep_action), target               :: jorek_stepper
type(type_edge_domain), allocatable, dimension(:) :: edge_domains
type(edge_elements)                               :: D_edge
type(write_particle_diagnostics)                  :: diag
type(phase_space_projection)                      :: fourD_dist, power_exchange_vpar_mu, power_exchange_vpar_mu_trap, power_exchange_vpar_mu_pass, power_exchange_Emu, power_exchange_Emu_trap, power_exchange_Emu_pass, power_exchange_EmuPphi, power_exchange_EmuPphi_trap, power_exchange_EmuPphi_pass, RZ_dist, PparMu_dist, psiN_dist, sqrtPsiN_dist, KE_dist, res_num_trap, res_num_pass, EPphi_dist, EPphiMu_dist, EMu_dist

real*8, parameter  :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)
real*8    :: target_time
real*8    :: physical_particles, weight, tstep_keep
real*8    :: oldtime, step_rest_time, particle_step_time, particle_start_time, diag_time, theta_part
real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, tstep_si, timesteps
real*8    :: v_kin_temp, E(3), B(3), psi, U, B_norm(3), Rbox(2), Zbox(2)
real*8    :: rescale_coef, T_axis(1), E_axis, E_hot, v2, tstart_jorek
real*8    :: dummy_real8_1, dummy_real8_2, dummy_real8_3, kT_mc2_ratio, p_star, gamma_star
real*8, allocatable :: real_gdf_param(:), real_pdf_param(:), psi_minmax_loc(:), phase_space_bounds_RZPhiP(:,:), phase_space_bounds_EPphi(:,:), phase_space_bounds_psiNThetaPhiE(:,:), phase_space_bounds_psiNThetaPhiP(:,:),  uniform_var(:)
integer, allocatable :: int_param_sampler(:)

! CHANGE THIS FOR DIFFERENT SIMULATION SETUPS
! ---- Wang2020-style setup (centrally peaked beta profile + isotropic Maxwellian) ----
! Distribution: f(R,Z,p) = exp(-(r/a/0.4)^2) * f_MJ(p), sampled in (R,Z,phi,p),
! isotropic momentum split beta = 2u-u^2, random phi in [0,2pi), isotropic pressure feedback.
! NOTE: T_EP_eV and RHO_EP are NOT hardcoded here - they are read from the namelist (&in1).
real*8, parameter :: &
                      A_MINOR          = 0.85          , &    ! in meter
                      Poloidalbound(2) = [0.d0, TWOPI], &
                      Phibound(2)      = [0.d0, TWOPI], &
                      ! Pbound(2)        = [1.5d4, 2.5d7],&  ! [AMU·m/s] for deuterium ions at T_EP = 100 keV
                      Pbound(2)        = [1.d2, 1.2d6], &    ! [AMU·m/s] relativistic momentum bound
                      vparbound(2)     = [-SPEED_OF_LIGHT, SPEED_OF_LIGHT] ! [m/s] parallel velocity bound
integer, parameter :: &
                      MASS_IDX              = -1,  & ! MASS INDEX for EP
                      ATOMIC_NUMBER         = -1,  & ! ATOMIC NUMBER for EP
                      CHARGE_EP             = -1,  & ! EP charge in units of elementary charge
                      ! N_PHI_PLANES          = 12   ! particles are spread over TWOPI * i/N_PHI_PLANES in phi to reduce noise when initial mode structure is used
                      N_PHI_PLANES          = 1  ! so particles' phi are uniformly random in [0, 2*PI)
! ---- Orbit-class enum (replaces bare 0/1/-1 magic numbers) ----
integer, parameter :: ORBIT_UNKNOWN =  0  !< classification not yet determined
integer, parameter :: ORBIT_TRAPPED =  1  !< banana (trapped) orbit
integer, parameter :: ORBIT_PASSING = -1  !< transit (passing) orbit
real   , parameter :: TAE_FREQ      = 300d3  ! Precomputed analytical TAE frequency [Hz]
integer, parameter :: N_TOR_RESONANCE      = 4   !< toroidal mode number n in resonance condition
integer, parameter :: MIN_PERIODS_RESONANCE = 5  !< minimum completed periods before computing L
logical, parameter :: &
                      use_CGL_pressure                      = .false.,   & ! .true.: full CGL anisotropic pressure tensor (needs use_pcs_full=.t.); .false.: isotropic (1/3)p^2 diagonal
                      use_trap_passing_for_PE               = .true.    ! whether to distinguish the power exchange due to trapped particles and that due to passing particles
integer   :: n_real_gdf_param, n_real_pdf_param, n_variables, n_int_param_sampler, RZ_OUTPUT_STEPS, POWER_EXCHANGE_STEPS, PARTICLE_OUTPUT_STEPS
integer   :: n_particles_local,n_reflect,ifail,dummy_int, ino
integer   :: i, j, k, l, m, n_steps, i_elm_old, i_diagno
integer   :: seed, i_rng, n_stream
character(len=50) :: part_output_filename
real*8    :: t_loop_start, t_loop_end, t_barrier_start, t_loop_time

! ---- Persistent per-particle state for trapped/passing classification ----
! theta_prev_arr  : poloidal angle at the end of the last call [0, 2pi)
! theta_accum_arr : net signed angle accumulated since the last completed orbit
! dtheta_prev_arr : last angular step (carries sign; 0 at initialisation skips reversal check)
! trap_pass_arr   : ORBIT_UNKNOWN / ORBIT_TRAPPED / ORBIT_PASSING
real*8,  allocatable :: theta_prev_arr(:)
real*8,  allocatable :: theta_accum_arr(:)
real*8,  allocatable :: dtheta_prev_arr(:)
integer, allocatable :: trap_pass_arr(:)
! ---- Additional persistent state for omega_theta / omega_phi computation ----
real*8,  allocatable :: phi_prev_arr(:)          ! raw phi at last step
real*8,  allocatable :: phi_unwrap_arr(:)        ! monotonically accumulating (unwrapped) toroidal angle
real*8,  allocatable :: t_ref_a_arr(:)           ! time of last upper turn (trapped) / crossing (passing)
real*8,  allocatable :: phi_ref_a_arr(:)         ! phi_unwrap at ref_a
real*8,  allocatable :: t_ref_b_arr(:)           ! time of last lower turn (trapped only)
real*8,  allocatable :: phi_ref_b_arr(:)         ! phi_unwrap at ref_b
integer, allocatable :: n_periods_arr(:)          ! completed period count
real*8,  allocatable :: omega_theta_buf_arr(:,:) ! circular buffer of per-period omega_theta (MIN_PERIODS_RESONANCE, n_particles)
real*8,  allocatable :: omega_phi_buf_arr(:,:)   ! circular buffer of per-period omega_phi   (MIN_PERIODS_RESONANCE, n_particles)

! Start up MPI, jorek
call sim%initialize(num_groups=1)

#ifndef __NVCOMPILER
if (sim%my_id .eq. 0)  write(*,*) '__NVCOMPILER not defined'
#endif

n_particles_local = int(n_particles/sim%n_cpu) 
timesteps         = tstep_particles
tstep_keep       = tstep

! Ensure a sane projection output cadence:
! - nout_projection presets to -1; fall back to nout (mirrors mod_project_particles) so the
!   cadence is correct even when it is not set in the namelist.
! - Clamp the step counts to >= 1 so mod(ino, POWER_EXCHANGE_STEPS) can never divide by zero.
if (nout_projection .le. 0) nout_projection = nout
POWER_EXCHANGE_STEPS =  max(nout_projection*2, 1)
RZ_OUTPUT_STEPS      =  max(nout_projection*5, 1)
PARTICLE_OUTPUT_STEPS = max(nout*50, 1)

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1)) ! read jorek_restart.h5 
! fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1, mode_divisor=1)) ! read jorek_restart.h5 and do not reduce its mode amplitudes
! fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1,mode_divisor=200)) ! read jorek_restart.h5 and reduce its mode amplitude
call with(sim, fieldreader)
! Check that the mesh element_list pointer is associated (must be after fieldreader loads the mesh)
if (.not. associated(sim%fields%element_list)) then
  if (sim%my_id .eq. 0) then
    write(*,*) 'ERROR: sim%fields%element_list is not associated after reading fields!'
    write(*,*) 'Aborting due to invalid mesh pointer.'
  endif
  call MPI_ABORT(MPI_COMM_WORLD, 2, 2)
end if

write(*,*) 'Finished reading JOREK fields. n_elements = ', sim%fields%element_list%n_elements

tstep = tstep_keep
write(*,*) 'main : t_start = ',t_start

if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)

call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

call update_equil_state(sim%my_id,sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)
if(sim%my_id .eq. 0) write(*,*) "ES%Psi_axis, ES%Psi_bnd:", ES%Psi_axis, ES%Psi_bnd

!! Compute minimum and maximum psi on each element, for particle initalisation later
allocate(psi_minmax_loc(2*sim%fields%element_list%n_elements))
if(sim%my_id .eq. 0) call extract_element_psi_minmax(sim%fields,psi_minmax_loc)
call MPI_BCAST(psi_minmax_loc, 2*sim%fields%element_list%n_elements, MPI_REAL8, 0, MPI_COMM_WORLD, ifail)
if(ifail .ne. MPI_SUCCESS) then
  write(*,*) 'ERROR: MPI_BCAST failed with code', ifail
  call MPI_ABORT(MPI_COMM_WORLD, ifail, ifail)
endif

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = central_mass * ATOMIC_MASS_UNIT * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek

tstep_si  = tstep * t_norm
n_steps   = floor(tstep_si / timesteps)
timesteps = tstep_si / n_steps
n_steps   = tstep_si / timesteps

if (sim%my_id .eq.0) then
  write(*,*) ' Relativistic TAE loop with particle ', trim(particle_pusher), ' model'
  write(*,*) ' Pressure projection, use_CGL_pressure = ', use_CGL_pressure
  write(*,*) ' adapt time step to be multiple of jorek time step'
  write(*,*) "tstep = ", tstep_si, n_steps, timesteps
  write(*,*) "check :", n_steps, tstep_si - n_steps*timesteps

  i_diagno =  sim%fields%node_list%n_nodes / 3 
  write(*,'(A,6f8.4)') ' probe at : ',sim%fields%node_list%node(i_diagno)%x(1,1,1:2)
  open(111,file='diagno.txt')
endif

if (.not. restart_particles) then
  ! Set up particles
  sim%groups(1)%Z    = ATOMIC_NUMBER
  sim%groups(1)%mass = atomic_weights(MASS_IDX) !< atomic mass units
  ! Set up parameters for sampling
  n_variables = 5 ! psiN, theta, phi, p, u_beta
  kT_mc2_ratio = T_EP_eV * EL_CHG / (sim%groups(1)%mass * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2)
  
  n_real_gdf_param   = 2*sim%fields%element_list%n_elements+3
  allocate(real_gdf_param(n_real_gdf_param))
  real_gdf_param(1:2*sim%fields%element_list%n_elements) = psi_minmax_loc
  
  ! if (sim%my_id .eq. 0) then
  !   call find_axis(sim%my_id, sim%fields%node_list, sim%fields%element_list, dummy_real8_1, &
  !   real_gdf_param(2*sim%fields%element_list%n_elements+1), & ! R_axis
  !   real_gdf_param(2*sim%fields%element_list%n_elements+2), & ! Z_axis
  !   dummy_int,dummy_real8_2,dummy_real8_3,ifail)
  ! endif
  ! call MPI_BCAST(real_gdf_param(2*sim%fields%element_list%n_elements+1), 2, MPI_REAL8, 0, MPI_COMM_WORLD, ifail)
  
  real_gdf_param(2*sim%fields%element_list%n_elements+1) = ES%R_axis
  real_gdf_param(2*sim%fields%element_list%n_elements+2) = ES%Z_axis
  real_gdf_param(2*sim%fields%element_list%n_elements+3) = sim%groups(1)%mass ! in [AMU]
  
  if (sim%my_id .eq. 0) write(*,*) "R_axis, Z_axis:", real_gdf_param(2*sim%fields%element_list%n_elements+1), real_gdf_param(2*sim%fields%element_list%n_elements+2)
  
  n_int_param_sampler = 1
  allocate(int_param_sampler(n_int_param_sampler))
  int_param_sampler(1) = CHARGE_EP  ! charge in units of elementary charge
  
  n_real_pdf_param = 7
  allocate(real_pdf_param(n_real_pdf_param))
  real_pdf_param(1) = kT_mc2_ratio
  real_pdf_param(2) = T_EP_eV ! in eV
  real_pdf_param(3) = sim%groups(1)%mass ! in [AMU]

  call p_star_gamma0(kT_mc2_ratio, sim%groups(1)%mass, p_star, gamma_star)
  real_pdf_param(4) = 2*log(p_star) - gamma_star/kT_mc2_ratio ! logfmax for the f_{MJ}/f_max for the rejection sampling
  real_pdf_param(5) = ES%R_axis
  real_pdf_param(6) = ES%Z_axis
  real_pdf_param(7) = amin

  
  call domain_bounding_box(sim%fields%node_list, sim%fields%element_list, Rbox(1), Rbox(2), Zbox(1), Zbox(2)) ! get the min and max of R and Z for particle initialisation
  
  select case(trim(particle_pusher))
  case('gc_rel')
    allocate(particle_gc_relativistic::sim%groups(1)%particles(n_particles_local))
  case('kinetic_rel')
    allocate(particle_kinetic_relativistic::sim%groups(1)%particles(n_particles_local))
  case default
    if (sim%my_id .eq. 0) write(*,*) 'ERROR: unknown particle_pusher = ', trim(particle_pusher)
    call MPI_ABORT(MPI_COMM_WORLD, 3, 3)
  end select

  ! Debug: print temperature diagnostics
  if (sim%my_id .eq. 0) then
    write(*,*) '=== Temperature Diagnostics ==='
    write(*,*) 'T_EP_eV (from namelist) = ', T_EP_eV
    write(*,*) 'kT_mc2_ratio = ', kT_mc2_ratio
    write(*,*) 'mass_AMU = ', sim%groups(1)%mass
    write(*,*) 'Lambda_peak = ', Lambda_peak
    write(*,*) 'delta_Lambda = ', delta_Lambda
    write(*,*) 'F0 = ', F0, ', R_geo = ', R_geo, ', |B0| = |F0|/R_geo = ', abs(F0)/R_geo
    write(*,*) 'ES%R_axis = ', ES%R_axis, ', ES%Psi_axis = ', ES%Psi_axis, ', ES%Psi_bnd = ', ES%Psi_bnd
    write(*,*) '=============================='
  endif
  
  ! CHANGE THIS FOR DIFFERENT INITIAL DISTRIBUTIONS:
  ! 1. gdf_sampler -> reject sampling with pdf -> repeat until all particles are sampled
  ! 2. convert with sampler to particle function  
    
  ! ********** WANG2020 in flux coordinates: f(psiN,p) = exp(-psiN/0.4^2) * f_MJ(p)/sup(f_MJ) **********
  ! r/a ~ sqrt(psiN) so the wang2020 spatial profile exp(-(r/a/0.4)^2) becomes exp(-psiN/0.16).
  ! NOTE: PsiNThetaPhiP_to_gc_relativistic_cb uses the isotropic split only when Lambda_peak <= 0.01,
  !       so the namelist MUST keep Lambda_peak = -1 (isotropic) for the wang2020 setup.
  allocate(phase_space_bounds_psiNThetaPhiP(5,2))
  phase_space_bounds_psiNThetaPhiP(:, 1) = [0.01d0, Poloidalbound(1), Phibound(1), Pbound(1), 0.d0] ! min for psiN, theta, phi, p [AMU·m/s] , u_beta
  phase_space_bounds_psiNThetaPhiP(:, 2) = [0.99d0, Poloidalbound(2), Phibound(2), Pbound(2), 1.d0] ! max for psiN, theta, phi, p [AMU·m/s] , u_beta

  select case(trim(particle_pusher))
  case('gc_rel')
    call initialise_particles_in_phase_space(n_variables, sim%groups(1)%particles, sim%fields, pcg32_rng(), pdf_psiN_p_wang2020, &
    particle_weight_one, gdf_uniform, gdf_PsiNThetaPhi_EorP_sampler, 1.001d0, 1.0d0, &
    PsiNThetaPhiP_to_gc_relativistic_cb, sim%groups(1)%mass, sim%time, phase_space_bounds_psiNThetaPhiP, &
    n_real_pdf_param_in=n_real_pdf_param, real_pdf_param_in=real_pdf_param,& ! parameters for pdf
    n_real_gdf_param_in=n_real_gdf_param, real_gdf_param_in=real_gdf_param,& ! parameters for gdf sampler
    n_real_samp_to_part_param_in=n_real_gdf_param, real_samp_to_part_param_in=real_gdf_param, n_int_samp_to_part_param_in=n_int_param_sampler, int_samp_to_part_param_in=int_param_sampler & ! parameters for sample_to_gc_relativistic
    )
  case('kinetic_rel')
    call initialise_particles_in_phase_space(n_variables, sim%groups(1)%particles, sim%fields, pcg32_rng(), pdf_psiN_p_wang2020, &
    particle_weight_one, gdf_uniform, gdf_PsiNThetaPhi_EorP_sampler, 1.001d0, 1.0d0, &
    PsiNThetaPhiP_to_kinetic_relativistic_cb, sim%groups(1)%mass, sim%time, phase_space_bounds_psiNThetaPhiP, &
    n_real_pdf_param_in=n_real_pdf_param, real_pdf_param_in=real_pdf_param,& ! parameters for pdf
    n_real_gdf_param_in=n_real_gdf_param, real_gdf_param_in=real_gdf_param,& ! parameters for gdf sampler
    n_real_samp_to_part_param_in=n_real_gdf_param, real_samp_to_part_param_in=real_gdf_param, n_int_samp_to_part_param_in=n_int_param_sampler, int_samp_to_part_param_in=int_param_sampler & ! parameters for sample_to_kinetic_relativistic
    )
  end select
  ! ********************************************************************************

  ! ********** distribution f(E, P_phi) as function of invariants via rejection sampling in (R,Z,phi,ppar,mu) this distr df/dt = 0 for a few hundreds of Alfven times
  
  ! allocate(phase_space_bounds_EPphi(5, 2))
  ! phase_space_bounds_EPphi(:,1) = [Rbox(1), Zbox(1), 0.d0,      -5d5, 0.0d0]     ! min for R,Z,phi, ppar [AMU·m/s], mu [AMU·m²/(T·s²)]
  ! phase_space_bounds_EPphi(:,2) = [Rbox(2), Zbox(2), TWOPI,      5d5, 9.d13]    ! max

  ! ! Allocate a temporary 4-element array because Intel Fortran does not allow
  ! ! passing an array section to an allocatable dummy argument.
  ! block
  !   real*8, allocatable :: pdf_param_ep(:)
  !   allocate(pdf_param_ep(4))
  !   pdf_param_ep = real_pdf_param(1:4)
  !   call initialise_particles_in_phase_space(5, sim%groups(1)%particles, sim%fields, pcg32_rng(), pdf_E_Pphi, &
  !     particle_weight_one, gdf_uniform, gdf_RZPhiPparMu_sampler, 1.001d0, 1.0d0, &
  !     RZPhiPparMu_to_gc_relativistic_cb, sim%groups(1)%mass, sim%time, phase_space_bounds_EPphi, &
  !     n_real_pdf_param_in=4, real_pdf_param_in=pdf_param_ep, &
  !     n_int_pdf_param_in=n_int_param_sampler, int_pdf_param_in=int_param_sampler, &
  !     n_int_samp_to_part_param_in=n_int_param_sampler, int_samp_to_part_param_in=int_param_sampler &
  !     )
  !   deallocate(pdf_param_ep)
  ! end block

  ! deallocate(phase_space_bounds_EPphi)
! ********************************************************************************

  call adjust_particle_weights(sim%groups(1)%particles, RHO_EP)
  if (sim%my_id .eq. 0) write(*,*) "Particle density was adjusted to:", RHO_EP, sim%groups(1)%particles(1:10)%weight

  
  allocate(uniform_var(1))
  call rng_part_init%initialize(1, random_seed(), sim%n_cpu, sim%my_id+1)
  if(N_PHI_PLANES .gt. 1) then
    ! Randomly spread particles over N_PHI_PLANES toroidal planes to reduce noise from initial mode structure
    ! while avoiding correlations between phi and particle ordering in memory
    do i=1,size(sim%groups(1)%particles,1)
      call rng_part_init%next(out=uniform_var)
      j = floor(uniform_var(1) * N_PHI_PLANES)
      ! j = mod(i-1, N_PHI_PLANES)
      sim%groups(1)%particles(i)%x(3) = TWOPI * real(j,8) / real(N_PHI_PLANES,8)
    enddo
  else 
    ! Randomly initialise particles' phi in [0, 2*PI)
    do i=1,size(sim%groups(1)%particles,1)
      call rng_part_init%next(out=uniform_var)
      sim%groups(1)%particles(i)%x(3) = TWOPI * uniform_var(1)
    enddo
  endif 
  deallocate(uniform_var)

  do i=1,min(300, size(sim%groups(1)%particles,1))
    ! print particle's info for debugging
    theta_part = atan2(sim%groups(1)%particles(i)%x(2) - ES%Z_axis, sim%groups(1)%particles(i)%x(1) - ES%R_axis)
    if (sim%my_id .eq. 0) write(*,*) ' particle ', i, ' phi = ', sim%groups(1)%particles(i)%x(3), ', theta = ', theta_part
  enddo

  ! Debug: print sample particle energies to verify the distribution
  if (sim%my_id .eq. 0) then
    write(*,*) '=== Sample particle energies (first 20) ==='
    do i = 1, min(20, size(sim%groups(1)%particles, 1))
      associate(p => sim%groups(1)%particles(i))
        select type(p)
        type is (particle_gc_relativistic)
          write(*,'(A,I4,A,ES12.4, A,ES12.4, A,ES12.4)') &
            '  p(', i, '): p_par=', p%p(1), ' mu=', p%p(2), ' weight=', p%weight
        type is (particle_kinetic_relativistic)
          write(*,'(A,I4,A,3ES12.4, A,ES12.4)') &
            '  p(', i, '): px,py,pz=', p%p(1), p%p(2), p%p(3), ' weight=', p%weight
        end select
      end associate
    end do
    write(*,*) '========================================='
  endif

  ! Write initial particle distribution to file
  partwriter = event(write_action(filename='part_restart.h5'))
  call with(sim, partwriter)
  
else  ! restarting particles

  if (sim%my_id .eq. 0) write(*,*) 'restarting particles: reading part_restart.h5'

  deallocate(sim%groups)
  allocate(sim%groups(0))

  partreader = event(read_action(filename='part_restart.h5'))
  call with(sim, partreader)

  ! Guard: verify that part_restart.h5 actually contains a particle group
  if (size(sim%groups,1) .lt. 1) then
    if (sim%my_id .eq. 0) write(*,*) 'ERROR: part_restart.h5 contains no particle groups'
    call MPI_ABORT(MPI_COMM_WORLD, 4, 4)
  endif

  ! Guard: verify that the particle type read from part_restart.h5 is consistent
  ! with the particle_pusher specified in the namelist. If not, abort the run.
  select type (p => sim%groups(1)%particles)
  type is (particle_gc_relativistic)
    if (trim(particle_pusher) .ne. 'gc_rel') then
      if (sim%my_id .eq. 0) write(*,*) &
        'ERROR: part_restart.h5 contains particle_gc_relativistic, but particle_pusher = ', trim(particle_pusher)
      call MPI_ABORT(MPI_COMM_WORLD, 4, 4)
    endif
  type is (particle_kinetic_relativistic)
    if (trim(particle_pusher) .ne. 'kinetic_rel') then
      if (sim%my_id .eq. 0) write(*,*) &
        'ERROR: part_restart.h5 contains particle_kinetic_relativistic, but particle_pusher = ', trim(particle_pusher)
      call MPI_ABORT(MPI_COMM_WORLD, 4, 4)
    endif
  class default
    if (sim%my_id .eq. 0) write(*,*) &
      'ERROR: part_restart.h5 contains an unsupported particle type'
    call MPI_ABORT(MPI_COMM_WORLD, 4, 4)
  end select

endif

! Initialise persistent trap/pass arrays for restarted particles
allocate(theta_prev_arr(n_particles_local),  source=0.d0)
allocate(theta_accum_arr(n_particles_local), source=0.d0)
allocate(dtheta_prev_arr(n_particles_local), source=0.d0)
allocate(trap_pass_arr(n_particles_local),   source=ORBIT_UNKNOWN)
allocate(phi_prev_arr(n_particles_local),        source=0.d0)
allocate(phi_unwrap_arr(n_particles_local),      source=0.d0)
allocate(t_ref_a_arr(n_particles_local),         source=0.d0)
allocate(phi_ref_a_arr(n_particles_local),       source=0.d0)
allocate(t_ref_b_arr(n_particles_local),         source=0.d0)
allocate(phi_ref_b_arr(n_particles_local),       source=0.d0)
allocate(n_periods_arr(n_particles_local),       source=0)
allocate(omega_theta_buf_arr(MIN_PERIODS_RESONANCE, n_particles_local)); omega_theta_buf_arr = 0.d0
allocate(omega_phi_buf_arr(MIN_PERIODS_RESONANCE, n_particles_local));   omega_phi_buf_arr   = 0.d0

do i = 1, n_particles_local
  theta_prev_arr(i) = atan2(sim%groups(1)%particles(i)%x(2) - ES%Z_axis, sim%groups(1)%particles(i)%x(1) - ES%R_axis)
  phi_prev_arr(i)   = modulo(sim%groups(1)%particles(i)%x(3), TWOPI)
  phi_unwrap_arr(i) = modulo(sim%groups(1)%particles(i)%x(3), TWOPI)
enddo

jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                filter    = filter_perp, filter_hyper    = filter_hyper, filter_parallel    = filter_par, &
                                filter_n0 = filter_perp, filter_hyper_n0 = filter_hyper, filter_parallel_n0 = filter_par_n0, &
                                calc_integrals=.false., to_vtk=.false., to_h5 = .false., basename='projections')

!Full pressure tensor has 6 elements due to symmetry.
!Components 1..6 in the order expected by model307 use_pcs_full:
!  1: Pi_RR, 2: Pi_ZZ, 3: Pi_phiphi, 4: Pi_RZ, 5: Pi_Rphi, 6: Pi_Zphi
allocate(jorek_feedback%rhs(n_degrees, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 6))

jorek_feedback%rhs = 0.d0

aux_node_list => jorek_feedback%node_list

if(sim%my_id .eq. 0) write(*,*) "T_EP_eV = ", T_EP_eV
call set_T_EP_eV(T_EP_eV)

! Set up 4D phase space projection: R, Z, P (normalized velocity), Energy
RZ_dist = new_phase_space_projection(ndim=2,res=[91, 91],start=[R_geo-A_MINOR, -A_MINOR], end=[R_geo+A_MINOR, A_MINOR], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_RZ, group=1),basename='RZ_dist')
psiN_dist = new_phase_space_projection(ndim=1,res=[201],start=[0.d0], end=[1.d0], &
                                       f_proj = proj_f(proj_one_weighted,1), f_grids=proj_ndim_f(f=proj_psiN, group=1),basename='psiN_dist')
KE_dist = new_phase_space_projection(ndim=1,res=[151],start=[0.d0], end=[1.d0], &
                                       f_proj = proj_f(proj_Ekin,1), f_grids=proj_ndim_f(f=proj_psiN, group=1),basename='KE_dist')
! fourD_dist = new_phase_space_projection(ndim=4,res=[31,31,41,81],start=[9.d0,-1.d0,-1.1d0, -20.d0],end=[11.0,1.d0,1.1d0,1700.d0], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_RZPE, group=1),basename='fourD_dist')
! fourD_dist = new_phase_space_projection(ndim=2,res=[51,81],start=[-1.1d0, -20.d0],end=[1.1d0,1700.d0], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_PparE, group=1),basename='PparE_dist')
!for electrons
PparMu_dist = new_phase_space_projection(ndim=2,res=[201,201],start=[-2.5, 0d0],end=[2.5, 1.d3], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_PparMu, group=1), basename='PparMu_dist') ! project particles on the grid (1) ppar/mc (2) mu [keV/T]

!for IONS
! PparMu_dist = new_phase_space_projection(ndim=2,res=[301,301],start=[-0.045, 0d0],end=[0.045, 2.5d2], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_PparMu, group=1), basename='PparMu_dist') ! project particles on the grid (1) ppar/mc (2) mu [keV/T]
call with(sim,PparMu_dist)
call output_phase_project(PparMu_dist,0,output_grids_in=.true.)
call with(sim,psiN_dist)
call output_phase_project(psiN_dist,0,output_grids_in=.true.)
call with(sim,KE_dist)
call output_phase_project(KE_dist,0,output_grids_in=.true.)
! call with(sim,RZ_dist)
! call output_phase_project(RZ_dist,0,output_grids_in=.true.)

! Set up (E, P_phi, mu) projection to directly verify invariance of f(E,P_phi,mu)
! dim 1: kinetic energy [keV], dim 2: canonical toroidal angular momentum P_phi/amu [AMU·m²/s], dim 3: magnetic moment [keV/T]
! EPphiMu_dist = new_phase_space_projection(ndim=3,res=[81,281,75], &
!                start=[0.d0, 5.d5, 0.d0], end=[8.d2, 1.0d7, 8.d2], &
!                f_proj = proj_f(proj_one,1), &
!                f_grids=proj_ndim_f(f=proj_EPphiMu, group=1), &
!                basename='EPphiMu_dist')
! call with(sim,EPphiMu_dist)
! call output_phase_project(EPphiMu_dist,0,output_grids_in=.true.)

! E in [keV], mu in [keV/T]]
EMu_dist = new_phase_space_projection(ndim=2,res=[151, 151], &
               start=[0.d0, 0.d0], end=[4.5d2, 4.5d2], &
               f_proj = proj_f(proj_one_weighted,1), &
               f_grids=proj_ndim_f(f=proj_EMu, group=1), &
               basename='EMu_dist')
call with(sim,EMu_dist)
call output_phase_project(EMu_dist,0,output_grids_in=.true.)

! The output is only done on the root process. To prevent a too big imbalance, barrier here.
call MPI_BARRIER(MPI_COMM_WORLD,ifail)
 

! Set up 3D power exchange. Each family is independent and enabled by its own flag:
!   output_pe_E_mu     -> (psi_N, E, mu)     E [keV], mu [eV/T]
!   output_pe_mu_Pphi  -> (psi_N, mu, P_phi) mu [eV/T], P_phi [AMU·m²/s]
!   use_pe_EmuPphi     -> (E, mu, P_phi)     E [keV], mu [eV/T], P_phi [AMU·m²/s]
! The families are NOT exclusive: any combination may be enabled. Only enabled
! families are allocated/output, to save memory and compute time.

! (psi_N, E, mu): dim 1 psi_N, dim 2 E=(gamma-1) m c^2 [keV], dim 3 mu [eV/T]
if (output_pe_E_mu) then
  power_exchange_Emu = new_phase_space_projection(ndim=3,res=[91,201,181], start=[0.d0,0.d0,0.d0], end=[1.0,1.d3,1.0d6],basename="power_exchange_psiN_EMu")
power_exchange_Emu%values = 0.0d0
  power_exchange_Emu_trap = new_phase_space_projection(ndim=3,res=[91,201,181], start=[0.d0,0.d0,0.d0], end=[1.0,1.d3,1.0d6],basename="power_exchange_psiN_EMu_trap")
power_exchange_Emu_trap%values = 0.0d0
  power_exchange_Emu_pass = new_phase_space_projection(ndim=3,res=[91,201,181], start=[0.d0,0.d0,0.d0], end=[1.0,1.d3,1.0d6],basename="power_exchange_psiN_EMu_pass")
power_exchange_Emu_pass%values = 0.0d0
  call output_phase_project(power_exchange_Emu, 0, output_grids_in=.true.)
end if

! (psi_N, mu, P_phi): dim 1 psi_N, dim 2 mu [eV/T], dim 3 P_phi [AMU·m²/s]
!   P_phi = q*psi + p_par_SI * R * Bphi/B   (SI, then /ATOMIC_MASS_UNIT)
!   Bphi = B(3) is the physical toroidal component from calc_EBpsiU
if (output_pe_mu_Pphi) then
  power_exchange_vpar_mu = new_phase_space_projection(ndim=3,res=[91,281,301], start=[0.d0,0.d0,7.d5], end=[1.0,1.0d6,0.8d7],basename="power_exchange_psiN_muPphi")
  power_exchange_vpar_mu%values = 0.0d0
  power_exchange_vpar_mu_trap = new_phase_space_projection(ndim=3,res=[91,281,301], start=[0.d0,0.d0,7.d5], end=[1.0,1.0d6,0.8d7],basename="power_exchange_psiN_muPphi_trap")
  power_exchange_vpar_mu_trap%values = 0.0d0
  power_exchange_vpar_mu_pass = new_phase_space_projection(ndim=3,res=[91,281,301], start=[0.d0,0.d0,7.d5], end=[1.0,1.0d6,0.8d7],basename="power_exchange_psiN_muPphi_pass")
  power_exchange_vpar_mu_pass%values = 0.0d0
  call output_phase_project(power_exchange_vpar_mu, 0, output_grids_in=.true.)
end if

! (E, mu, P_phi): the 3 invariants of motion.
!   dim 1: E     (kinetic energy, (gamma-1) m c^2)      [keV]
!   dim 2: mu    (magnetic moment)                      [eV/T]
!   dim 3: P_phi (canonical toroidal angular momentum)  [AMU·m²/s]
! Grid: moderate memory, high resolution on P_phi.
if (use_pe_EmuPphi) then
  power_exchange_EmuPphi = new_phase_space_projection(ndim=3,res=[101,121,401], start=[0.d0,0.d0,7.d5], end=[1.d3,1.0d6,0.8d7],basename="power_exchange_EmuPphi")
  power_exchange_EmuPphi%values = 0.0d0
  power_exchange_EmuPphi_trap = new_phase_space_projection(ndim=3,res=[101,121,401], start=[0.d0,0.d0,7.d5], end=[1.d3,1.0d6,0.8d7],basename="power_exchange_EmuPphi_trap")
  power_exchange_EmuPphi_trap%values = 0.0d0
  power_exchange_EmuPphi_pass = new_phase_space_projection(ndim=3,res=[101,121,401], start=[0.d0,0.d0,7.d5], end=[1.d3,1.0d6,0.8d7],basename="power_exchange_EmuPphi_pass")
  power_exchange_EmuPphi_pass%values = 0.0d0
  call output_phase_project(power_exchange_EmuPphi, 0, output_grids_in=.true.)
end if

! Project particles on the grid (1) r (2) L = (omega_0 - n omega_phi) / omega_theta
! res_num_trap = new_phase_space_projection(ndim=2,res=[91,151],start=[0.d0, -13.d0],end=[amin, +13.d0], basename='res_num_trap_r')  ! for r
res_num_trap = new_phase_space_projection(ndim=2,res=[81,161],start=[0.d0, -13.d0],end=[1.d0, +13.d0], basename='res_num_trap_psiN') 
res_num_trap%values = 0.d0
call output_phase_project(res_num_trap, 0, output_grids_in=.true.)
! res_num_pass = new_phase_space_projection(ndim=2,res=[91,151],start=[0.d0, -13.d0],end=[amin, +13.d0], basename='res_num_pass_r') ! for r
res_num_pass = new_phase_space_projection(ndim=2,res=[81,161],start=[0.d0, -13.d0],end=[1.d0, +13.d0], basename='res_num_pass_psiN') 
res_num_pass%values = 0.d0
call output_phase_project(res_num_pass, 0, output_grids_in=.true.)

! For proper timestepping, the projections need to be defined before the jorek timestepper
jorek_stepper = new_jorek_timestep_action(jorek_feedback%node_list)

diag = write_particle_diagnostics(filename='diag.h5', append=.true.)

if (restart) then
   tstart_jorek = sim%time + tstep_si
else
   tstart_jorek = sim%time
endif

if (sim%my_id .eq. 0) write(*,*) 'tstart_jorek : ',tstart_jorek


diag_time = timesteps
events = [new_event_ptr(jorek_feedback,  start = tstart_jorek),  &
          new_event_ptr(jorek_stepper,   start = tstart_jorek), &
!          new_event_ptr(project_density, start = tstart_jorek, step=diag_time), &
!          new_event_ptr(project_current, start = tstart_jorek, step=diag_time), &
!          event(write_action(), step=diag_time), &
!          event(diag, step=diag_time) &
          event(stop_action(), start=1d12)  ]

jorek_stepper%extra_event => events(1)


! Set up random numbers for ionisation probability
seed = random_seed()
n_stream = 1
!$ n_stream = omp_get_max_threads()
allocate(rng(n_stream))
do i=1,n_stream
  call rng(i)%initialize(1, seed, n_stream, i)
end do

! Call events at sim%time once to help event scheduler, before entering particle loop
step_rest_time = 0.d0
ino = 0
if (.not. restart_particles) then
  call MPI_BARRIER(MPI_COMM_WORLD,ifail)
  call with(sim, events, at=sim%time)
endif
! call with(sim, events, at=sim%time) ! need to call this before the particle loop to properly initialize the jorek_stepper's internal time and timestep.

!call with(sim, project_density)

do while (.not. sim%stop_now)
  if(sim%my_id .eq. 0) write(*,*) "fields%static = " , sim%fields%static

  target_time = next_event_at(sim, events) 
  particle_start_time = (sim%time - step_rest_time)
  particle_step_time  = target_time - particle_start_time
  n_steps             = particle_step_time/timesteps
  step_rest_time      = particle_step_time - real(n_steps,8) * timesteps

  if (sim%my_id .eq. 0) then
     if (n_steps < 10) write(*,*) 'low n_steps,', n_steps
     write(*,*) 'Time difference between particles and jorek: ', step_rest_time
     write(*,*) "PARTICLE : target time         : ",target_time
     write(*,*) "PARTICLE : timesteps           : ",timesteps
     write(*,*) "PARTICLE : sim%time            : ",sim%time
     write(*,*) "PARTICLE : particle_start_time : ",particle_start_time
     write(*,*) "PARTICLE : particle_step_time  : ",particle_step_time
     write(*,*) "PARTICLE : n_steps             : ",n_steps
     write(*,*) "PARTICLE : step_rest_time      : ",step_rest_time
  endif
  
  t_loop_start = MPI_WTIME()
  select case(trim(particle_pusher))
  case('gc_rel')
    call loop_particle_gc_relativistic(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time, power_exchange_vpar_mu, power_exchange_vpar_mu_trap, power_exchange_vpar_mu_pass, power_exchange_Emu, power_exchange_Emu_trap, power_exchange_Emu_pass, power_exchange_EmuPphi, power_exchange_EmuPphi_trap, power_exchange_EmuPphi_pass, &
                                       theta_prev_arr, theta_accum_arr, dtheta_prev_arr, trap_pass_arr, &
                                       phi_prev_arr, phi_unwrap_arr, &
                                       t_ref_a_arr, phi_ref_a_arr, t_ref_b_arr, phi_ref_b_arr, &
                                       n_periods_arr, omega_theta_buf_arr, omega_phi_buf_arr)
  case('kinetic_rel')
    call loop_particle_kinetic_relativistic(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time, power_exchange_vpar_mu, power_exchange_Emu, power_exchange_EmuPphi, &
                                       n_periods_arr, omega_theta_buf_arr, omega_phi_buf_arr)
  end select
  t_loop_end = MPI_WTIME()
  t_loop_time = t_loop_end - t_loop_start

  ! Per-rank timing: each rank reports immediately (before any barrier) so the
  ! slowest rank is visible even if others are already waiting at MPI_BARRIER.
  write(*,'(A,I5,A,I5,A,F12.3,A)') &
    "RANK ", sim%my_id, " step ", ino + index_start, " particle loop: ", t_loop_time, " s"
  call MPI_Reduce(t_loop_time, dummy_real8_1, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ifail)
  call MPI_Reduce(t_loop_time, dummy_real8_2, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ifail)
  if(sim%my_id .eq. 0) write(*,'(A,I5,A,F12.3,A,F12.3,A)') &
    "SUMMARY step ", ino + index_start, &
    " particle loop: max=", dummy_real8_1, " s  min=", dummy_real8_2, " s"

  ino = ino + 1
  ! Output power exchange every POWER_EXCHANGE_STEPS timesteps
  if(mod(ino,POWER_EXCHANGE_STEPS) .eq. 0) then
    ! (psi_N, E, mu)
    if (output_pe_E_mu) then
    if(use_trap_passing_for_PE .and. trim(particle_pusher) == 'gc_rel') then
        call output_phase_project(power_exchange_Emu_trap, ino + index_start, output_grids_in=.false.)
        call output_phase_project(power_exchange_Emu_pass, ino + index_start, output_grids_in=.false.)
      else
        call output_phase_project(power_exchange_Emu, ino + index_start, output_grids_in=.false.)
      end if
      power_exchange_Emu%values = 0.d0
      power_exchange_Emu_trap%values = 0.d0
      power_exchange_Emu_pass%values = 0.d0
    end if
    ! (psi_N, mu, P_phi)
    if (output_pe_mu_Pphi) then
      if(use_trap_passing_for_PE .and. trim(particle_pusher) == 'gc_rel') then
        call output_phase_project(power_exchange_vpar_mu_trap, ino + index_start, output_grids_in=.false.)
        call output_phase_project(power_exchange_vpar_mu_pass, ino + index_start, output_grids_in=.false.)
      else
        call output_phase_project(power_exchange_vpar_mu, ino + index_start, output_grids_in=.false.)
    end if
    power_exchange_vpar_mu%values = 0.d0
    power_exchange_vpar_mu_trap%values = 0.d0
    power_exchange_vpar_mu_pass%values = 0.d0
    end if
    ! (E, mu, P_phi)
    if (use_pe_EmuPphi) then
      if(use_trap_passing_for_PE .and. trim(particle_pusher) == 'gc_rel') then
        call output_phase_project(power_exchange_EmuPphi_trap, ino + index_start, output_grids_in=.false.)
        call output_phase_project(power_exchange_EmuPphi_pass, ino + index_start, output_grids_in=.false.)
      else
        call output_phase_project(power_exchange_EmuPphi, ino + index_start, output_grids_in=.false.)
      end if
      power_exchange_EmuPphi%values = 0.d0
      power_exchange_EmuPphi_trap%values = 0.d0
      power_exchange_EmuPphi_pass%values = 0.d0
    end if
  endif

  ! Plot r v.s. L = (omega_0 - n omega_phi) / omega_theta for resonant particles every POWER_EXCHANGE_STEPS timesteps
  ! Only applicable for GC particles (trap/pass classification)
  ! if(mod(ino,POWER_EXCHANGE_STEPS*2) .eq. 0 .and. trim(particle_pusher) == 'gc_rel') then
  !   if(use_trap_passing_for_PE) then
  !     call output_phase_project(res_num_trap, ino + index_start, output_grids_in=.false.)
  !     call output_phase_project(res_num_pass, ino + index_start, output_grids_in=.false.)
  !   end if
  !   res_num_trap%values = 0.d0
  !   res_num_pass%values = 0.d0
  ! endif
  
  ! Output 4D distribution function periodically
  if(mod(ino,RZ_OUTPUT_STEPS) .eq. 0) then 
    ! call with(sim,RZ_dist)
    ! call output_phase_project(RZ_dist,ino + index_start,output_grids_in=.false.)
    call with(sim,psiN_dist)
    call output_phase_project(psiN_dist,ino + index_start,output_grids_in=.false.)
    ! call with(sim,KE_dist)
    ! call output_phase_project(KE_dist,ino + index_start,output_grids_in=.false.)
    ! call with(sim,EPphi_dist)
    ! call output_phase_project(EPphi_dist,ino + index_start,output_grids_in=.false.)
    ! call with(sim,EPphiMu_dist)
    ! call output_phase_project(EPphiMu_dist,ino + index_start,output_grids_in=.false.)
    call MPI_BARRIER(MPI_COMM_WORLD,ifail)
  endif

  sim%time = target_time 
  
  t_barrier_start = MPI_WTIME()
  call MPI_BARRIER(MPI_COMM_WORLD,ifail)
  if(sim%my_id .eq. 0) write(*,'(A,F10.3,A)') "Barrier wait time: ", MPI_WTIME()-t_barrier_start, " seconds"
  
  call with(sim, events, at=sim%time)

  ! Save particles every PARTICLE_OUTPUT_STEPS timesteps AFTER JOREK timestepping
  ! This ensures particles are consistent with the JOREK fields at this timestep
  if (mod(ino, PARTICLE_OUTPUT_STEPS) .eq. 0) then
    write(part_output_filename, '(A,I6.6,A)') 'part_', ino + index_start, '.h5'
    if (sim%my_id .eq. 0) write(*,*) 'Saving particles to: ', trim(part_output_filename)
    partwriter = event(write_action(filename=trim(part_output_filename)))
    call with(sim, partwriter)
  endif

end do

!call write_simulation_hdf5(sim, 'part_restart.h5')

! partwriter = event(write_action(filename='part_restart.h5'))
! call with(sim, partwriter)

call sim%finalize

if (sim%my_id == 0) close(111)

contains

subroutine loop_particle_kinetic_relativistic(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time, power_exchange_vpar_mu, power_exchange_Emu, power_exchange_EmuPphi, &
                                          n_periods_arr, omega_theta_buf_arr, omega_phi_buf_arr)
use mod_project_particles
use mod_random_seed
use mod_interp, only: mode_moivre
use mod_basisfunctions
use mod_particle_types, only: copy_particle
use mod_kinetic_relativistic, only: runge_kutta_fixed_dt_relativistic_particle_push_jorek, &
                                     relativistic_kinetic_to_relativistic_gc
use mod_import_experimental_dist, only: calculate_B 

implicit none

class(particle_sim), target, intent(inout)                :: sim
type(projection), target, intent(inout)                   :: jorek_feedback
type(phase_space_projection), target, intent(inout)       :: power_exchange_vpar_mu, power_exchange_Emu, power_exchange_EmuPphi
type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
type(particle_kinetic_relativistic)                       :: particle_tmp, particle_mid
type(particle_gc_relativistic)                            :: gc_mid

! Dummy arrays for interface compatibility (unused for kinetic — trap/pass disabled)
integer, intent(inout) :: n_periods_arr(:)
real*8,  intent(inout) :: omega_theta_buf_arr(:,:), omega_phi_buf_arr(:,:)

real*8, intent(in)     :: timesteps, particle_start_time 
real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm, tstep_si
real*8    :: v_kin_temp, E(3), B(3), psi, psiN, psiN_mid, U, n_e, T_e, rz_old(2), st_old(2), r, v, tmp2(2) 
real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
real*8    :: gamma, p_mag, gamma_mid, p_mag_mid, target_time, dt_local, tmid, t
real*8    :: fitsinevpar(4), fitsinemu(4), fitsineE(4), fitsineB(4), omega, mumid, vparmid, pparmid, Pphi_mid, E_mid_keV, mu_start, vpar_start, B_start, mu_end, vpar_end, B_end, E_diff, E_i, avgB, B_mag
real*8    :: rcontainer(n_steps)
logical   :: midpoint_stored
real*8, parameter :: AMU_C2_EV = ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 / EL_CHG ! rest-mass energy of 1 amu in eV. 
!$ real*8 :: w0, w1, mmm(3)

integer, intent(in)   :: n_steps
integer   :: i, j, k, l, m, i_elm_old, i_elm 
integer   :: seed, n_stream, ierr, nthreads
integer   :: i_tor, index_lm, i_elm_temp
integer   :: n_particles, ifail
real*8,allocatable :: feedback_rhs(:,:,:,:,:)
real*8, allocatable :: phase_proj_1(:) ! power exchange projection (no trap/pass split for kinetic)
real*8, allocatable :: phase_proj_1_E(:) ! (psi_N, mu, E) power exchange projection (kinetic)
real*8, allocatable :: phase_proj_1_EP(:) ! (E, mu, P_phi) power exchange projection (kinetic)

!$ w0 = omp_get_wtime()

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = central_mass * ATOMIC_MASS_UNIT * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation

jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps
allocate(feedback_rhs,source=jorek_feedback%rhs)
! Allocate with correct shape, but initialise to 0 (not source values).
! Allocate accumulators per enabled family (size-1 dummies keep OMP reductions valid).
if (output_pe_mu_Pphi) then
allocate(phase_proj_1,    mold=power_exchange_vpar_mu%values);      phase_proj_1    = 0.d0
else
  allocate(phase_proj_1(1), source=0.d0)
end if
if (output_pe_E_mu) then
allocate(phase_proj_1_E,  mold=power_exchange_Emu%values);          phase_proj_1_E  = 0.d0
else
  allocate(phase_proj_1_E(1), source=0.d0)
end if
if (use_pe_EmuPphi) then
  allocate(phase_proj_1_EP,  mold=power_exchange_EmuPphi%values);      phase_proj_1_EP  = 0.d0
else
  allocate(phase_proj_1_EP(1), source=0.d0)
end if

jorek_feedback%rhs = 0.d0
feedback_rhs       = 0.d0

! Calculate time values before parallel region
target_time = particle_start_time + n_steps * timesteps
tmid = particle_start_time + 0.5d0 * (n_steps * timesteps)


select type (particles => sim%groups(1)%particles)
type is (particle_kinetic_relativistic)
#ifdef __GFORTRAN__
 !$omp parallel do default(shared) & 
#else
 !$omp parallel do default(none) &
 !$omp shared(sim, particles, n_steps, timesteps, rng, particle_start_time,        &
 !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, target_time, tmid,       &
 !$omp jorek_feedback, power_exchange_vpar_mu, power_exchange_Emu, power_exchange_EmuPphi, &
 !$omp output_pe_E_mu, output_pe_mu_Pphi, use_pe_EmuPphi, ES)                              &
#endif
 !$omp private(particle_tmp, particle_mid, gc_mid, i,j,k,l,m, t, E, B, psi, psiN, psiN_mid, U, rz_old, st_old, r, v, tmp2,    &
 !$omp i_elm_old, i_elm, n_e, T_e, E_diff, avgB, B_mag,                                               & 
 !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, ifail, gamma, p_mag, gamma_mid, p_mag_mid, &
 !$omp omega, vparmid, pparmid, mumid, Pphi_mid, E_mid_keV, mu_start, vpar_start, B_start, mu_end, vpar_end, B_end, &
 !$omp dt_local, midpoint_stored) &
 !$omp schedule(dynamic,10) &
 !$omp reduction(+:feedback_rhs)&
 !$omp reduction(+:phase_proj_1)&
 !$omp reduction(+:phase_proj_1_E)&
 !$omp reduction(+:phase_proj_1_EP)
  do j=1,size(particles,1)
    call copy_particle(particle_tmp, particles(j))

    ! Lost particles (i_elm < 0) and Uninitialized/error particles (i_elm = 0)
    if (particle_tmp%i_elm .le. 0) then
      cycle  ! Skip particles with invalid element indices
    endif

    call sim%fields%calc_EBpsiU(particle_start_time, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
    B_mag = norm2(B)
    ! Initial kinetic energy (for E_diff computation at end)
    p_mag = sqrt(particle_tmp%p(1)**2 + particle_tmp%p(2)**2 + particle_tmp%p(3)**2)
    gamma = sqrt(1.d0 + (p_mag / (sim%groups(1)%mass * SPEED_OF_LIGHT))**2)
    E_diff = (gamma - 1.d0) * sim%groups(1)%mass * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2
    
    ! Initialize time stepping
    t = particle_start_time
    dt_local = timesteps
    midpoint_stored = .false.
    k=1
    
    ! Push the particle with fixed time stepping (adaptive not available for kinetic)
    do while ((t .lt. target_time) .and. (particle_tmp%i_elm .gt. 0))
      ! Store midpoint data for phase space projection
      if ((.not. midpoint_stored) .and. (t .ge. tmid)) then
        call copy_particle(particle_mid, particle_tmp)
        midpoint_stored = .true.
      endif

      ! Push the particle — fixed dt only for kinetic particles
      call runge_kutta_fixed_dt_relativistic_particle_push_jorek(sim%fields, t, dt_local, sim%groups(1)%mass, particle_tmp)
        
      t = t + dt_local
      k = k + 1
    end do ! time steps

    call copy_particle(particles(j), particle_tmp)

    i_elm = particle_tmp%i_elm

    if (i_elm .gt. 0 .and. i_elm .le. sim%fields%element_list%n_elements) then

      call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
      call mode_moivre(particle_tmp%x(3), HZ)
      call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
      B_mag = norm2(B)

      ! precompute the quantities for the feedback projection
      p_mag = sqrt(particle_tmp%p(1)**2 + particle_tmp%p(2)**2 + particle_tmp%p(3)**2)
      gamma = sqrt(1.d0 + (p_mag / (sim%groups(1)%mass * SPEED_OF_LIGHT))**2)

      ! Safety: skip feedback and projections for runaway particles.
      if (gamma >= 1.d2) cycle

      do l=1,n_vertex_max
        do m=1,n_degrees

          index_lm = (l-1)*n_degrees + m

          v = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m) 

          do i_tor=1,n_tor
            ! Isotropic pressure tensor: (1/3) p^2/(gamma m) on the diagonal,
            ! zero off-diagonal (components 4..6 remain 0). This keeps kinetic_rel
            ! correct also when use_pcs_full is enabled.
            feedback_rhs(m,l,i_elm,i_tor,1) = feedback_rhs(m,l,i_elm,i_tor,1) &
                                                  + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                  * (1.d0/3.d0) * (p_mag**2 / gamma / sim%groups(1)%mass) * mu_zero !Pi_RR
            feedback_rhs(m,l,i_elm,i_tor,2) = feedback_rhs(m,l,i_elm,i_tor,2) &
                                                  + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                  * (1.d0/3.d0) * (p_mag**2 / gamma / sim%groups(1)%mass) * mu_zero !Pi_ZZ
            feedback_rhs(m,l,i_elm,i_tor,3) = feedback_rhs(m,l,i_elm,i_tor,3) &
                                                  + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                  * (1.d0/3.d0) * (p_mag**2 / gamma / sim%groups(1)%mass) * mu_zero !Pi_phiphi
          enddo

        enddo   !< order
      enddo     !< vertex

      ! Compute power exchange only when the particle is still in the domain
      psiN = (psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis)
      if(psiN > 0.01 .and. psiN < 0.99 .and. midpoint_stored) then
        ! Final kinetic energy for E_diff
        p_mag = sqrt(particle_tmp%p(1)**2 + particle_tmp%p(2)**2 + particle_tmp%p(3)**2)
        gamma = sqrt(1.d0 + (p_mag / (sim%groups(1)%mass * SPEED_OF_LIGHT))**2)
        E_diff = -E_diff + (gamma - 1.d0) * sim%groups(1)%mass * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2

        ! Convert kinetic particle to GC state at midpoint for diagnostics
        call sim%fields%calc_EBpsiU(tmid, particle_mid%i_elm, particle_mid%st, particle_mid%x(3), E, B, psi, U)
        B_mag = norm2(B)
        
        gc_mid = relativistic_kinetic_to_relativistic_gc(sim%fields%node_list, &
                  sim%fields%element_list, particle_mid, sim%groups(1)%mass, B)
        p_mag_mid = sqrt(particle_mid%p(1)**2 + particle_mid%p(2)**2 + particle_mid%p(3)**2)
        gamma_mid = sqrt(1.d0 + (p_mag_mid / (sim%groups(1)%mass * SPEED_OF_LIGHT))**2)
        
        psiN_mid = (psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis)
        mumid   = gc_mid%p(2) * ATOMIC_MASS_UNIT / EL_CHG  ! in [eV/Tesla]
        vparmid = gc_mid%p(1) / (gamma_mid * sim%groups(1)%mass) ! in [m/s]
        ! P_phi (canonical toroidal angular momentum) at midpoint, in [AMU·m²/s]
        !   P_phi = q*psi + p_par_SI * R * Bphi/B   (SI), Bphi = B(3) physical toroidal component
        Pphi_mid = ( gc_mid%x(1) * gc_mid%p(1) * ATOMIC_MASS_UNIT * (B(3)/B_mag) &
                     + dble(gc_mid%q) * EL_CHG * psi ) / ATOMIC_MASS_UNIT
        ! E = (gamma-1) m c^2 at midpoint, in [keV]
        E_mid_keV = (gamma_mid - 1.d0) * sim%groups(1)%mass * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 / EL_CHG / 1.d3
        
        if (sim%my_id .eq. 0) then
          if(mod(j,500) .eq. 0) then 
            write(*,*) 'DEBUG: kinetic particle ', j, ' at psiN_mid= ', psiN_mid, ' mumid [eV/T] = ', mumid, ' Pphi_mid = ', Pphi_mid, ' E_mid_keV = ', E_mid_keV, ' E_diff [J] = ', E_diff
          endif
        endif
        
        ! No trap/pass split for kinetic — project to combined power exchange only
        if (output_pe_E_mu)     call project_single_particle_x(power_exchange_Emu,[psiN_mid,E_mid_keV,mumid],phase_proj_1_E,E_diff)
        if (output_pe_mu_Pphi)  call project_single_particle_x(power_exchange_vpar_mu,[psiN_mid,mumid,Pphi_mid],phase_proj_1,E_diff)
        if (use_pe_EmuPphi)     call project_single_particle_x(power_exchange_EmuPphi,[E_mid_keV,mumid,Pphi_mid],phase_proj_1_EP,E_diff)
  
      endif ! psiN and midpoint check for projection 
    endif ! i_elm > 0

  end do   ! particles
  !$omp end parallel do
  
  if (sim%my_id .eq. 0) write(*,*) "End of the kinetic particle loop"

end select

jorek_feedback%rhs = feedback_rhs
! No trap/pass split for kinetic — only combined power exchange
if (output_pe_E_mu)     power_exchange_Emu%values    = power_exchange_Emu%values    + phase_proj_1_E
if (output_pe_mu_Pphi)  power_exchange_vpar_mu%values = power_exchange_vpar_mu%values + phase_proj_1
if (use_pe_EmuPphi)     power_exchange_EmuPphi%values = power_exchange_EmuPphi%values + phase_proj_1_EP
deallocate(phase_proj_1)
deallocate(phase_proj_1_E)
deallocate(phase_proj_1_EP)
deallocate(feedback_rhs)

end subroutine loop_particle_kinetic_relativistic


subroutine loop_particle_gc_relativistic(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time, power_exchange_vpar_mu, power_exchange_vpar_mu_trap, power_exchange_vpar_mu_pass, power_exchange_Emu, power_exchange_Emu_trap, power_exchange_Emu_pass, power_exchange_EmuPphi, power_exchange_EmuPphi_trap, power_exchange_EmuPphi_pass, &
                                          theta_prev_arr, theta_accum_arr, dtheta_prev_arr, trap_pass_arr, &
                                          phi_prev_arr, phi_unwrap_arr, &
                                          t_ref_a_arr, phi_ref_a_arr, t_ref_b_arr, phi_ref_b_arr, &
                                          n_periods_arr, omega_theta_buf_arr, omega_phi_buf_arr)
use mod_project_particles
use mod_random_seed
use mod_interp, only: mode_moivre
use mod_basisfunctions
use mod_particle_types, only: copy_particle
use mod_gc_relativistic, only: compute_relativistic_factor, rel_kinetic_energy, runge_kutta_fixed_dt_gc_push_jorek
use mod_import_experimental_dist, only: calculate_B 

implicit none

class(particle_sim), target, intent(inout)                :: sim
type(projection), target, intent(inout)                   :: jorek_feedback
type(phase_space_projection), target, intent(inout)       :: power_exchange_vpar_mu, power_exchange_vpar_mu_trap, power_exchange_vpar_mu_pass, power_exchange_Emu, power_exchange_Emu_trap, power_exchange_Emu_pass, power_exchange_EmuPphi, power_exchange_EmuPphi_trap, power_exchange_EmuPphi_pass
type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
type(particle_gc_relativistic)                            :: particle_tmp, particle_mid

! Persistent classification state (one entry per local particle, survives across calls)
! theta_prev_arr  : poloidal angle at end of previous call, in [0, 2pi)
! theta_accum_arr : signed net angle accumulated since the last completed orbit
! dtheta_prev_arr : last angular step (0 at first call => reversal check naturally skipped)
! trap_pass_arr   : 0 = undetermined, 1 = trapped, -1 = passing
real*8,  intent(inout) :: theta_prev_arr(:)
real*8,  intent(inout) :: theta_accum_arr(:)
real*8,  intent(inout) :: dtheta_prev_arr(:)
integer, intent(inout) :: trap_pass_arr(:)
real*8,  intent(inout) :: phi_prev_arr(:), phi_unwrap_arr(:)
real*8,  intent(inout) :: t_ref_a_arr(:), phi_ref_a_arr(:)
real*8,  intent(inout) :: t_ref_b_arr(:), phi_ref_b_arr(:)
integer, intent(inout) :: n_periods_arr(:)
real*8,  intent(inout) :: omega_theta_buf_arr(:,:), omega_phi_buf_arr(:,:)

real*8, intent(in)     :: timesteps, particle_start_time 
real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm, tstep_si
real*8    :: v_kin_temp, E(3), B(3), psi, psiN, psiN_mid, U, n_e, T_e, rz_old(2), st_old(2), r, v, tmp2(2) 
real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
real*8    :: gamma, target_time, dt_local, tmid, t
real*8    :: b_norm_r, b_norm_z, b_norm_phi, p2_par, p2_perp, p_atrop
logical   :: fb_valid
real*8    :: fitsinevpar(4), fitsinemu(4), fitsineE(4), fitsineB(4), omega, mumid, vparmid, pparmid, Pphi_mid, E_mid_keV, mu_start, vpar_start, B_start, mu_end, vpar_end, B_end, E_diff, E_i, avgB, B_mag
real*8    :: rcontainer(n_steps)
logical   :: midpoint_stored
real*8, parameter :: AMU_C2_EV = ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 / EL_CHG ! rest-mass energy of 1 amu in eV. 
!$ real*8 :: w0, w1, mmm(3)

! Per-particle (private) poloidal-angle / frequency tracking variables
real*8  :: theta_curr_p, theta_prev_p, theta_accum_p, dtheta_p, dtheta_prev_p
real*8  :: phi_unwrap_p, phi_prev_p
real*8  :: t_ref_a_p, phi_ref_a_p, t_ref_b_p, phi_ref_b_p
integer :: n_periods_p
real*8  :: omega_theta_buf_p(MIN_PERIODS_RESONANCE), omega_phi_buf_p(MIN_PERIODS_RESONANCE)
real*8  :: omega_theta_avg_p, omega_phi_avg_p, L_res
integer :: n_use_p

integer, intent(in)   :: n_steps
integer   :: i, j, k, l, m, i_elm_old, i_elm 
integer   :: seed, n_stream, ierr, nthreads
integer   :: i_tor, index_lm, i_elm_temp
integer   :: n_particles, ifail
integer   :: log_unit
character(len=50) :: log_filename
real*8,allocatable :: feedback_rhs(:,:,:,:,:)
real*8, allocatable :: phase_proj_1(:), phase_proj_trap(:) ! phase_proj_1 is the projection used for all particles or passing particles (depending on use_trap_passing_for_PE)
real*8, allocatable :: phase_proj_1_E(:), phase_proj_trap_E(:) ! (psi_N, mu, E) power exchange accumulators
real*8, allocatable :: phase_proj_1_EP(:), phase_proj_trap_EP(:) ! (E, mu, P_phi) power exchange accumulators
real*8, allocatable :: E_diff_arr(:) ! per-particle energy difference for res_num projection

!$ w0 = omp_get_wtime()

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = central_mass * ATOMIC_MASS_UNIT * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation

jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps
allocate(feedback_rhs,source=jorek_feedback%rhs)
! Allocate with correct shape, but initialise to 0 (not source values).
! Using source= here would seed the OMP reduction with the old accumulated
! values, causing those values to be counted twice when added back below.
! Allocate accumulators per enabled family. Disabled families get size-1 dummies so
! the OMP reduction clauses stay valid even though they are never written.
! (psi_N, mu, P_phi)
if (output_pe_mu_Pphi) then
allocate(phase_proj_1,    mold=power_exchange_vpar_mu%values);      phase_proj_1    = 0.d0
allocate(phase_proj_trap, mold=power_exchange_vpar_mu_trap%values);  phase_proj_trap = 0.d0
else
  allocate(phase_proj_1(1), source=0.d0)
  allocate(phase_proj_trap(1), source=0.d0)
end if
! (psi_N, E, mu)
if (output_pe_E_mu) then
allocate(phase_proj_1_E,    mold=power_exchange_Emu%values);        phase_proj_1_E  = 0.d0
allocate(phase_proj_trap_E, mold=power_exchange_Emu_trap%values);    phase_proj_trap_E = 0.d0
else
  allocate(phase_proj_1_E(1), source=0.d0)
  allocate(phase_proj_trap_E(1), source=0.d0)
end if
! (E, mu, P_phi)
if (use_pe_EmuPphi) then
  allocate(phase_proj_1_EP,    mold=power_exchange_EmuPphi%values);      phase_proj_1_EP    = 0.d0
  allocate(phase_proj_trap_EP, mold=power_exchange_EmuPphi_trap%values);  phase_proj_trap_EP = 0.d0
else
  allocate(phase_proj_1_EP(1), source=0.d0)
  allocate(phase_proj_trap_EP(1), source=0.d0)
end if

allocate(E_diff_arr(size(sim%groups(1)%particles, 1)));  E_diff_arr = 0.d0

jorek_feedback%rhs = 0.d0
feedback_rhs       = 0.d0

! Calculate time values before parallel region
target_time = particle_start_time + n_steps * timesteps
tmid = particle_start_time + 0.5d0 * (n_steps * timesteps)

! Open per-rank trap/pass log file (append so successive calls accumulate)
! write(log_filename, '(A,I4.4,A)') 'part_trap_pass_', sim%my_id, '.log'
! open(newunit=log_unit, file=trim(log_filename), status='unknown', position='append')

select type (particles => sim%groups(1)%particles)
type is (particle_gc_relativistic)
#ifdef __GFORTRAN__
 !$omp parallel do default(shared) & 
#else
 !$omp parallel do default(none) &
 !$omp shared(sim, particles, n_steps, timesteps, rng, particle_start_time,        &
 !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, target_time, tmid,       &
 !$omp jorek_feedback, power_exchange_vpar_mu, power_exchange_vpar_mu_trap, power_exchange_vpar_mu_pass,                &
 !$omp power_exchange_Emu, power_exchange_Emu_trap, power_exchange_Emu_pass,        &
 !$omp power_exchange_EmuPphi, power_exchange_EmuPphi_trap, power_exchange_EmuPphi_pass, &
 !$omp output_pe_E_mu, output_pe_mu_Pphi, use_pe_EmuPphi,                           &
 !$omp theta_prev_arr, theta_accum_arr, dtheta_prev_arr, trap_pass_arr,            &
 !$omp phi_unwrap_arr, phi_prev_arr,                                               &
 !$omp t_ref_a_arr, phi_ref_a_arr, t_ref_b_arr, phi_ref_b_arr,                     &
 !$omp n_periods_arr, omega_theta_buf_arr, omega_phi_buf_arr, ES,                    &
 !$omp log_unit, E_diff_arr)                              &
#endif
 !$omp private(particle_tmp, particle_mid, i,j,k,l,m, t, E, B, psi, psiN, psiN_mid, U, rz_old, st_old, r, v, tmp2,    &
 !$omp i_elm_old, i_elm, n_e, T_e, E_diff, avgB, B_mag,                                               & 
 !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, ifail, gamma,                &
 !$omp omega, vparmid, pparmid, mumid, Pphi_mid, E_mid_keV, mu_start, vpar_start, B_start, mu_end, vpar_end, B_end, &
 !$omp dt_local, midpoint_stored, fb_valid,                                                            &
 !$omp b_norm_r, b_norm_z, b_norm_phi, p2_par, p2_perp, p_atrop,                                      &
 !$omp theta_curr_p, theta_prev_p, theta_accum_p, dtheta_p, dtheta_prev_p,                              &
 !$omp phi_unwrap_p, phi_prev_p, t_ref_a_p, phi_ref_a_p, t_ref_b_p, phi_ref_b_p,                       &
 !$omp n_periods_p, omega_theta_buf_p, omega_phi_buf_p, n_use_p) &
 !$omp schedule(dynamic,10) &
 !$omp reduction(+:feedback_rhs)&
 !$omp reduction(+:phase_proj_1)&
 !$omp reduction(+:phase_proj_trap)&
 !$omp reduction(+:phase_proj_1_E)&
 !$omp reduction(+:phase_proj_trap_E)&
 !$omp reduction(+:phase_proj_1_EP)&
 !$omp reduction(+:phase_proj_trap_EP)
  do j=1,size(particles,1)
    call copy_particle(particle_tmp, particles(j))

    ! Lost particles (i_elm < 0) and Uninitialized/error particles (i_elm = 0)
    if (particle_tmp%i_elm .le. 0) then
      cycle  ! Skip particles with invalid element indices
    endif

    call sim%fields%calc_EBpsiU(particle_start_time, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
    B_mag = norm2(B)
    E_diff = rel_kinetic_energy(particle_tmp, sim%groups(1)%mass, B_mag)

    ! ---- Initialise per-particle tracking from persistent state ----
    theta_prev_p      = theta_prev_arr(j)
    theta_accum_p     = theta_accum_arr(j)
    dtheta_prev_p     = dtheta_prev_arr(j)
    phi_prev_p        = phi_prev_arr(j)
    phi_unwrap_p      = phi_unwrap_arr(j)
    t_ref_a_p         = t_ref_a_arr(j)
    phi_ref_a_p       = phi_ref_a_arr(j)
    t_ref_b_p         = t_ref_b_arr(j)
    phi_ref_b_p       = phi_ref_b_arr(j)
    n_periods_p       = n_periods_arr(j)
    omega_theta_buf_p = omega_theta_buf_arr(:,j)
    omega_phi_buf_p   = omega_phi_buf_arr(:,j)
    
    ! Initialize time stepping
    t = particle_start_time
    dt_local = timesteps
    midpoint_stored = .false.
    ! rcontainer = 0.d0
    
    ! Push the particle over n_steps substeps. Bounded loop (as in tae_loop_CGL):
    ! executes exactly n_steps iterations for confined particles, which matches
    ! the /n_steps averaging below, and avoids floating-point drift issues of
    ! accumulating t in a do-while condition.
    do k=1,n_steps
      if (particle_tmp%i_elm .le. 0) exit

      ! Time at the start of this substep: computed from k to avoid accumulating
      ! rounding error in t across many substeps.
      t = particle_start_time + real(k-1,8) * dt_local

      ! if(k <= n_steps) rcontainer(k) = sqrt((particle_tmp%x(1)-ES%R_axis)**2 + (particle_tmp%x(2) - ES%Z_axis)**2)
      ! Store midpoint data for phase space projection
      if ((.not. midpoint_stored) .and. (t .ge. tmid)) then
        call copy_particle(particle_mid, particle_tmp)
        midpoint_stored = .true.
      endif

      ! ---- Evaluate fields at the substep start for the CGL feedback ----
      ! (particle_tmp%i_elm > 0 is guaranteed by the loop condition)
      call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
      B_mag = norm2(B)
      gamma = compute_relativistic_factor(particle_tmp, sim%groups(1)%mass, B_mag)
      ! Safety: do not accumulate feedback for runaway particles.
      ! Normal EP have gamma ~ 1.00016 (300 keV D); gamma >= 100 indicates an
      ! unphysical particle that would destabilise GMRES and corrupt diagnostics.
      fb_valid = (gamma .lt. 1.d2)
      ! Normalised b vector (physical R,Z,phi components)
      b_norm_r   = B(1)/B_mag
      b_norm_z   = B(2)/B_mag
      b_norm_phi = B(3)/B_mag
      ! Parallel & perpendicular momenta squared (per AMU)
      p2_par  = particle_tmp%p(1)**2
      p2_perp = 2.d0 * sim%groups(1)%mass * B_mag * particle_tmp%p(2)
      ! Coefficient of b_i b_j in the CGL tensor (divided by gamma*mass below)
      p_atrop = p2_par - 0.5d0 * p2_perp

      ! Push the particle and determine it's new location (fixed time step).
      call runge_kutta_fixed_dt_gc_push_jorek(sim%fields, t, dt_local, sim%groups(1)%mass, particle_tmp)

      ! ---- Trapped/passing classification + omega_theta/omega_phi ----
      call classify_trap_pass(particle_tmp%x(1), particle_tmp%x(2), particle_tmp%x(3), &
                              particle_tmp%i_elm, t, &
                              trap_pass_arr(j), theta_prev_p, theta_accum_p, dtheta_prev_p, &
                              phi_unwrap_p, phi_prev_p, &
                              t_ref_a_p, phi_ref_a_p, t_ref_b_p, phi_ref_b_p, &
                              n_periods_p, omega_theta_buf_p, omega_phi_buf_p)

      ! DEBUG:  Trapped/passing diagnostics (written at every substep) ----
      ! if(mod(j,200) .eq. 0) then
      !   !$omp critical (trap_pass_log)
      !   call write_trap_pass_diag(log_unit, j, particle_tmp%x(1), particle_tmp%x(2), trap_pass_arr(j), t)
      !   !$omp end critical (trap_pass_log)
      ! endif

      ! ---- Pressure-tensor feedback, accumulated per substep ----
      ! B, p2_par and p2_perp are evaluated at the substep start (pre-push), while
      ! the projection location (HH, HZ) is the post-push position - the same
      ! pattern as tae_loop_CGL. Components are ordered as model307 use_pcs_full
      ! expects them: 1: Pi_RR, 2: Pi_ZZ, 3: Pi_phiphi, 4: Pi_RZ, 5: Pi_Rphi, 6: Pi_Zphi.
      !
      ! use_CGL_pressure = .true. : gyrotropic CGL tensor
      !     Pi_ij = [0.5*p_perp^2*delta_ij + (p_par^2 - 0.5*p_perp^2)*b_i b_j] / (gamma*m)
      !   (requires use_pcs_full = .t. in the namelist)
      ! use_CGL_pressure = .false.: isotropic tensor
      !     Pi_ij = (1/3) p^2/(gamma*m) * delta_ij   (off-diagonal components stay 0)
      if (fb_valid .and. particle_tmp%i_elm .gt. 0 .and. particle_tmp%i_elm .le. sim%fields%element_list%n_elements) then
        call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
        call mode_moivre(particle_tmp%x(3), HZ)

        do l=1,n_vertex_max
          do m=1,n_degrees

            v = HH(l,m) * sim%fields%element_list%element(particle_tmp%i_elm)%size(l,m)

            do i_tor=1,n_tor

              if (use_CGL_pressure) then
                feedback_rhs(m,l,particle_tmp%i_elm,i_tor,1) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,1) &
                                                              + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                              * (0.5d0*p2_perp + b_norm_r**2 * p_atrop) / (gamma * sim%groups(1)%mass) * mu_zero !Pi_RR
                feedback_rhs(m,l,particle_tmp%i_elm,i_tor,2) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,2) &
                                                              + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                              * (0.5d0*p2_perp + b_norm_z**2 * p_atrop) / (gamma * sim%groups(1)%mass) * mu_zero !Pi_ZZ
                feedback_rhs(m,l,particle_tmp%i_elm,i_tor,3) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,3) &
                                                              + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                              * (0.5d0*p2_perp + b_norm_phi**2 * p_atrop) / (gamma * sim%groups(1)%mass) * mu_zero !Pi_phiphi
                feedback_rhs(m,l,particle_tmp%i_elm,i_tor,4) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,4) &
                                                              + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                              * (b_norm_r*b_norm_z * p_atrop) / (gamma * sim%groups(1)%mass) * mu_zero !Pi_RZ
                feedback_rhs(m,l,particle_tmp%i_elm,i_tor,5) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,5) &
                                                              + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                              * (b_norm_r*b_norm_phi * p_atrop) / (gamma * sim%groups(1)%mass) * mu_zero !Pi_Rphi
                feedback_rhs(m,l,particle_tmp%i_elm,i_tor,6) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,6) &
                                                              + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                              * (b_norm_z*b_norm_phi * p_atrop) / (gamma * sim%groups(1)%mass) * mu_zero !Pi_Zphi
              else
                ! Isotropic pressure: (1/3) p^2/(gamma m) on the diagonal,
                ! zero off-diagonal (components 4..6 remain 0). This keeps the loop
                ! correct also when use_pcs_full is enabled.
                feedback_rhs(m,l,particle_tmp%i_elm,i_tor,1) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,1) &
                                                              + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                              * (1.d0/3.d0) * ((p2_par + p2_perp) / gamma / sim%groups(1)%mass) * mu_zero !Pi_RR
                feedback_rhs(m,l,particle_tmp%i_elm,i_tor,2) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,2) &
                                                              + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                              * (1.d0/3.d0) * ((p2_par + p2_perp) / gamma / sim%groups(1)%mass) * mu_zero !Pi_ZZ
                feedback_rhs(m,l,particle_tmp%i_elm,i_tor,3) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,3) &
                                                              + HZ(i_tor) * v * particle_tmp%weight * atomic_mass_unit &
                                                              * (1.d0/3.d0) * ((p2_par + p2_perp) / gamma / sim%groups(1)%mass) * mu_zero !Pi_phiphi
              end if

            enddo !< tor harmonic
          enddo   !< order
        enddo     !< vertex
      endif
        
    end do ! n_steps substeps

    ! Write back persistent state for next call
    theta_prev_arr(j)      = theta_prev_p
    theta_accum_arr(j)     = theta_accum_p
    dtheta_prev_arr(j)     = dtheta_prev_p
    phi_prev_arr(j)        = phi_prev_p
    phi_unwrap_arr(j)      = phi_unwrap_p
    t_ref_a_arr(j)         = t_ref_a_p
    phi_ref_a_arr(j)       = phi_ref_a_p
    t_ref_b_arr(j)         = t_ref_b_p
    phi_ref_b_arr(j)       = phi_ref_b_p
    n_periods_arr(j)       = n_periods_p
    omega_theta_buf_arr(:,j) = omega_theta_buf_p
    omega_phi_buf_arr(:,j)   = omega_phi_buf_p

    call copy_particle(particles(j), particle_tmp)

    i_elm = particle_tmp%i_elm

    if (i_elm .gt. 0 .and. i_elm .le. sim%fields%element_list%n_elements) then

      ! Time after the last substep (= target_time for confined particles)
      t = particle_start_time + real(n_steps,8) * dt_local

      call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
      call mode_moivre(particle_tmp%x(3), HZ)
      call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
      B_mag = norm2(B)
      gamma = compute_relativistic_factor(particle_tmp, sim%groups(1)%mass, B_mag)

      ! Safety: skip feedback and projections for runaway particles.
      ! Normal EP have gamma ~ 1.00016 (300 keV D); gamma >= 100 indicates an
      ! unphysical particle that would destabilise GMRES and corrupt diagnostics.
      ! (The anisotropic pressure-tensor feedback is accumulated per substep
      ! inside the time-stepping loop above, mirroring tae_loop_CGL.)
      if (gamma >= 1.d2) cycle

      ! Compute power exchange only when the particle is still in the domain
      psiN = (psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis)
      if(psiN > 0.01 .and. psiN < 0.99 .and. midpoint_stored) then
        E_diff = -E_diff + rel_kinetic_energy(particle_tmp, sim%groups(1)%mass, B_mag)
        E_diff_arr(j) = E_diff

        ! Calculate gamma at midpoint for correct vparmid calculation
        call sim%fields%calc_EBpsiU(tmid, particle_mid%i_elm, particle_mid%st, particle_mid%x(3), E, B, psi, U)
        B_mag = norm2(B)
        gamma = compute_relativistic_factor(particle_mid, sim%groups(1)%mass, B_mag)

        ! psiN from midpoint (consistent with mumid, Pphi_mid)
        psiN_mid = (psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis)
        
        mumid   = particle_mid%p(2) * ATOMIC_MASS_UNIT / EL_CHG  ! in [eV/Tesla]
        vparmid = particle_mid%p(1) / (gamma * sim%groups(1)%mass) ! in [m/s]
        ! P_phi (canonical toroidal angular momentum) at midpoint, in [AMU·m²/s]
        !   P_phi = q*psi + p_par_SI * R * Bphi/B   (SI), Bphi = B(3) physical toroidal component
        Pphi_mid = ( particle_mid%x(1) * particle_mid%p(1) * ATOMIC_MASS_UNIT * (B(3)/B_mag) &
                     + dble(particle_mid%q) * EL_CHG * psi ) / ATOMIC_MASS_UNIT
        ! E = (gamma-1) m c^2 at midpoint, in [keV]
        E_mid_keV = (gamma - 1.d0) * sim%groups(1)%mass * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 / EL_CHG / 1.d3
        if (sim%my_id .eq. 0) then
          if(mod(j,500) .eq. 0) then 
            write(*,*) 'DEBUG: particle ', j, ' at psiN_mid= ', psiN_mid, ' mumid [eV/T] = ', mumid, ' Pphi_mid = ', Pphi_mid, ' E_mid_keV = ', E_mid_keV, ' E_diff [J] = ', E_diff, 'trap flag = ', trap_pass_arr(j)
          endif
        endif
        ! (psi_N, E, mu)
        if (output_pe_E_mu) then
        if(use_trap_passing_for_PE) then 
          if (trap_pass_arr(j) .eq. ORBIT_TRAPPED) then
              call project_single_particle_x(power_exchange_Emu_trap,[psiN_mid,E_mid_keV,mumid],phase_proj_trap_E,E_diff)
          else if (trap_pass_arr(j) .eq. ORBIT_PASSING) then
              call project_single_particle_x(power_exchange_Emu_pass,[psiN_mid,E_mid_keV,mumid],phase_proj_1_E,E_diff)
          end if
        else
            call project_single_particle_x(power_exchange_Emu,[psiN_mid,E_mid_keV,mumid],phase_proj_1_E,E_diff)
        end if ! use_trap_passing_for_PE 
        end if
        ! (psi_N, mu, P_phi)
        if (output_pe_mu_Pphi) then
          if(use_trap_passing_for_PE) then
            if (trap_pass_arr(j) .eq. ORBIT_TRAPPED) then
              call project_single_particle_x(power_exchange_vpar_mu_trap,[psiN_mid,mumid,Pphi_mid],phase_proj_trap,E_diff)
            else if (trap_pass_arr(j) .eq. ORBIT_PASSING) then
              call project_single_particle_x(power_exchange_vpar_mu_pass,[psiN_mid,mumid,Pphi_mid],phase_proj_1,E_diff)
            end if
          else
            call project_single_particle_x(power_exchange_vpar_mu,[psiN_mid,mumid,Pphi_mid],phase_proj_1,E_diff)
          end if ! use_trap_passing_for_PE
        end if
        ! (E, mu, P_phi)
        if (use_pe_EmuPphi) then
          if(use_trap_passing_for_PE) then
            if (trap_pass_arr(j) .eq. ORBIT_TRAPPED) then
              call project_single_particle_x(power_exchange_EmuPphi_trap,[E_mid_keV,mumid,Pphi_mid],phase_proj_trap_EP,E_diff)
            else if (trap_pass_arr(j) .eq. ORBIT_PASSING) then
              call project_single_particle_x(power_exchange_EmuPphi_pass,[E_mid_keV,mumid,Pphi_mid],phase_proj_1_EP,E_diff)
            end if
          else
            call project_single_particle_x(power_exchange_EmuPphi,[E_mid_keV,mumid,Pphi_mid],phase_proj_1_EP,E_diff)
          end if ! use_trap_passing_for_PE
        end if
  
      endif ! psiN and midpoint check for projection 
    endif ! i_elm > 0

  end do   ! particles
  !$omp end parallel do
  
  if (sim%my_id .eq. 0) write(*,*) "End of the particle loop"

  ! ---- Project resonance number L = (omega_0 - n*omega_phi) / omega_theta vs r ----
  ! Done serially after the OMP region to avoid reduction overhead on res_num arrays.
  do j = 1, size(particles, 1)
    if (particles(j)%i_elm .le. 0) cycle
    if (n_periods_arr(j) < MIN_PERIODS_RESONANCE) cycle
    
    ! r = sqrt((particles(j)%x(1) - ES%R_axis)**2 + (particles(j)%x(2) - ES%Z_axis)**2)
    call sim%fields%calc_EBpsiU(tmid, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
    psiN = (psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis)

    n_use_p = min(n_periods_arr(j), MIN_PERIODS_RESONANCE)
    omega_theta_avg_p = sum(omega_theta_buf_arr(1:n_use_p, j)) / dble(n_use_p)
    omega_phi_avg_p   = sum(omega_phi_buf_arr(1:n_use_p, j))   / dble(n_use_p)
    if (abs(omega_theta_avg_p) < 1.d-10) cycle
    L_res = (TWOPI * dble(TAE_FREQ) - dble(N_TOR_RESONANCE) * omega_phi_avg_p) / omega_theta_avg_p
    if (trap_pass_arr(j) .eq. ORBIT_TRAPPED) then
      ! call project_single_particle_x(res_num_trap, [r, L_res], res_num_trap%values, E_diff_arr(j)) ! for L vs r
      call project_single_particle_x(res_num_trap, [psiN, L_res], res_num_trap%values, E_diff_arr(j)) ! for L vs psiN
    else if (trap_pass_arr(j) .eq. ORBIT_PASSING) then
      ! call project_single_particle_x(res_num_pass, [psiN, L_res], res_num_pass%values, E_diff_arr(j)) ! for L vs r
      call project_single_particle_x(res_num_pass, [psiN, L_res], res_num_pass%values, E_diff_arr(j)) ! for L vs psiN
    end if
  end do
end select

! close(log_unit)

! Average the per-substep accumulation over the JOREK step (same as tae_loop_CGL)
jorek_feedback%rhs = feedback_rhs / n_steps
! (psi_N, mu, P_phi)
  if (output_pe_mu_Pphi) then
  if(use_trap_passing_for_PE) then
    power_exchange_vpar_mu_trap%values = power_exchange_vpar_mu_trap%values + phase_proj_trap
    power_exchange_vpar_mu_pass%values = power_exchange_vpar_mu_pass%values + phase_proj_1
  else
    power_exchange_vpar_mu%values = power_exchange_vpar_mu%values + phase_proj_1
  end if
end if
! (psi_N, E, mu)
if (output_pe_E_mu) then
  if(use_trap_passing_for_PE) then
    power_exchange_Emu_trap%values = power_exchange_Emu_trap%values + phase_proj_trap_E
    power_exchange_Emu_pass%values = power_exchange_Emu_pass%values + phase_proj_1_E
  else
    power_exchange_Emu%values = power_exchange_Emu%values + phase_proj_1_E
  end if
end if
! (E, mu, P_phi)
if (use_pe_EmuPphi) then
  if(use_trap_passing_for_PE) then
    power_exchange_EmuPphi_trap%values = power_exchange_EmuPphi_trap%values + phase_proj_trap_EP
    power_exchange_EmuPphi_pass%values = power_exchange_EmuPphi_pass%values + phase_proj_1_EP
  else
    power_exchange_EmuPphi%values = power_exchange_EmuPphi%values + phase_proj_1_EP
  end if
end if
deallocate(phase_proj_1)
deallocate(phase_proj_trap)
deallocate(phase_proj_1_E)
deallocate(phase_proj_trap_E)
deallocate(phase_proj_1_EP)
deallocate(phase_proj_trap_EP)
deallocate(feedback_rhs)
deallocate(E_diff_arr)

end subroutine


!> Write one line of trap/pass diagnostics for particle j at time t.
!> Columns: particle_index, time, poloidal_angle [rad], trap_pass_flag
!>   flag:  0 = undetermined,  1 = trapped,  -1 = passing
subroutine write_trap_pass_diag(unit, j, R, Z, flag, t)
  implicit none
  integer, intent(in) :: unit, j, flag
  real*8,  intent(in) :: R, Z, t
  real*8 :: theta
  theta = atan2(Z - ES%Z_axis, R - ES%R_axis) 
  write(unit, '(I8, 2(1X,ES16.8), 1X,I3)') j, t, theta, flag
end subroutine write_trap_pass_diag

subroutine extract_element_psi_minmax(fields,psi_minmax_list_1d)
  use mod_fields, only: fields_base
  implicit none
  !> inputs:
  class(fields_base),intent(in)  :: fields
  !> outputs:
  real*8,dimension(2*fields%element_list%n_elements),intent(out) :: psi_minmax_list_1d
  !> variables:
  integer :: ii
  real*8,dimension(fields%element_list%n_elements,2) :: psi_minmax_list
  !> initialisation
  psi_minmax_list(:,1) = 1d10; psi_minmax_list(:,2) = -1d10;
  !> extract maximum and minimun
  !$omp parallel do default(shared) private(ii)
  do ii=1,fields%element_list%n_elements
    call psi_minmax(fields%node_list,fields%element_list,ii,psi_minmax_list(ii,1),&
    psi_minmax_list(ii,2))
  enddo
  !$omp end parallel do
  psi_minmax_list_1d(1:fields%element_list%n_elements) = psi_minmax_list(:,1)
  psi_minmax_list_1d(fields%element_list%n_elements+1:2*fields%element_list%n_elements) = &
  psi_minmax_list(:,2)
end subroutine extract_element_psi_minmax

subroutine p_star_gamma0(theta, mass_AMU, p_star_out, gamma_0)
  implicit none
  real*8, intent(in) :: theta, mass_AMU
  real*8, intent(out) :: p_star_out, gamma_0
  real*8 :: x, y
  y = 2*theta*(theta + sqrt(1+theta**2))
  x = sqrt(y)
  p_star_out = mass_AMU*SPEED_OF_LIGHT*x
  gamma_0 = sqrt(1+x**2)
end subroutine p_star_gamma0


!> Classify a guiding-centre particle as trapped or passing and
!> compute omega_theta (poloidal/bounce frequency) and omega_phi (toroidal precession frequency).
!> All inout arguments are persistent per-particle state that must survive between calls.
subroutine classify_trap_pass(R_p, Z_p, phi_p, i_elm, t, &
                              trap_pass, theta_prev_p, theta_accum_p, dtheta_prev_p, &
                              phi_unwrap_p, phi_prev_p, &
                              t_ref_a_p, phi_ref_a_p, t_ref_b_p, phi_ref_b_p, &
                              n_periods_p, omega_theta_buf_p, omega_phi_buf_p)
  implicit none
  ! INPUTS and INOUTS
  real*8,  intent(in)    :: R_p, Z_p, phi_p !< particle (R, Z, phi) position
  integer, intent(in)    :: i_elm           !< element index (> 0 means inside mesh)
  real*8,  intent(in)    :: t               !< current time [s]
  integer, intent(inout) :: trap_pass       !< 0 = unknown, 1 = trapped, -1 = passing
  real*8,  intent(inout) :: theta_prev_p, theta_accum_p, dtheta_prev_p
  real*8,  intent(inout) :: phi_unwrap_p, phi_prev_p
  real*8,  intent(inout) :: t_ref_a_p, phi_ref_a_p   !< upper-turn (trapped) / crossing (passing)
  real*8,  intent(inout) :: t_ref_b_p, phi_ref_b_p   !< lower-turn (trapped only)
  integer, intent(inout) :: n_periods_p               !< completed-period counter
  real*8,  intent(inout) :: omega_theta_buf_p(MIN_PERIODS_RESONANCE) !< circular buffer of per-period omega_theta
  real*8,  intent(inout) :: omega_phi_buf_p(MIN_PERIODS_RESONANCE)   !< circular buffer of per-period omega_phi
  ! Local variables
  real*8  :: theta_curr_p, dtheta_p, dphi, T_per, dphi_per
  integer :: old_class, idx_buf

  if (i_elm .le. 0) return

  ! ---- Phase 0: poloidal angle theta(t) = atan2(Z-Z_axis, R-R_axis) ----
  theta_curr_p = atan2(Z_p - ES%Z_axis, R_p - ES%R_axis)  ! range (-pi, pi]
  dtheta_p = theta_curr_p - theta_prev_p
  ! ---- Wrap dtheta to remove 2pi jumps ----
  if (dtheta_p >  0.5d0*TWOPI) dtheta_p = dtheta_p - TWOPI
  if (dtheta_p < -0.5d0*TWOPI) dtheta_p = dtheta_p + TWOPI

  ! ---- Wrap dphi to remove 2pi jumps ----
  dphi = phi_p - phi_prev_p
  if (dphi >  0.5d0*TWOPI) dphi = dphi - TWOPI
  if (dphi < -0.5d0*TWOPI) dphi = dphi + TWOPI
  phi_unwrap_p = phi_unwrap_p + dphi

  old_class = trap_pass

  ! ---- Phase 1 + 2 + 3: classify & compute frequencies ----
  if (dtheta_prev_p * dtheta_p < 0.d0) then
    ! Sign reversal in dtheta -> trapped (banana orbit)
    trap_pass     = ORBIT_TRAPPED
    theta_accum_p = 0.d0

    ! Reset frequency state on classification change
    if (old_class .ne. ORBIT_TRAPPED) then
      t_ref_a_p = 0.d0;  t_ref_b_p = 0.d0
      n_periods_p = 0;  omega_theta_buf_p = 0.d0;  omega_phi_buf_p = 0.d0
    end if

    ! Turning-point type: dtheta_prev > 0 -> was moving upward -> upper turning point
    if (dtheta_prev_p > 0.d0) then
      if (t_ref_a_p > 0.d0) then
        T_per     = t - t_ref_a_p
        dphi_per  = phi_unwrap_p - phi_ref_a_p
        if (T_per > 0.d0) then
          n_periods_p = n_periods_p + 1
          idx_buf = mod(n_periods_p - 1, MIN_PERIODS_RESONANCE) + 1
          omega_theta_buf_p(idx_buf) = TWOPI / T_per
          omega_phi_buf_p(idx_buf)   = dphi_per / T_per
        end if
      end if
      t_ref_a_p = t;  phi_ref_a_p = phi_unwrap_p
    else
      ! Lower turning point
      if (t_ref_b_p > 0.d0) then
        T_per     = t - t_ref_b_p
        dphi_per  = phi_unwrap_p - phi_ref_b_p
        if (T_per > 0.d0) then
          n_periods_p = n_periods_p + 1
          idx_buf = mod(n_periods_p - 1, MIN_PERIODS_RESONANCE) + 1
          omega_theta_buf_p(idx_buf) = TWOPI / T_per
          omega_phi_buf_p(idx_buf)   = dphi_per / T_per
        end if
      end if
      t_ref_b_p = t;  phi_ref_b_p = phi_unwrap_p
    end if

  else
    theta_accum_p = theta_accum_p + dtheta_p
    if (abs(theta_accum_p) >= TWOPI) then
      ! Completed a full poloidal transit -> passing
      trap_pass     = ORBIT_PASSING
      theta_accum_p = theta_accum_p - sign(TWOPI, theta_accum_p)

      ! Reset frequency state on classification change
      if (old_class .ne. ORBIT_PASSING) then
        t_ref_a_p = 0.d0
        n_periods_p = 0;  omega_theta_buf_p = 0.d0;  omega_phi_buf_p = 0.d0
      end if

      ! Transit crossing event
      if (t_ref_a_p > 0.d0) then
        T_per     = t - t_ref_a_p
        dphi_per  = phi_unwrap_p - phi_ref_a_p
        if (T_per > 0.d0) then
          n_periods_p = n_periods_p + 1
          idx_buf = mod(n_periods_p - 1, MIN_PERIODS_RESONANCE) + 1
          omega_theta_buf_p(idx_buf) = TWOPI / T_per
          omega_phi_buf_p(idx_buf)   = dphi_per / T_per
        end if
      end if
      t_ref_a_p = t;  phi_ref_a_p = phi_unwrap_p
    end if
  end if

  dtheta_prev_p = dtheta_p
  theta_prev_p  = theta_curr_p
  phi_prev_p    = phi_p
end subroutine classify_trap_pass

end program tae_loop
