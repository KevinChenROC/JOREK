!> Testing the coupling of the projections of particles to JOREK
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
use phys_module, only: tstep, restart, t_start, restart_particles, nout, nout_projection
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint, amin, F0, R_geo
use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
use phys_module, only: n_mode_families

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
use mod_import_experimental_dist, only: calculate_B
use mod_gc_variational, only: convert_gc_to_gc_vpar, convert_leapfrog_to_gc_vpar
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
type(phase_space_projection)                      :: fourD_dist, power_exchange_vpar_mu, RZ_dist, PparMu_dist

real*8, parameter  :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)
real*8    :: target_time
real*8    :: physical_particles, weight, tstep_keep
real*8    :: oldtime, step_rest_time, particle_step_time, particle_start_time, diag_time
real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, tstep_si, timesteps
real*8    :: v_kin_temp, E(3), B(3), psi, U, B_norm(3)
real*8    :: rescale_coef, T_axis(1), E_axis, E_hot, rho_part, v2, tstart_jorek
real*8    :: dummy_real8_1, dummy_real8_2, dummy_real8_3, kT_mc2_ratio, p_star, gamma_star
real*8, allocatable :: real_gdf_param(:), real_pdf_param(:), psi_minmax_loc(:), phase_space_bounds(:,:), uniform_var(:)
real*8, parameter :: &
                      Poloidalbound(2) = [0.d0, TWOPI], &
                      Phibound(2)      = [0.d0, TWOPI], &
                      Pbound(2)        = [1.5d4, 2.5d7], & ! [AMU·m/s]
                      T_EP_eV          = 1.0d5      ! EP temperature in eV
integer, parameter :: &
                      N_PHI_PLANES         = 12,  & ! particles are spread over TWOPI * i/N_PHI_PLANES in phi to reduce noise when initial mode structure is used
                      CHARGE_EP            = 1       ! EP charge in units of elementary charge
!$ real*8 :: w0, w1, mmm(3)
integer   :: n_real_gdf_param, n_int_gdf_param, n_int_pdf_param, n_real_pdf_param, n_variables, RZ_OUTPUT_STEPS, POWER_EXCHANGE_STEPS, PARTICLE_OUTPUT_STEPS
integer   :: n_particles_local,n_reflect,ifail,dummy_int, ino
integer   :: i, j, k, l, m, n_steps, i_elm_old, i_diagno
integer   :: seed, i_rng, n_stream
character(len=50) :: part_output_filename

! Start up MPI, jorek
call sim%initialize(num_groups=1)

rho_part    = 1.195d19 !(corrected value to obtain density=1.441e17 (as in benchmark, for original profile with toroidal flux) 

n_particles_local = int(n_particles/sim%n_cpu) 
timesteps         = tstep_particles
tstep_keep       = tstep

POWER_EXCHANGE_STEPS =  nout_projection
RZ_OUTPUT_STEPS      =  nout_projection*5
PARTICLE_OUTPUT_STEPS = nout*10

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1)) ! If a restart file with the n=6 mode structure is used, divide to obtain the mode at noise-level
! fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1,mode_divisor=200)) ! If a restart file with the n=6 mode structure is used, divide to obtain the mode at noise-level
call with(sim, fieldreader)

tstep = tstep_keep
write(*,*) 'main : t_start = ',t_start

if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)

call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

call update_equil_state(sim%my_id,sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = central_mass * ATOMIC_MASS_UNIT * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek

tstep_si  = tstep * t_norm
n_steps   = floor(tstep_si / timesteps)
timesteps = tstep_si / n_steps
n_steps   = tstep_si / timesteps

if (sim%my_id .eq.0) then
  write(*,*) ' Nonrelativistic GC particle loop with JOREK coupling'
  write(*,*) ' adapt time step to be multiple of jorek time step'
  write(*,*) "tstep = ", tstep_si, n_steps, timesteps
  write(*,*) "check :", n_steps, tstep_si - n_steps*timesteps

  i_diagno =  sim%fields%node_list%n_nodes / 3 
  write(*,'(A,6f8.4)') ' probe at : ',sim%fields%node_list%node(i_diagno)%x(1,1,1:2)
  open(111,file='diagno.txt')
endif

if (.not. restart_particles) then
  ! Set up particles
  sim%groups(1)%Z    = 1
  sim%groups(1)%mass = atomic_weights(-2) !< atomic mass units
                    
  if(sim%my_id .eq. 0) write(*,*) "phase_space_bounds(:,1):", phase_space_bounds(:,1), "\n (:,2):", phase_space_bounds(:,2)
  allocate(particle_gc_vpar::sim%groups(1)%particles(n_particles_local))
  ! There are logical errors if I use initialise_particles_H_mu_psi_phiplanes, when N_PHI_PLANES > 1, so I use this subroutine and initialise phi separately
  call initialise_particles_H_mu_psi(sim%groups(1)%particles, sim%fields, pcg32_rng(),sim%groups(1)%mass, &
  uniform_space=.true., uniform_space_rej_f=f_toroidal_flux, &
  uniform_space_rej_vars=[1], charge = CHARGE_EP, T_maxwell = T_EP_eV)
  
  
  call adjust_particle_weights(sim%groups(1)%particles, rho_part)
  if (sim%my_id .eq. 0) write(*,*) "Particle density was adjusted to:", rho_part, sim%groups(1)%particles(1:10)%weight

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
  endif
  deallocate(uniform_var)

  select type (part=>sim%groups(1)%particles)
    type is (particle_gc_vpar)
    do i=1, 300 
        ! print particle's info for debugging
        if (sim%my_id .eq. 0) write(*,*) ' particle ', i, &
        ' phi = ', part(i)%x(3), & 
        ", theta = ", atan2(part(i)%x(2) - ES%Z_axis, part(i)%x(1)-ES%R_axis), & 
        " q = ", part(i)%q
    enddo
    if(sim%my_id .eq. 0) then
      do i=1, size(part,1)
        if(part(i)%mu < 0.d0) then
          write(*,*) ' Warning: negative mu encountered during particle initialisation. Setting mu=0 for particle ', i
          part(i)%mu = 0.d0
        endif
      end do
    endif
  end select
  ! Write initial particle distribution to file
  partwriter = event(write_action(filename='part_restart.h5'))
  call with(sim, partwriter)
  
else  ! restarting particles

  if (sim%my_id .eq. 0) write(*,*) 'restarting particles: reading part_restart.h5'

  deallocate(sim%groups)
  allocate(sim%groups(0))

  partreader = event(read_action(filename='part_restart.h5'))
  call with(sim, partreader)

endif

!do i=1,sim%fields%node_list%n_nodes
!  sim%fields%node_list%node(i)%values(2:3,:,:) = 1.d-2 * sim%fields%node_list%node(i)%values(2:3,:,:)
!enddo

jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                filter    = filter_perp, filter_hyper    = filter_hyper, filter_parallel    = filter_par, &
                                filter_n0 = filter_perp, filter_hyper_n0 = filter_hyper, filter_parallel_n0 = filter_par_n0, &
                                calc_integrals=.false., to_vtk=.false., to_h5 = .false., basename='projections')

allocate(jorek_feedback%rhs(n_degrees, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 1))

jorek_feedback%rhs = 0.d0

aux_node_list => jorek_feedback%node_list

! Set the EP temperature for projection PparMu normalization
call set_T_EP_eV(T_EP_eV)

! Set up 4D phase space projection: R, Z, P (normalized velocity), Energy
RZ_dist = new_phase_space_projection(ndim=2,res=[91, 91],start=[9.d0, -1.d0], end=[11.d0, 1.d0], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_RZ, group=1),basename='RZ_dist')
! fourD_dist = new_phase_space_projection(ndim=4,res=[31,31,41,81],start=[9.d0,-1.d0,-1.1d0, -20.d0],end=[11.0,1.d0,1.1d0,1700.d0], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_RZPE, group=1),basename='fourD_dist')
! fourD_dist = new_phase_space_projection(ndim=2,res=[51,81],start=[-1.1d0, -20.d0],end=[1.1d0,1700.d0], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_PparE, group=1),basename='PparE_dist')
PparMu_dist = new_phase_space_projection(ndim=2,res=[301,301],start=[-0.045, 0d0],end=[0.045, 2.5d2], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_PparMu, group=1), basename='PparMu_dist') ! project particles on the grid (1) ppar/mc (2) mu [keV/T]
call with(sim,PparMu_dist)
call output_phase_project(PparMu_dist,0,output_grids_in=.true.)
call with(sim,RZ_dist)
call output_phase_project(RZ_dist,0,output_grids_in=.true.)

! The output is only done on the root process. To prevent a too big imbalance, barrier here.
call MPI_BARRIER(MPI_COMM_WORLD,ifail)


! Set up 2D power exchange: v_parallel (in m/s) vs mu (magnetic moment in eV/Tesla)
power_exchange_vpar_mu = new_phase_space_projection(ndim=2,res=[300,300],start=[-2.d7,-0.1d6],end=[2.d7,1.5d6],basename="power_exchange",bandwidths=[0.15d7,0.08d6])
power_exchange_vpar_mu%values = 0.d0
call output_phase_project(power_exchange_vpar_mu, 0, output_grids_in=.true.)

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
if (.not. restart_particles) call with(sim, events, at=sim%time)

!call with(sim, project_density)

do while (.not. sim%stop_now)

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
  
  call loop_particle_gc(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time, power_exchange_vpar_mu)
  
  if(sim%my_id .eq. 0) write(*,*) "Finished particle loop for ino=", ino

  ino = ino + 1
  ! Output power exchange projection
  if(mod(ino,POWER_EXCHANGE_STEPS) .eq. 0) then ! a few more fluid steps to accumulate power exchange
    call output_phase_project(power_exchange_vpar_mu, ino + index_start, output_grids_in=.false.)
    power_exchange_vpar_mu%values = 0.d0 
  endif
  
  ! Output 4D distribution function periodically
  if(mod(ino,RZ_OUTPUT_STEPS) .eq. 0) then 
    call with(sim,RZ_dist)
    call output_phase_project(RZ_dist, ino + index_start, output_grids_in=.false.)
    call MPI_BARRIER(MPI_COMM_WORLD,ifail)
  endif

  sim%time = target_time 
  
  call with(sim, events, at=sim%time)

  ! Save particles every PARTICLE_OUTPUT_STEPS timesteps AFTER JOREK timestepping
  ! This ensures particles are consistent with the JOREK fields at this timestep
  if (mod(ino, PARTICLE_OUTPUT_STEPS) .eq. 0) then
    write(part_output_filename, '(A,I6.6,A)') 'part_', ino + index_start, '.h5'
    if (sim%my_id .eq. 0) write(*,*) 'Saving particles to: ', trim(part_output_filename)
    partwriter = event(write_action(filename=trim(part_output_filename)))
    call with(sim, partwriter)
  endif

  if (sim%my_id == 0) write(111,'(12e20.12)') sim%time, sim%fields%node_list%node(i_diagno)%values(1:n_tor,1,1:2) 

end do

!call write_simulation_hdf5(sim, 'part_restart.h5')

! partwriter = event(write_action(filename='part_restart.h5'))
call with(sim, partwriter)

call sim%finalize

if (sim%my_id == 0) close(111)

contains

subroutine loop_particle_gc(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time, test_phase)
use mod_project_particles
use mod_random_seed
use mod_interp, only: mode_moivre
use mod_basisfunctions
use mod_particle_types, only: copy_particle
use mod_gc_variational, only: push_gc_rk4

implicit none

class(particle_sim), target, intent(inout)                :: sim
type(projection), target, intent(inout)                   :: jorek_feedback
type(phase_space_projection), target, intent(inout)       :: test_phase
type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
type(particle_gc_vpar)                                    :: particle_tmp

real*8, intent(in)     :: timesteps, particle_start_time 
real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
real*8    :: t, E(3), B(3), psi, U, n_e, T_e, rz_old(2), st_old(2), r, v, tmp2(2)
real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
real*8    :: xcontainer(n_steps), mucontainer(n_steps), pparcontainer(n_steps), vparcontainer(n_steps), Econtainer(n_steps)
real*8    :: gamma, p_mag, t0,  tend, tmid
real*8    :: fitsinevpar(4), fitsinemu(4), fitsineE(4), fitsineB(4), omega, mumid, vparmid, pparmid, mu_start, vpar_start, B_start, mu_end, vpar_end, B_end, E_diff, E_i, avgB, B_mag
real*8, parameter :: AMU_C2_EV = ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 / EL_CHG ! rest-mass energy of 1 amu in eV. 
!$ real*8 :: w0, w1, mmm(3)

integer, intent(in)   :: n_steps
integer   :: i, j, k, l, m, i_elm_old, i_elm 
integer   :: seed, n_stream, ierr, nthreads
integer   :: i_tor, index_lm, i_elm_temp
integer   :: n_particles, ifail
real*8,allocatable :: feedback_rhs(:,:,:,:,:)
real*8, allocatable :: phase_proj(:)

!$ w0 = omp_get_wtime()

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = central_mass * ATOMIC_MASS_UNIT * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation

jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps
allocate(feedback_rhs,source=jorek_feedback%rhs)
allocate(phase_proj,source=test_phase%values)

jorek_feedback%rhs = 0.d0
feedback_rhs       = 0.d0

! Calculate time values before parallel region
t0   = timesteps
tend = n_steps * timesteps
tmid = 0.5d0*(t0+tend)
xcontainer = (/(i,i=1,n_steps)/) ! treated as constants

select type (particles => sim%groups(1)%particles)
type is (particle_gc_vpar)
#ifdef __GFORTRAN__
 !$omp parallel do default(shared) & 
#else
 !$omp parallel do default(none) &
 !$omp shared(sim, particles, n_steps, timesteps, rng, particle_start_time,        &
 !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, t0, tmid, tend,           &
 !$omp jorek_feedback, CENTRAL_DENSITY, CENTRAL_MASS, test_phase, xcontainer)      &
#endif
 !$omp private(particle_tmp, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old, r, v, tmp2,    &
 !$omp i_elm_old, i_elm, n_e, T_e, E_diff, avgB, B_mag,                                               & 
 !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, ifail, gamma, p_mag,         &
 !$omp fitsinevpar, fitsinemu, fitsineE, fitsineB, omega, vparmid, pparmid, mumid, mu_start, vpar_start, B_start, mu_end, vpar_end, B_end, &
 !$omp mucontainer, pparcontainer, vparcontainer, Econtainer) & 
 !$omp schedule(dynamic,10) &
 !$omp reduction(+:feedback_rhs)&
 !$omp reduction(+:phase_proj)
  do j=1,size(particles,1)
    call copy_particle(particle_tmp, particles(j))

    ! Lost particles (i_elm < 0) and Uninitialized/error particles (i_elm = 0)
    if (particle_tmp%i_elm .le. 0) cycle ! skip to next particle

    call sim%fields%calc_EBpsiU(particle_start_time, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
    B_mag = norm2(B)
    E_diff = particle_energy(particle_tmp, B_mag, sim%groups(1)%mass) ! initial energy in J

    ! push from start to mid
    call push_gc_rk4(sim%fields, particle_tmp, sim%groups(1)%mass, timesteps, n_steps/2, n_gyro_phases=0, time_start=particle_start_time) ! no gyro-averaging. No FLR effect

    if (particle_tmp%i_elm .le. 0) cycle ! skip to next particle if lost during push

    vparmid = particle_tmp%vpar ! in [m/s]
    mumid   = abs(particle_tmp%mu) * (sim%groups(1)%mass * ATOMIC_MASS_UNIT) / EL_CHG ! convert from [m^2/s^2/T] to [eV/T]
    
    ! push from mid to end (use n_steps - n_steps/2 to handle odd n_steps correctly)
    call push_gc_rk4(sim%fields, particle_tmp, sim%groups(1)%mass, timesteps, n_steps - n_steps/2, n_gyro_phases=0, time_start=particle_start_time + (n_steps/2)*timesteps)

    call copy_particle(particles(j), particle_tmp)

    i_elm = particle_tmp%i_elm

    if (i_elm .gt. 0 .and. i_elm .le. sim%fields%element_list%n_elements) then
      call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
      call mode_moivre(particle_tmp%x(3), HZ)
      ! t = particle_start_time + (n_steps)*timesteps
      ! call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
      B_mag = particle_tmp%B_norm
      
      
      do l=1,n_vertex_max
        do m=1,n_degrees

          index_lm = (l-1)*n_degrees + m

          v = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m) 

          do i_tor=1,n_tor
            feedback_rhs(m,l,i_elm,i_tor,1) = feedback_rhs(m,l,i_elm,i_tor,1) &
                                                  
                                                  + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * atomic_mass_unit &
            
                                                  * (1.d0/3.d0) * (particle_tmp%vpar**2 + 2 * abs(particle_tmp%mu) * B_mag) * mu_zero
          enddo

        enddo   !< order
      enddo     !< vertex

      ! Project only when the particle is still in the domain and near the n=6 TAE resonance surface r=0.5
      r = sqrt((particle_tmp%x(1)-ES%R_axis)**2 + (particle_tmp%x(2) - ES%Z_axis)**2)
      if(r>0.3 .and. r < 0.7) then
        E_diff = -E_diff + particle_energy(particle_tmp, B_mag, sim%groups(1)%mass) ! energy difference in J
        ! if (sim%my_id .eq. 0) then
        !   if(mod(j,2000) .eq. 0) write(*,*) 'DEBUG: particle ', j, ' at r= ', r, ' mumid [eV/T] = ', mumid, ' vparmid [m/s] = ', vparmid, ' E_diff [J] = ', E_diff
        ! endif
        call project_single_particle_x(test_phase,[vparmid,mumid],phase_proj,E_diff)
      endif
    endif ! i_elm > 0

  end do   ! particles
  !$omp end parallel do
  
  if (sim%my_id .eq. 0) write(*,*) "End of the particle loop"

end select

jorek_feedback%rhs = feedback_rhs
test_phase%values = test_phase%values + phase_proj
deallocate(phase_proj)
deallocate(feedback_rhs)

end subroutine

pure function f_adapted(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:4)
  real*4 :: f

  coeff(0)=0.53
  coeff(1)=0.3
  coeff(2)=0.2
  coeff(3)=0.52
  coeff(4)=0.26

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

  f = (f - coeff(4)) / (1.d0 - coeff(4))

end function f_adapted

pure function f_original(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_original

pure function f_toroidal_flux(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 :: s, psi_norm, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  psi_norm = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  s = 0.957 * psi_norm + 0.043 * psi_norm**2 

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_toroidal_flux

pure function f_toroidal_flux_Lu2023(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 :: s, psi_norm, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.46623
  coeff(1)=0.17042
  coeff(2)=0.11357
  coeff(3)=0.521298

  psi_norm = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  s = 0.957 * psi_norm + 0.043 * psi_norm**2 

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_toroidal_flux_lu2023

pure function particle_energy(particle, B_mag, mass_AMU)
  implicit none
  type(particle_gc_vpar), intent(in) :: particle
  real*8, intent(in) :: B_mag ! [T]
  real*8, intent(in) :: mass_AMU ! [AMU]
  real*8 :: particle_energy

  ! Here mu is 0.5vper^2/B in [m^2/s^2/T]
  particle_energy = 0.5d0 * mass_AMU * ATOMIC_MASS_UNIT * (particle%vpar**2 + 2 * abs(particle%mu) * B_mag ) ! in J
  
end function particle_energy

end program tae_loop
