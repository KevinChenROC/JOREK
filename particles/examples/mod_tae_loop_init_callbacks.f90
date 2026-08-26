module mod_tae_loop_init_callbacks
  use constants, only: SPEED_OF_LIGHT, TWOPI, ATOMIC_MASS_UNIT, EL_CHG
  use equil_info
  use mod_fields, only: fields_base
  use mod_particle_types, only: particle_base, particle_gc_relativistic, particle_kinetic_relativistic
  use phys_module, only: pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness,T_EP_eV, F0, RHO_EP, Lambda_peak, delta_Lambda, R_geo
  use mod_bessel, only: bessel_k2exp
  use mod_sampling, only: c_erfinv

  implicit none
  private
  public gdf_uniform, gdf_sqrtPsiN_p_sampler, gdf_RZPhiP_sampler_cb, pdf_RZ_p_HL2A_flat_decay, pdf_r_p_wang2020, pdf_sqrtPsiN_p_HL2A_parametric, pdf_psiN_p_Relativistic_MJ, pdf_psiN_E_parametric, gdf_PsiNThetaPhi_EorP_sampler, pdf_RZ_p_relativistic_MJ, pdf_psiN_p_wang2020
  public RZPhiP_to_gc_relativistic, RZPhiPparMu_to_gc_relativistic_cb, RZPhi_to_gc_relativistic_temp_gradient, PsiNThetaPhiP_to_gc_relativistic_cb, PsiNThetaPhiP_to_kinetic_relativistic_cb
  public sqrtPsiN_Theta_Phi_P_to_gc_relativistic, PsiNThetaPhiE_to_gc_relativistic_cb, particle_weight_one, log_f_MJ, gamma_relativistic
  public pdf_E_Pphi, gdf_RZPhiPparMu_sampler

  real*8, parameter :: LAMBDA_PEAK_MIN = 0.01d0 ! minimum Lambda_peak to use the peak distribution, otherwise use isotropic distribution
contains

function gamma_relativistic(ppar, mu, Bnorm, mass_AMU) result(out)
  implicit none
  ! Input
  real*8, intent(in) :: ppar ! parallel momentum in AMU·m/s
  real*8, intent(in) :: mu   ! magnetic moment in AMU·m^2/(T·s^2)
  real*8, intent(in) :: Bnorm    ! magnetic field in T
  real*8, intent(in) :: mass_AMU ! mass in AMU
  ! Output
  real*8 :: out
  ! Variables
  real*8 :: p_perp_sqr, p_total
  p_perp_sqr = 2.d0 * mass_AMU * mu * Bnorm
  p_total = sqrt(ppar**2 + p_perp_sqr)
  out = sqrt( (p_total/(mass_AMU*SPEED_OF_LIGHT))**2 + 1.d0 )
end function


function log_f_MJ(p, theta, mass_AMU) result(out)
  implicit none
  real*8, intent(in) :: p ! p in AMU·m/s
  real*8, intent(in) :: theta ! dimensionless temperature kB*T/(m*c^2)
  real*8, intent(in) :: mass_AMU ! mass in AMU
  real*8 :: out
  real*8 :: g
  ! p in AMU·m/s, m in AMU
  ! log of unnormalized MJ: log f = 2 ln p - gamma/θ
  g = sqrt( (p/(mass_AMU*SPEED_OF_LIGHT))**2 + 1.d0 )
  out = 2*log(p) - g/theta
end function log_f_MJ

!> Uniform distribution in (psi, theta, phi, p) or (R, Z, phi, p)
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random state to accept
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   fields:       (fields_base) jorek MHD fields
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution weight (normally mass/volume)
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   gdf: (real8) value of the sampler probability density 
function gdf_uniform(n_x,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(gdf)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: n_x,i_elm
  integer,intent(in)                          :: n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(n_x),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)               :: fields
  !> Outputs:
  real*8 :: gdf
  gdf = 1d0
end function gdf_uniform

subroutine gdf_RZPhiP_sampler_cb(n_x,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param)
  use mod_fields, only: fields_base
  implicit none
  integer,intent(in) :: n_x,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(n_x),intent(in)            :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)               :: fields
  integer,intent(inout)               :: i_elm
  real*8,dimension(2),intent(inout)   :: st
  real*8,dimension(n_x),intent(inout) :: x
  integer :: ifail
  real*8  :: DUMMY_R, DUMMY_Z

  x(1:n_x) = x_min(1:n_x)+(x_max(1:n_x)-x_min(1:n_x))*x(1:n_x)
  if (i_elm .le. 0 .or. i_elm .gt. fields%element_list%n_elements) i_elm = 1

  call find_RZ(fields%node_list, fields%element_list, x(1), x(2), DUMMY_R, DUMMY_Z, i_elm, st(1), st(2), ifail) ! find s,t,i_elm from R,Z
  if (ifail .ne. 0) then
    i_elm = 0
    st    = -1.d0
  endif
end subroutine gdf_RZPhiP_sampler_cb

subroutine RZPhiP_to_gc_relativistic(p_inout,n_x,x,time,fields,n_real_param, real_param,n_int_param,int_param)
  use mod_particle_types, only: particle_base
  use mod_fields,         only: fields_base
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x
  integer,intent(in)               :: n_int_param,n_real_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                :: time
  real*8,dimension(n_x),intent(in) :: x
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)      :: fields
  !> Variables:
  real*8 :: p, u_beta, beta, p_perp, u_ppar_sign, E(3), B(3), psi, electric_potential, normB, mass_AMU, DUMMY_R, DUMMY_Z, p_par, mu, Lambda_max, erf_a, erf_b, KE_eV, B0, gamma, Lambda
  integer :: ifail, i_elm
  !> inputs-outputs:
  class(particle_base),intent(inout) :: p_inout

  ! --- pre setup ---
  p = x(4) ! in [AMU·m/s]
  mass_AMU = real_param(2*fields%element_list%n_elements+3)
  u_beta = 2*mod(x(5), 0.5d0) ! Use the same uniform(0,1) for computing the beta(1,1/2) and the sign for p_par.
  u_ppar_sign = x(5) - 0.5d0
  B0 = abs(F0/R_geo)
  ! --- pre setup ---
  
  call fields%calc_EBpsiU(time,p_inout%i_elm,p_inout%st,p_inout%x(3),E,B,psi,electric_potential); normB = norm2(B);

  !! Compute mu and p_par
  if(Lambda_peak > LAMBDA_PEAK_MIN) then
    ! Maximum allowed Lambda at this location: ensures μB < E
    Lambda_max = 0.995d0 * B0 / normB
    if (Lambda_max <= 0.d0) then
      p_inout%i_elm = 0
      p_inout%st = -1.d0
      return
    end if
    gamma = sqrt(1 + (p/(mass_AMU*SPEED_OF_LIGHT))**2) ! relativistic factor
    KE_eV = (gamma - 1) * mass_AMU * (ATOMIC_MASS_UNIT / EL_CHG * SPEED_OF_LIGHT**2)  ! in [eV]
    ! Sample Λ from truncated Normal(Λ_peak, (ΔΛ/√2)²) with support [0, Λ_max]
    ! via inverse transform of the truncated Normal CDF:
    ! Λ = μ + σ√2 · erfinv( erf((a-μ)/(σ√2)) + u·(erf((b-μ)/(σ√2)) - erf((a-μ)/(σ√2))) )
    ! where μ=Λ_peak, σ=ΔΛ/√2, a=0, b=Λ_max, σ√2=ΔΛ
    erf_a = erf(-Lambda_peak / delta_Lambda)
    erf_b = erf((Lambda_max - Lambda_peak) / delta_Lambda)
    Lambda = Lambda_peak + delta_Lambda * c_erfinv(erf_a + u_beta * (erf_b - erf_a))
    ! Λ is guaranteed to be in [0, Λ_max] — no rejection needed
    mu = Lambda * (KE_eV * EL_CHG / ATOMIC_MASS_UNIT) / B0 ! in [AMU·m^2/(T·s^2)]
    ! Relativistic p_par: p_par^2 = p_total^2 - p_perp^2 = p^2 - 2*mass_AMU*mu*normB
    p_par = sign(sqrt(max(p**2 - 2.d0 * mass_AMU * mu * normB, 0.d0)), u_ppar_sign)
    if(abs(p_par) <= 1.d-8) write(*,*) "Warning: |p_par| <= 1.d-8, KE [eV] = ", KE_eV, " mu B [eV] = ", (mu * normB * ATOMIC_MASS_UNIT / EL_CHG), " p_par [AMU·m/s] = ", p_par
  else !isotropic momentum distribution
    beta = 2.d0*u_beta - u_beta**2
    p_par = sign(p*sqrt(1.d0 - beta), u_ppar_sign)
    ! compute the magnetic moment [AMU·m^2/(T·s^2)]
    p_perp = p * sqrt(beta)
    mu = p_perp**2 / (2.d0 * mass_AMU * normB)  
  endif


  select type( particle => p_inout)
    type is (particle_gc_relativistic)
      ! compute the parallel momentum [AMU·m/s]
      particle%p(1) = p_par      
      particle%p(2) = mu
      if(allocated(int_param)) then 
        particle%q    = int_param(1)
      else 
        write(*,*) "Error: int_param not allocated in sample_to_gc_relativistic for particle charge assignment"
      endif 
    class default
      write(*,*) "Error: particle type not supported in sample_to_gc_relativistic"
      stop
  end select

end subroutine RZPhiP_to_gc_relativistic

subroutine gdf_sqrtPsiN_p_sampler(n_x,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param)
  use constants,          only: MU_ZERO,EL_CHG,ATOMIC_MASS_UNIT
  use phys_module,        only: central_density
  use mod_model_settings, only: var_T,var_Vpar
  use mod_sampling,       only: sample_chi_squared_3
  use mod_fields,         only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in) :: n_x,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(n_x),intent(in)            :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)               :: fields
  !> Inputs-Outputs:
  integer,intent(inout)               :: i_elm
  real*8,dimension(2),intent(inout)   :: st
  real*8,dimension(n_x),intent(inout) :: x
  !> Variables:
  real*8 :: psiN,psi,normB,electric_potential,u,temperature_ev
  real*8 :: R,R_s,R_t,Z,Z_s,Z_t
  real*8,dimension(1) :: P,P_s,P_t,P_phi
  real*8,dimension(3) :: B,E

  !> Sample (sqrt_PsiN,theta,phi,p) from the uniform distribution
  x(1:4) = x_min(1:4)+(x_max(1:4)-x_min(1:4))*x(1:4)
  psiN = x(1)**2
  psi  = psiN * (ES%Psi_bnd - ES%Psi_axis) + ES%Psi_axis


  call find_theta_psi(fields%node_list,fields%element_list,& 
  [real_param(1:fields%element_list%n_elements), real_param(fields%element_list%n_elements+1:2*fields%element_list%n_elements)],&
  x(2), psi, x(3),&
  real_param(2*fields%element_list%n_elements+1), real_param(2*fields%element_list%n_elements+2),&
  i_elm,st(1),st(2),R,Z)
end subroutine gdf_sqrtPsiN_p_sampler


function pdf_sqrtPsiN_p_HL2A_parametric(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  !! Parametric skewed-normal PDF using namelist parameters (pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness).
  !! Allows scanning s_max without recompilation.
  use mod_fields, only: fields_base
  implicit none
  integer,intent(in)                          :: nx,i_elm,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  real*8 :: pdf
  pdf = pdf_RZ_p_HL2A_template(nx,x,st,time,i_elm,fields,x_min,x_max,n_real_param,real_param,n_int_param,int_param, &
  pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness)
end function pdf_sqrtPsiN_p_HL2A_parametric

function pdf_RZ_p_HL2A_template(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param, &
A, xi, omega, alpha, flatness) result(pdf)
  !! Same functional form as Wang2020 but with the radial profile f_psi(psi) to maximally excite the targeted TAE mode.
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  real*8,intent(in)                           :: A, xi, omega, alpha, flatness ! parameters for skew Gaussian f(sqrt(psi_N))=f(R,Z)
  !> Outputs:
  real*8 :: pdf
  !> Variables: 
  real *8 :: psi, psi_norm, s, p, f_psi, f_MJ, T_eV, mass_AMU, theta, logfmax
  real *8 :: E(3), B(3), electric_potential

  ! Here x = [R,Z,phi,p,u_beta]
  ! Designed skew Gaussian f(s) with sup f=1 and maximum negative gradient at sqrt(psi_N) ~ 0.54 to maximally excite the n=4, m=6,7 TAE mode in HL-2A.
  
  call fields%calc_EBpsiU(time,i_elm,st,x(3),E,B,psi,electric_potential)
  
  psi_norm = (psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis)
  s = sqrt(psi_norm)

  f_psi = skewed_normal(s, A, xi, omega, alpha, flatness)

  ! Compute the Maxwell-Juttner distribution function f_MJ(p) = p^2 * exp(-gamma/theta) 
  ! here p has unit [AMU·m/s]
  theta = real_param(1) ! kT/(m*c^2)
  T_eV = real_param(2) ! in [eV] 
  mass_AMU = real_param(3) ! in [AMU]
  logfmax = real_param(4)

  p = x(4) ! in [AMU·m/s]
  f_MJ = exp(log_f_MJ(p, theta, mass_AMU) - logfmax)

  if (f_MJ .lt. 0d0) write (*,*) "Error: f_MJ < 0 ", f_MJ, " p=", p, " theta=", theta, " mass_AMU=", mass_AMU

  pdf = f_psi * f_MJ
end function pdf_RZ_p_HL2A_template


function pdf_psiN_p_Relativistic_MJ(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  !! Parametric skewed-normal PDF using namelist parameters (pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness).
  !! Allows scanning s_max without recompilation.
  use mod_fields, only: fields_base
  implicit none
  integer,intent(in)                          :: nx,i_elm,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  !> Variables: 
  real *8 :: psi, psi_norm, p, f_psi, f_MJ, mass_AMU, theta, logfmax
  real *8 :: E(3), B(3), electric_potential

  ! Here x = [psiN,theta,phi,p,u_beta]
  psi_norm = x(1)
  if (psi_norm < 0d0 .or. psi_norm > 1d0) then
    pdf = 0d0
    return
  end if
  
  f_psi = skewed_normal(psi_norm, pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness)

  ! Compute the Maxwell-Juttner distribution function f_MJ(p) = p^2 * exp(-gamma/theta) 
  ! here p has unit [AMU·m/s]
  theta = real_param(1) ! kT/(m*c^2)
  mass_AMU = real_param(3) ! in [AMU]
  logfmax = real_param(4)

  p = x(4) ! in [AMU·m/s]
  f_MJ = exp(log_f_MJ(p, theta, mass_AMU) - logfmax)

  if (f_MJ .lt. 0d0) then
    write (*,*) "Error: f_MJ < 0 ", f_MJ, " p=", p, " theta=", theta, " mass_AMU=", mass_AMU
    pdf = 0d0
    return
  end if

  pdf = f_psi * f_MJ
end function pdf_psiN_p_Relativistic_MJ


function pdf_psiN_E_parametric(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  !! Parametric skewed-normal PDF using namelist parameters (pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness).
  !! Allows scanning s_max without recompilation.
  use mod_fields, only: fields_base
  implicit none
  integer,intent(in)                          :: nx,i_elm,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  !> Variables: 
  real *8 ::  psi_norm, f_psi, f_MB, T_eV, sup_fMB, Estar, E_eV

  ! Here x = [psiN,theta,phi,E,u_beta]  
  psi_norm = x(1)
  if (psi_norm < 0d0 .or. psi_norm > 1d0) then
    pdf = 0d0
    return
  end if
  E_eV = x(4) ! eV

  f_psi = skewed_normal(psi_norm, pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness)

  ! Compute the Maxwell-Boltzmann energy distribution f_MB(E) \propto sqrt(E) * exp(-E/kT)
  T_eV = real_param(2) ! in [eV] 
  Estar = T_eV / 2
  sup_fMB = sqrt(Estar) * exp(- 0.5d0) 
  f_MB =  1/sup_fMB * sqrt(E_eV) * exp(-E_eV/T_eV)

  if (f_MB .lt. 0d0) then
    write (*,*) "Error: f_MB < 0 ", f_MB, " E [eV]=", E_eV, " psiN=", psi_norm
    pdf = 0d0
    return
  end if

  pdf = f_psi * f_MB
end function pdf_psiN_E_parametric

function pdf_RZ_p_relativistic_MJ(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  !! Parametric skewed-normal PDF using namelist parameters (pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness).
  !! Allows scanning s_max without recompilation.
  use mod_fields, only: fields_base
  implicit none
  integer,intent(in)                          :: nx,i_elm,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  !> Variables: 
  real *8 ::  psi_norm, f_psi, f_MJ, T_eV, sup_fMB, Estar, E_eV, E(3), B(3), psi, electric_potential
  real *8 ::  theta, mass_AMU, logfmax, p

  ! Here x = [R,Z,phi,p,u_beta] and assuming i_elm is valid
  call fields%calc_EBpsiU(time,i_elm,st,x(3),E,B,psi,electric_potential)  
  psi_norm = (psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis)
  if (psi_norm < 0.01d0 .or. psi_norm > 0.99d0) then
    pdf = 0d0
    return
  end if

  f_psi = skewed_normal(psi_norm, pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness)

  ! Compute the Maxwell-Juttner distribution function f_MJ(p) = p^2 * exp(-gamma/theta) 
  ! here p has unit [AMU·m/s]
  theta = real_param(1) ! kT/(m*c^2)
  mass_AMU = real_param(3) ! in [AMU]
  logfmax = real_param(4)

  p = x(4) ! in [AMU·m/s]
  f_MJ = exp(log_f_MJ(p, theta, mass_AMU) - logfmax)

  if (f_MJ .lt. 0d0) then
    write (*,*) "Error: f_MJ < 0 ", f_MJ, " p=", p, " theta=", theta, " mass_AMU=", mass_AMU
    pdf = 0d0
    return
  end if

  pdf = f_psi * f_MJ
end function pdf_RZ_p_relativistic_MJ

function pdf_RZ_p_HL2A_flat_decay(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  !! Same functional form as Wang2020 but with the radial profile f_psi(psi) to maximally excite the targeted TAE mode.
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  !> Variables: 
  real *8 :: psi, psi_norm, s, p, f_psi, f_MJ, T_eV, mass_AMU, theta, logfmax
  real *8 :: width, s0 ! parameters for skew Gaussian f(sqrt(psi_N))=f(R,Z)
  real *8 :: E(3), B(3), electric_potential

  ! Here x = [R,Z,phi,p,u_beta]
  ! Designed skew Gaussian f(s) with sup f=1 and maximum negative gradient at sqrt(psi_N) ~ 0.54 to maximally excite the n=4, m=6,7 TAE mode in HL-2A.
  
  call fields%calc_EBpsiU(time,i_elm,st,x(3),E,B,psi,electric_potential)
  
  psi_norm = max((psi - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)
  s = sqrt(psi_norm)
  s0 = 0.54
  width = 0.020833
  f_psi = 0.5 * (1.0 - tanh((s - s0) / (2.0 * width)))

  ! Compute the Maxwell-Juttner distribution function f_MJ(p) = p^2 * exp(-gamma/theta) 
  ! here p has unit [AMU·m/s]
  theta = real_param(1) ! kT/(m*c^2)
  T_eV = real_param(2) ! in [eV] 
  mass_AMU = real_param(3) ! in [AMU]
  logfmax = real_param(4)

  p = x(4) ! in [AMU·m/s]
  f_MJ = exp(log_f_MJ(p, theta, mass_AMU) - logfmax)

  if (f_MJ .lt. 0d0) write (*,*) "Error: f_MJ < 0 ", f_MJ, " p=", p, " theta=", theta, " mass_AMU=", mass_AMU

  pdf = f_psi * f_MJ
end function pdf_RZ_p_HL2A_flat_decay

!> probability density function f(psi, p) = f_ITPA(psi) * f_MJ(p)/sup(f_MJ)
!> inputs:
!>   n_x:          (integer) number of variables:
!>   x:            (real8)(n_x) particle phase space coordinates 
!>                  1) poloidal flux 2) poloidal angle 3) toroidal angle 4) relativistic momentum
!>   st:           (real8)(2) local particle coordinates 
!>   time:         (real8) current simulation time
!>   i_elm:        (integer) jorek element
!>   fields:       (fields_base) JOREK MHD field class
!>   x_min:        (real8)(n_x) lower bound uniform sampling, 1) poloidal flux
!>                 2) poloidal angle, 3) toroidal angle, 4) relativistic momentum
!>   x_max:        (real8)(n_x) upper bound uniform sampling, 1) poloidal flux
!>                 2) poloidal angle, 3) toroidal angle, 4) relativistic momentum
!>   n_real_param: (integer) size of the real parameters
!>   real_param:   (real8)(2*n_elements+3) real gdf_sampler parameters:
!>                 1:              <- theta = kB*T/(m*c^2)
!>                 2:              <- temperature in eV
!>                 3:              <- particle mass in AMU
!>   n_int_param:  (integer) size of the integer parameters
!>   int_param:    (integer)(n_int_param) integer parameters:
!> outputs:
!>   pdf: (real8) value of the sampler probability density
function pdf_psi_p_itpa(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  !> Variables: 
  real *8 :: psi_norm, s, p, coeff(0:3), f_ITPA, f_MJ, T_eV, mass_AMU, gamma, theta, logfmax

  ! The ITPA distribution function f_ITPA(psi) = n_ITPA/n0 gives the spatial profile
  ! central densiy n0 should be 1.44131x10^17
  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  psi_norm = max((x(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)
  s = 0.957 * psi_norm + 0.043 * psi_norm**2 
  f_ITPA = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

  ! Compute the Maxwell-Juttner distribution function f_MJ(p) = p^2 * exp(-gamma/theta) 
  ! here p has unit [AMU·m/s]
  theta = real_param(1) ! kT/(m*c^2)
  T_eV = real_param(2) ! in [eV] 
  mass_AMU = real_param(3) ! in [AMU]
  logfmax = real_param(4)

  p = x(4) ! in [AMU·m/s]
  f_MJ = exp(log_f_MJ(p, theta, mass_AMU) - logfmax)

  if (f_MJ .lt. 0d0) write (*,*) "Error: f_MJ < 0 ", f_MJ, " p=", p, " theta=", theta, " mass_AMU=", mass_AMU

  pdf = f_ITPA * f_MJ
end function pdf_psi_p_itpa

function pdf_r_p_wang2020(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  !> Variables: 
  real *8 :: psi_norm, r, s, p, coeff(0:3), f_psi, f_MJ, T_eV, mass_AMU, theta, logfmax, R_axis, Z_axis, amin

  ! Here x = [R,Z,phi,p,u_beta]
  R_axis = real_param(5)
  Z_axis = real_param(6)
  amin   = real_param(7)
  r = sqrt((x(1)-R_axis)**2 + (x(2)-Z_axis)**2)
  ! Wang uses beta(r/a) = C*n(r/a) = 0.025 * np.exp(-(r_over_amin / 0.4) ** 2)
  ! Here assume sqrt(psi_norm) ~ r/a to have a similar shape
  if(r > amin) then
    f_psi = 0.d0
  else
    f_psi = exp(-(r/amin / 0.4) ** 2)
  endif
  ! Compute the Maxwell-Juttner distribution function f_MJ(p) = p^2 * exp(-gamma/theta) 
  ! here p has unit [AMU·m/s]
  theta = real_param(1) ! kT/(m*c^2)
  T_eV = real_param(2) ! in [eV] 
  mass_AMU = real_param(3) ! in [AMU]
  logfmax = real_param(4)

  p = x(4) ! in [AMU·m/s]
  f_MJ = exp(log_f_MJ(p, theta, mass_AMU) - logfmax)

  if (f_MJ .lt. 0d0) write (*,*) "Error: f_MJ < 0 ", f_MJ, " p=", p, " theta=", theta, " mass_AMU=", mass_AMU

  pdf = f_psi * f_MJ
end function pdf_r_p_wang2020

!> Wang2020 spatial profile expressed in flux coordinates: f(psiN, p) = f_psi(psiN) * f_MJ(p)/sup(f_MJ)
!> with r/a ~ sqrt(psiN), so f_psi(psiN) = exp(-(sqrt(psiN)/0.4)^2) = exp(-psiN/0.16).
!> inputs:
!>   x:            (real8)(nx) particle phase space coordinates
!>                 1) psiN 2) poloidal angle 3) toroidal angle 4) relativistic momentum 5) u_beta
!>   real_param:   (real8)(n_real_param) pdf parameters:
!>                 1: <- theta = kB*T/(m*c^2)
!>                 2: <- temperature in eV
!>                 3: <- particle mass in AMU
!>                 4: <- logfmax for the f_MJ/f_max rejection sampling
function pdf_psiN_p_wang2020(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  !> Variables:
  real *8 :: psi_norm, p, f_psi, f_MJ, mass_AMU, theta, logfmax

  ! Here x = [psiN,theta,phi,p,u_beta]
  psi_norm = x(1)
  if (psi_norm < 0d0 .or. psi_norm > 1d0) then
    pdf = 0d0
    return
  end if

  ! r/a ~ sqrt(psiN): f(r/a) = exp(-(r/a/0.4)^2) -> f(psiN) = exp(-psiN/0.4^2)
  f_psi = exp(-psi_norm / 0.4d0**2)

  ! Compute the Maxwell-Juttner distribution function f_MJ(p) = p^2 * exp(-gamma/theta)
  ! here p has unit [AMU·m/s]
  theta = real_param(1) ! kT/(m*c^2)
  mass_AMU = real_param(3) ! in [AMU]
  logfmax = real_param(4)

  p = x(4) ! in [AMU·m/s]
  f_MJ = exp(log_f_MJ(p, theta, mass_AMU) - logfmax)

  if (f_MJ .lt. 0d0) write (*,*) "Error: f_MJ < 0 ", f_MJ, " p=", p, " theta=", theta, " mass_AMU=", mass_AMU

  pdf = f_psi * f_MJ
end function pdf_psiN_p_wang2020


subroutine sqrtPsiN_Theta_Phi_P_to_gc_relativistic(p_inout,n_x,x,time,fields,n_real_param, real_param,n_int_param,int_param)
  use mod_particle_types, only: particle_base
  use mod_fields,         only: fields_base
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x
  integer,intent(in)               :: n_int_param,n_real_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                :: time
  real*8,dimension(n_x),intent(in) :: x
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)      :: fields
  !> Variables:
  real*8 :: psi, theta, phi, p, u_beta, beta, R, Z, p_perp, u_ppar_sign, E(3), B(3), electric_potential, normB, mass_AMU
  !> inputs-outputs:
  class(particle_base),intent(inout) :: p_inout

  psi = x(1)
  theta = x(2)
  phi = x(3)
  p = x(4) ! in [AMU·m/s]
  u_beta = 2*mod(x(5), 0.5d0) ! Use the same uniform(0,1) for computing the beta(1,1/2) and the sign for p_par.
  u_ppar_sign = x(5) - 0.5d0

  beta = 2.d0*u_beta - u_beta**2
  mass_AMU = real_param(2*fields%element_list%n_elements+3)

  ! find s,t, i_elm, and (R,Z) from (psi, theta, phi)
  call find_theta_psi(fields%node_list,fields%element_list,& 
  [real_param(1:fields%element_list%n_elements), real_param(fields%element_list%n_elements+1:2*fields%element_list%n_elements)],&
  theta, psi, phi,&
  real_param(2*fields%element_list%n_elements+1), real_param(2*fields%element_list%n_elements+2),&
  p_inout%i_elm, p_inout%st(1), p_inout%st(2), R, Z)

  p_inout%x = [R, Z, phi] ! set all three coordinates (R, Z, phi) - phi is needed for correct B-field in mu calculation


  select type( particle => p_inout)
    type is (particle_gc_relativistic)
      ! compute the parallel momentum [AMU·m/s]
      particle%p(1) = sign(p*sqrt(1.d0 - beta), u_ppar_sign)
    
      ! compute the magnetic moment [AMU·m^2/(T·s^2)]
      p_perp = p * sqrt(beta)
      
      call fields%calc_EBpsiU(time,particle%i_elm,particle%st,particle%x(3), E,B,psi,electric_potential); normB = norm2(B);
      particle%p(2) = p_perp**2 / (2.d0 * mass_AMU * normB)
      if(allocated(int_param)) then 
        particle%q    = int_param(1)
      else 
        write(*,*) "Error: int_param not allocated in sample_to_gc_relativistic for particle charge assignment"
      endif 
    class default
      write(*,*) "Error: particle type not supported in sample_to_gc_relativistic"
      stop
  end select

end subroutine

function particle_weight_one(nx,x,st,time,i_elm,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param) result(weight_res)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: nx,i_elm,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(2),intent(in)              :: st
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: weight_res
  ! Do not give weight for i_elm <= 0
  if (i_elm .le. 0 .or. i_elm .gt. fields%element_list%n_elements) then
    weight_res = 0.d0 
    return
  else
    weight_res = 1d0
  endif
end function particle_weight_one

function skewed_normal(s, A, xi, omega, alpha, flatness) result(f)
  real*8, intent(in) :: s, A, xi, omega, alpha, flatness
  real*8 :: f, z
  z = (s - xi) / omega
  f = A * exp(-0.5d0 * abs(z)**flatness) * (1.0d0 + erf(alpha * z / sqrt(2.0d0)))
end function skewed_normal

!> Uniform sampler for (R, Z, phi, ppar, mu) phase space.
!> Maps all 5 uniform random numbers linearly to [x_min, x_max],
!> then calls find_RZ to locate the JOREK element.
subroutine gdf_RZPhiPparMu_sampler(n_x,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param)
  use mod_fields, only: fields_base
  implicit none
  integer,intent(in) :: n_x,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(n_x),intent(in)            :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)               :: fields
  integer,intent(inout)               :: i_elm
  real*8,dimension(2),intent(inout)   :: st
  real*8,dimension(n_x),intent(inout) :: x
  integer :: ifail
  real*8  :: DUMMY_R, DUMMY_Z

  x(1:n_x) = x_min(1:n_x) + (x_max(1:n_x) - x_min(1:n_x)) * x(1:n_x)
  if (i_elm .le. 0 .or. i_elm .gt. fields%element_list%n_elements) i_elm = 1

  call find_RZ(fields%node_list, fields%element_list, x(1), x(2), DUMMY_R, DUMMY_Z, i_elm, st(1), st(2), ifail)
  if (ifail .ne. 0) then
    i_elm = 0
    st    = -1.d0
  endif
end subroutine gdf_RZPhiPparMu_sampler

subroutine gdf_PsiNThetaPhi_EorP_sampler(n_x,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param)
  use mod_fields, only: fields_base
  implicit none
  integer,intent(in) :: n_x,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(n_x),intent(in)            :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)               :: fields
  integer,intent(inout)               :: i_elm
  real*8,dimension(2),intent(inout)   :: st
  real*8,dimension(n_x),intent(inout) :: x
  !> Locals
  real*8 :: psi, R, Z

  x(1:n_x) = x_min(1:n_x) + (x_max(1:n_x) - x_min(1:n_x)) * x(1:n_x)
  ! Validate (psiN, theta) via find_theta_psi so only valid pairs pass
  ! rejection. Failed pairs set i_elm=0 and are rejected immediately.
  psi = x(1) * (ES%Psi_bnd - ES%Psi_axis) + ES%Psi_axis
  call find_theta_psi(fields%node_list, fields%element_list,&
  [real_param(1:fields%element_list%n_elements),&
   real_param(fields%element_list%n_elements+1:2*fields%element_list%n_elements)],&
  x(2), psi, x(3),&
  real_param(2*fields%element_list%n_elements+1),&
  real_param(2*fields%element_list%n_elements+2),&
  i_elm, st(1), st(2), R, Z)
  ! On failure, set sentinel st values so the rejection handler in the
  ! initialiser can reliably reject the particle (consistent with the
  ! error handling in gdf_RZPhiP_sampler_cb).
  if (i_elm .le. 0) then
    i_elm = 0
    st    = -1.d0
  endif
end subroutine gdf_PsiNThetaPhi_EorP_sampler

!> Equilibrium PDF: f = f_MJ(p) * g(P_phi)
!>   f_MJ = p^2 * exp(-gamma/theta)  (normalized to sup=1 via logfmax)
!>   g(P_phi) = exp(-(P_phi - Pphi0)^2 / (2*sigma^2))
!>   P_phi = R * ppar_SI * Bphi/B + q*e * psi
!> Pphi0 and sigma are hardcoded from the d13 sqrt(psi_N) profile.
!> x = [R, Z, phi, ppar, mu]
function pdf_E_Pphi(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  use mod_fields, only: fields_base
  use mod_gc_relativistic, only: compute_relativistic_factor
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  !> Variables:
  real*8 :: ppar, mu, Bnorm, Bphi, psi, p_total, theta, mass_AMU, logfmax
  real*8 :: f_MJ_norm, P_phi, g_Pphi, Rpos, q_charge
  real*8 :: E(3), B(3), electric_potential
  real*8 :: target_s, sigma_s, psi_target, Pphi0, sigma_Pphi

  ! x = [R, Z, phi, ppar, mu]
  Rpos = x(1)
  ppar  = x(4)  ! in [AMU·m/s]
  mu    = x(5)  ! pperp / 2mB in [AMU·m²/(T·s²)] 

  theta    = real_param(1)   ! kT/(m*c^2)
  ! real_param(2) = T_eV (unused here)
  mass_AMU = real_param(3)   ! in [AMU]
  logfmax  = real_param(4)   ! log of sup(f_MJ)
  q_charge = dble(int_param(1))  ! charge in units of e

  call fields%calc_EBpsiU(time,i_elm,st,x(3),E,B,psi,electric_potential)
  Bnorm = norm2(B)
  Bphi  = B(3)

  ! Total momentum and relativistic MJ factor
  p_total = sqrt(ppar**2 + 2.d0 * mass_AMU * mu * Bnorm) 
  f_MJ_norm = exp(log_f_MJ(p_total, theta, mass_AMU) - logfmax)

  ! Toroidal canonical angular momentum in SI [kg·m²/s]
  P_phi = Rpos * ppar * ATOMIC_MASS_UNIT * Bphi / Bnorm  + dble(q_charge) * EL_CHG * psi

  ! -- CHANGE THIS FOR DIFFERENT GAUSSIAN SHAPE of f(Pphi) --
  ! Gaussian g(P_phi) — designed so max |df/ds| occurs at target_s 
  !   P_phi ≈ qe * psi(s),  psi(s) = psi_axis + s^2 * (psi_bnd - psi_axis)
  !   Gaussian in P_phi centered at s0=0.35, width sigma_s=0.1626
  !   gives f''(s)=0 at s=0.5 (max gradient location)
  target_s = 0.470000d0
  sigma_s  = 0.045000d0
  ! ---------------------------------------------------------

  psi_target = ES%Psi_axis + target_s**2 * (ES%Psi_bnd - ES%Psi_axis)
  Pphi0      = dble(q_charge) * EL_CHG * psi_target ! q_charge is in units of e
  sigma_Pphi = abs(dble(q_charge) * EL_CHG * 2.0d0 * target_s * sigma_s * (ES%Psi_bnd - ES%Psi_axis))

  g_Pphi = exp(-(P_phi - Pphi0)**2 / (2.0d0 * sigma_Pphi**2))

  pdf = f_MJ_norm * g_Pphi

  if (pdf .lt. 0.d0) write (*,*) "Error: pdf_E_Pphi < 0 ", pdf
end function pdf_E_Pphi

!> Converter: directly copies ppar -> p(1), mu -> p(2).
!> No unit conversion needed — input is already in [AMU·m/s] and [AMU·m²/(T·s²)].
subroutine RZPhiPparMu_to_gc_relativistic_cb(p_inout,n_x,x,time,fields,n_real_param,real_param,n_int_param,int_param)
  use mod_particle_types, only: particle_base
  use mod_fields,         only: fields_base
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x
  integer,intent(in)               :: n_int_param,n_real_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                :: time
  real*8,dimension(n_x),intent(in) :: x
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)      :: fields
  !> Variables:
  real*8 :: R, Z, phi, ppar, mu, DUMMY_R, DUMMY_Z
  integer :: ifail
  !> inputs-outputs:
  class(particle_base),intent(inout) :: p_inout

  R = x(1)
  Z = x(2)
  phi = x(3)
  ppar = x(4)  ! already in [AMU·m/s]
  mu   = x(5)  ! already in [AMU·m²/(T·s²)]

  call find_RZ(fields%node_list, fields%element_list, R, Z, DUMMY_R, DUMMY_Z, p_inout%i_elm, p_inout%st(1), p_inout%st(2), ifail)
  p_inout%x = [R, Z, phi]

  select type( particle => p_inout)
    type is (particle_gc_relativistic)
      particle%p(1) = ppar
      particle%p(2) = mu
      if(allocated(int_param)) then
        particle%q    = int_param(1)
      else
        write(*,*) "Error: int_param not allocated in RZPhiPparMu_to_gc_relativistic_cb"
      endif
    class default
      write(*,*) "Error: particle type not supported in RZPhiPparMu_to_gc_relativistic_cb"
      stop
  end select
end subroutine RZPhiPparMu_to_gc_relativistic_cb

subroutine PsiNThetaPhiE_to_gc_relativistic_cb(p_inout,n_x,x,time,fields,n_real_param,real_param,n_int_param,int_param)
  use mod_particle_types, only: particle_base
  use mod_fields,         only: fields_base
  use mod_interp,        only: interp_RZ
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x
  integer,intent(in)               :: n_int_param,n_real_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                :: time
  real*8,dimension(n_x),intent(in) :: x
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)      :: fields
  !> Variables:
  real*8 :: R, Z, phi, p, ppar, mu, u_beta, u_ppar_sign, beta, p_par, p_perp, normB, KE_eV
  real*8 :: E(3), B(3), electric_potential, mass_AMU, gamma, psi, Lambda, B0
  real*8 :: Lambda_max, erf_a, erf_b
  !> inputs-outputs:
  class(particle_base),intent(inout) :: p_inout

  mass_AMU = real_param(2*fields%element_list%n_elements+3)
  B0       = abs(F0/R_geo)

  ! GDF sampler already validated i_elm>0 and st via find_theta_psi.
  ! Compute R,Z from known (st, i_elm) — fast, no element search.
  call interp_RZ(fields%node_list, fields%element_list, p_inout%i_elm, &
                  p_inout%st(1), p_inout%st(2), R, Z)
  p_inout%x = [R, Z, x(3)]

  call fields%calc_EBpsiU(time,p_inout%i_elm,p_inout%st,x(3),E,B,psi,electric_potential)
  normB = norm2(B)
  
  KE_eV = x(4) ! in [eV]
  gamma = 1.d0 + KE_eV * EL_CHG / (mass_AMU * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2)
  p     = mass_AMU * SPEED_OF_LIGHT * sqrt(max(gamma**2 - 1.d0, 0.d0))
  u_ppar_sign = x(5) - 0.5d0
  u_beta      = 2*mod(x(5), 0.5d0) ! Use the same uniform(0,1) for computing the beta(1,1/2) and the sign for p_par.
  
  if(Lambda_peak > LAMBDA_PEAK_MIN) then
    ! Maximum allowed Lambda at this location: ensures μB < 0.99·E
    Lambda_max = 0.9999d0 * B0 / normB
    if (Lambda_max <= 0.d0) then
      p_inout%i_elm = 0
      p_inout%st = -1.d0
      return
    end if
    ! Sample Λ from truncated Normal(Λ_peak, (ΔΛ/√2)²) with support [0, Λ_max]
    ! via inverse transform of the truncated Normal CDF:
    ! Λ = μ + σ√2 · erfinv( erf((a-μ)/(σ√2)) + u·(erf((b-μ)/(σ√2)) - erf((a-μ)/(σ√2))) )
    ! where μ=Λ_peak, σ=ΔΛ/√2, a=0, b=Λ_max, σ√2=ΔΛ
    erf_a = erf(-Lambda_peak / delta_Lambda)
    erf_b = erf((Lambda_max - Lambda_peak) / delta_Lambda)
    Lambda = Lambda_peak + delta_Lambda * c_erfinv(erf_a + u_beta * (erf_b - erf_a))
    ! Λ is guaranteed to be in [0, Λ_max] — no rejection needed
    mu = Lambda * (KE_eV * EL_CHG / ATOMIC_MASS_UNIT) / B0 ! in [AMU·m^2/(T·s^2)]
    ! Relativistic p_par: p_par^2 = p_total^2 - p_perp^2 = p^2 - 2*mass_AMU*mu*normB
    p_par = sign(sqrt(max(p**2 - 2.d0 * mass_AMU * mu * normB, 0.d0)), u_ppar_sign)
    if(abs(p_par) <= 1.d-8) write(*,*) "Warning: |p_par| <= 1.d-8, KE [eV] = ", KE_eV, " mu B [eV] = ", (mu * normB * ATOMIC_MASS_UNIT / EL_CHG), " p_par [AMU·m/s] = ", p_par
  else
    beta = 2.d0*u_beta - u_beta**2
    p_par = sign(p*sqrt(1.d0 - beta), u_ppar_sign)
    ! compute the magnetic moment [AMU·m^2/(T·s^2)]
    p_perp = p * sqrt(beta)
    mu = p_perp**2 / (2.d0 * mass_AMU * normB)  
  endif
  
  select type( particle => p_inout)
    type is (particle_gc_relativistic)
      particle%p(1) = p_par
      particle%p(2) = mu
      if(allocated(int_param)) then
        particle%q    = int_param(1)
      else
        write(*,*) "Error: int_param not allocated in RZPhiPparMu_to_gc_relativistic_cb"
      endif
    class default
      write(*,*) "Error: particle type not supported in RZPhiPparMu_to_gc_relativistic_cb"
      stop
  end select
end subroutine PsiNThetaPhiE_to_gc_relativistic_cb

subroutine PsiNThetaPhiP_to_gc_relativistic_cb(p_inout,n_x,x,time,fields,n_real_param,real_param,n_int_param,int_param)
  use mod_particle_types, only: particle_base
  use mod_fields,         only: fields_base
  use mod_interp,        only: interp_RZ
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x
  integer,intent(in)               :: n_int_param,n_real_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                :: time
  real*8,dimension(n_x),intent(in) :: x
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)      :: fields
  !> Variables:
  real*8 :: R, Z, phi, p, mu, u_beta, u_ppar_sign, beta, p_par, p_perp, normB, KE_eV
  real*8 :: E(3), B(3), electric_potential, mass_AMU, gamma, psi, Lambda, B0
  real*8 :: Lambda_max, erf_a, erf_b
  !> inputs-outputs:
  class(particle_base),intent(inout) :: p_inout

  mass_AMU = real_param(2*fields%element_list%n_elements+3)
  B0       = abs(F0/R_geo)

  ! GDF sampler already validated i_elm>0 and st via find_theta_psi.
  ! Compute R,Z from known (st, i_elm) — fast, no element search.
  call interp_RZ(fields%node_list, fields%element_list, p_inout%i_elm, &
                  p_inout%st(1), p_inout%st(2), R, Z)
  p_inout%x = [R, Z, x(3)]

  call fields%calc_EBpsiU(time,p_inout%i_elm,p_inout%st,x(3),E,B,psi,electric_potential)
  normB = norm2(B)
  
  p           = x(4)               ! in [AMU·m/s]
  gamma       = sqrt(1 + (p/(mass_AMU*SPEED_OF_LIGHT))**2) ! relativistic factor
  u_ppar_sign = x(5) - 0.5d0
  u_beta      = 2*mod(x(5), 0.5d0) ! Use the same uniform(0,1) for computing the beta(1,1/2) and the sign for p_par.
  
  if(Lambda_peak > LAMBDA_PEAK_MIN) then
    ! Maximum allowed Lambda at this location: ensures μB < E
    Lambda_max = 0.995d0 * B0 / normB
    if (Lambda_max <= 0.d0) then
      p_inout%i_elm = 0
      p_inout%st = -1.d0
      return
    end if
    KE_eV = (gamma - 1) * mass_AMU * (ATOMIC_MASS_UNIT / EL_CHG * SPEED_OF_LIGHT**2)  ! in [eV]
    ! Sample Λ from truncated Normal(Λ_peak, (ΔΛ/√2)²) with support [0, Λ_max]
    ! via inverse transform of the truncated Normal CDF:
    ! Λ = μ + σ√2 · erfinv( erf((a-μ)/(σ√2)) + u·(erf((b-μ)/(σ√2)) - erf((a-μ)/(σ√2))) )
    ! where μ=Λ_peak, σ=ΔΛ/√2, a=0, b=Λ_max, σ√2=ΔΛ
    erf_a = erf(-Lambda_peak / delta_Lambda)
    erf_b = erf((Lambda_max - Lambda_peak) / delta_Lambda)
    Lambda = Lambda_peak + delta_Lambda * c_erfinv(erf_a + u_beta * (erf_b - erf_a))
    ! Λ is guaranteed to be in [0, Λ_max] — no rejection needed
    mu = Lambda * (KE_eV * EL_CHG / ATOMIC_MASS_UNIT) / B0 ! in [AMU·m^2/(T·s^2)]
    ! Relativistic p_par: p_par^2 = p_total^2 - p_perp^2 = p^2 - 2*mass_AMU*mu*normB
    p_par = sign(sqrt(max(p**2 - 2.d0 * mass_AMU * mu * normB, 0.d0)), u_ppar_sign)
    if(abs(p_par) <= 1.d-8) write(*,*) "Warning: |p_par| <= 1.d-8, KE [eV] = ", KE_eV, " mu B [eV] = ", (mu * normB * ATOMIC_MASS_UNIT / EL_CHG), " p_par [AMU·m/s] = ", p_par
  else !isotropic momentum distribution
    beta = 2.d0*u_beta - u_beta**2
    p_par = sign(p*sqrt(1.d0 - beta), u_ppar_sign)
    ! compute the magnetic moment [AMU·m^2/(T·s^2)]
    p_perp = p * sqrt(beta)
    mu = p_perp**2 / (2.d0 * mass_AMU * normB)  
  endif
  

  select type( particle => p_inout)
    type is (particle_gc_relativistic)
      particle%p(1) = p_par
      particle%p(2) = mu
      if(allocated(int_param)) then
        particle%q    = int_param(1)
      else
        write(*,*) "Error: int_param not allocated in RZPhiPparMu_to_gc_relativistic_cb"
      endif
    class default
      write(*,*) "Error: particle type not supported in RZPhiPparMu_to_gc_relativistic_cb"
      stop
  end select
end subroutine PsiNThetaPhiP_to_gc_relativistic_cb

subroutine RZPhi_to_gc_relativistic_temp_gradient(p_inout,n_x,x,time,fields,n_real_param,real_param,n_int_param,int_param)
  use mod_particle_types, only: particle_base
  use mod_fields,         only: fields_base
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x
  integer,intent(in)               :: n_int_param,n_real_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                :: time
  real*8,dimension(n_x),intent(in) :: x
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)      :: fields
  !> Variables:
  real*8 :: R, Z, phi, ppar, mu, DUMMY_R, DUMMY_Z, f_val, gamma, mass_AMU
  real*8 :: u_beta, u_ppar_sign, beta, p_perp, normB
  integer :: ifail
  real *8 :: psi, psi_norm, p, f_psi, T_eV, kinetic_energy
  real *8 :: E(3), B(3), electric_potential
  !> inputs-outputs:
  class(particle_base),intent(inout) :: p_inout
  
  mass_AMU = real_param(2*fields%element_list%n_elements+3)
  R = x(1)
  Z = x(2)
  phi = x(3)
  p_inout%x = [R, Z, phi] ! set all three coordinates (R, Z, phi) - phi is needed for correct B-field in mu calculation

  ! find s,t, i_elm from (R, Z)
  call find_RZ(fields%node_list, fields%element_list, R, Z, DUMMY_R, DUMMY_Z, p_inout%i_elm, p_inout%st(1), p_inout%st(2), ifail)

  call fields%calc_EBpsiU(time,p_inout%i_elm,p_inout%st,phi,E,B,psi,electric_potential)
  normB = norm2(B)

  ! Compute the momentum magnitude from the particle's temperature 
  psi_norm = (psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis)
  f_psi = skewed_normal(psi_norm, pdf_A, pdf_xi, pdf_omega, pdf_alpha, pdf_flatness) ! sup(f_psi) = 1
  T_eV = f_psi * T_EP_eV ! in [eV]
  kinetic_energy = 3.0/2.0 * T_eV * EL_CHG ! in [J]
  gamma = 1 + kinetic_energy / (ATOMIC_MASS_UNIT * mass_AMU * SPEED_OF_LIGHT**2) ! gamma = 1 + E_kin/(m*c^2), mass in [AMU], E_kin in [eV], c in [m/s]
  p  = sqrt(gamma ** 2  - 1)* mass_AMU * SPEED_OF_LIGHT ! in [AMU·m/s]

  u_beta = 2*mod(x(5), 0.5d0) ! Use the same uniform(0,1) for computing the beta(1,1/2) and the sign for p_par.
  u_ppar_sign = x(5) - 0.5d0

  beta = 2.d0*u_beta - u_beta**2

  select type( particle => p_inout)
    type is (particle_gc_relativistic)
      ! compute the parallel momentum [AMU·m/s]
      particle%p(1) = sign(p*sqrt(1.d0 - beta), u_ppar_sign)
      ! compute the magnetic moment [AMU·m^2/(T·s^2)]
      p_perp = p * sqrt(beta)
      particle%p(2) = p_perp**2 / (2.d0 * mass_AMU * normB)
      if(allocated(int_param)) then 
        particle%q    = int_param(1)
      else 
        write(*,*) "Error: int_param not allocated in sample_to_gc_relativistic for particle charge assignment"
      endif 
    class default
      write(*,*) "Error: particle type not supported in sample_to_gc_relativistic"
      stop
  end select
end subroutine RZPhi_to_gc_relativistic_temp_gradient

!> Convert (psiN, theta, phi, p_total, u_beta) to a particle_kinetic_relativistic
!> Same as PsiNThetaPhiP_to_gc_relativistic_cb but constructs Cartesian momentum
!> with a random gyro-angle for the perpendicular component.
subroutine PsiNThetaPhiP_to_kinetic_relativistic_cb(p_inout,n_x,x,time,fields,n_real_param,real_param,n_int_param,int_param)
  use mod_particle_types, only: particle_base
  use mod_fields,         only: fields_base
  use mod_interp,        only: interp_RZ
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_pusher_tools, only: get_orthonormals
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x
  integer,intent(in)               :: n_int_param,n_real_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                :: time
  real*8,dimension(n_x),intent(in) :: x
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)      :: fields
  !> Variables:
  real*8 :: R, Z, phi, p, mu, beta, u_beta, u_ppar_sign, p_par, p_perp, normB, KE_eV
  real*8 :: E(3), B(3), electric_potential, mass_AMU, gamma, psi, Lambda, B0
  real*8 :: Lambda_max, erf_a, erf_b
  real*8 :: B_hat(3), e1(3), e2(3), gyro_angle, harvest(1)
  real*8 :: p_cyl(3)
  !> inputs-outputs:
  class(particle_base),intent(inout) :: p_inout

  mass_AMU = real_param(2*fields%element_list%n_elements+3)
  B0       = abs(F0/R_geo)

  ! GDF sampler already validated i_elm>0 and st via find_theta_psi.
  ! Compute R,Z from known (st, i_elm) — fast, no element search.
  call interp_RZ(fields%node_list, fields%element_list, p_inout%i_elm, &
                  p_inout%st(1), p_inout%st(2), R, Z)
  p_inout%x = [R, Z, x(3)]

  call fields%calc_EBpsiU(time,p_inout%i_elm,p_inout%st,x(3),E,B,psi,electric_potential)
  normB = norm2(B)
  
  p           = x(4)               ! in [AMU·m/s]
  gamma       = sqrt(1 + (p/(mass_AMU*SPEED_OF_LIGHT))**2) ! relativistic factor
  u_ppar_sign = x(5) - 0.5d0
  u_beta      = 2*mod(x(5), 0.5d0) ! Use the same uniform(0,1) for computing the beta and the sign for p_par.
  
  if(Lambda_peak > LAMBDA_PEAK_MIN) then
    ! Maximum allowed Lambda at this location: ensures μB < E
    Lambda_max = 0.995d0 * B0 / normB
    if (Lambda_max <= 0.d0) then
      p_inout%i_elm = 0
      p_inout%st = -1.d0
      return
    end if
    KE_eV = (gamma - 1) * mass_AMU * (ATOMIC_MASS_UNIT / EL_CHG * SPEED_OF_LIGHT**2)  ! in [eV]
    ! Sample Λ from truncated Normal(Λ_peak, (ΔΛ/√2)²) with support [0, Λ_max]
    erf_a = erf(-Lambda_peak / delta_Lambda)
    erf_b = erf((Lambda_max - Lambda_peak) / delta_Lambda)
    Lambda = Lambda_peak + delta_Lambda * c_erfinv(erf_a + u_beta * (erf_b - erf_a))
    ! Λ is guaranteed to be in [0, Λ_max] — no rejection needed
    mu = Lambda * (KE_eV * EL_CHG / ATOMIC_MASS_UNIT) / B0 ! in [AMU·m^2/(T·s^2)]
    ! Relativistic p_par: p_par^2 = p_total^2 - p_perp^2 = p^2 - 2*mass_AMU*mu*normB
    p_par = sign(sqrt(max(p**2 - 2.d0 * mass_AMU * mu * normB, 0.d0)), u_ppar_sign)
    if(abs(p_par) <= 1.d-8) write(*,*) "Warning: |p_par| <= 1.d-8, KE [eV] = ", KE_eV, " mu B [eV] = ", (mu * normB * ATOMIC_MASS_UNIT / EL_CHG), " p_par [AMU·m/s] = ", p_par
  else !isotropic momentum distribution
    beta = 2.d0*u_beta - u_beta**2
    p_par = sign(p*sqrt(1.d0 - beta), u_ppar_sign)
    ! compute the magnetic moment [AMU·m^2/(T·s^2)]
    p_perp = p * sqrt(beta)
    mu = p_perp**2 / (2.d0 * mass_AMU * normB)
  endif
  
  ! Construct orthonormal basis perpendicular to B
  if (normB > 0.d0) then
    B_hat = B / normB
    call get_orthonormals(B_hat, e1, e2)
    ! Random gyro-angle for perpendicular momentum direction
    call random_number(harvest)
    gyro_angle = TWOPI * harvest(1)
    ! Recompute p_perp from mu if using Lambda distribution (mu already set above)
    if (Lambda_peak > LAMBDA_PEAK_MIN) then
      p_perp = sqrt(max(2.d0 * mass_AMU * mu * normB, 0.d0))
    end if
    ! Build cylindrical momentum: p_cyl = p_par * B_hat + p_perp * (cos(α)*e1 + sin(α)*e2)
    p_cyl = p_par * B_hat + p_perp * (cos(gyro_angle) * e1 + sin(gyro_angle) * e2)
  else
    ! Degenerate case: B=0, put all momentum in z-direction
    p_cyl = [0.d0, 0.d0, p]
  end if

  select type( particle => p_inout)
    type is (particle_kinetic_relativistic)
      ! Convert cylindrical momentum to Cartesian
      particle%p = vector_cylindrical_to_cartesian(particle%x(3), p_cyl)
      if(allocated(int_param)) then
        particle%q    = int_param(1)
      else
        write(*,*) "Error: int_param not allocated in PsiNThetaPhiP_to_kinetic_relativistic_cb"
      endif
    class default
      write(*,*) "Error: particle type not supported in PsiNThetaPhiP_to_kinetic_relativistic_cb"
      stop
  end select
end subroutine PsiNThetaPhiP_to_kinetic_relativistic_cb

end module mod_tae_loop_init_callbacks