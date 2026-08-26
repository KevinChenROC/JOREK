!> Module containing projection functions for phase-space diagnostics
!> These functions are used with mod_phase_space_project to project
!> particle distributions onto various grids
module mod_projection_functions_phase
  use mod_particle_sim
  use mod_particle_types
  use mod_import_experimental_dist, only: calculate_B
  use mod_math_operators, only: cross_product
  use constants, only: ATOMIC_MASS_UNIT, EL_CHG, SPEED_OF_LIGHT
  use phys_module, only: F0, R_geo
  use equil_info, only: ES
  
  implicit none
  private
  public :: proj_RZPE, proj_PparE, proj_PparMu, proj_rminor, proj_RZ, proj_psiN, proj_sqrtPsiN, proj_EPphi, proj_EPphiMu, proj_EMu
  public :: set_T_EP_eV
  
  !> EP temperature in eV - can be set via set_T_EP_eV subroutine
  real*8, save :: T_EP_eV = 1.0d5
  
contains

  !> Set the EP temperature for proj_PparMu normalization
  subroutine set_T_EP_eV(T_EP_eV_in)
    real*8, intent(in) :: T_EP_eV_in
    T_EP_eV = T_EP_eV_in
  end subroutine set_T_EP_eV


  
  !> Project R, Z, parallel normed velocity vpar/|v|, and energy [keV]
  pure function proj_RZPE(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    real*8                                   :: B_tmp(3)
    real*8                                   :: proj_RZPE(ndim)
    B_tmp = calculate_B(sim%fields, particle%i_elm, particle%st(1), particle%st(2), particle%x(3))
    select type(p => particle)
      type is (particle_kinetic_leapfrog)
        proj_RZPE = [particle%x(1), particle%x(2), &
                     dot_product(p%v, B_tmp)/norm2(B_tmp)/norm2(p%v), &
                     0.5d0*sim%groups(1)%mass*ATOMIC_MASS_UNIT*dot_product(p%v, p%v)/1000.d0/EL_CHG]
    end select  
  end function proj_RZPE

  !> Project parallel velocity (normalized) and energy (in keV)
  pure function proj_PparE(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    real*8                                   :: B_tmp(3)
    real*8                                   :: proj_PparE(ndim)
    B_tmp = calculate_B(sim%fields, particle%i_elm, particle%st(1), particle%st(2), particle%x(3))
    select type(p => particle)
      type is (particle_kinetic_leapfrog)
        ! Project parallel velocity (normalized) and energy (in keV)
        proj_PparE = [dot_product(p%v, B_tmp)/norm2(B_tmp)/norm2(p%v), &
                      0.5d0*sim%groups(1)%mass*ATOMIC_MASS_UNIT*dot_product(p%v, p%v)/1000.d0/EL_CHG]
    end select  
  end function proj_PparE

  !> Project normalized parallel momentum ppar/mc 
  !> and magnetic moment [keV/T]
  pure function proj_PparMu(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    !> Output
    real*8                                   :: proj_PparMu(ndim)
    !> Variables:
    real*8                                   :: B_vec(3), B0, B_mag, mass_AMU, mu_norm
    mass_AMU = sim%groups(group)%mass
    select type(p => particle)
      type is (particle_kinetic_leapfrog)
        ! For non-relativistic kinetic leapfrog with phase space coordinates (x, v_R, v_Z, v_phi)
        B_vec = calculate_B(sim%fields, particle%i_elm, particle%st(1), particle%st(2), particle%x(3))
        B_mag = norm2(B_vec)
        proj_PparMu = [dot_product(p%v, B_vec/B_mag)/SPEED_OF_LIGHT, &
                      0.5d0*mass_AMU*ATOMIC_MASS_UNIT*norm2(cross_product(p%v, B_vec/B_mag))**2/B_mag/1000.d0/EL_CHG]
      type is (particle_gc_relativistic)
        ! For relativistic gc particles with phase space coordinates (x, ppar, mu)
        proj_PparMu = [p%p(1)/(mass_AMU * SPEED_OF_LIGHT), &
                      p%p(2)*ATOMIC_MASS_UNIT/1000.d0/EL_CHG] ! p%p(2) is mu in [AMU m^2/(s^2·T)]
      type is (particle_gc_vpar)
        ! For nonrelativistic gc particles with phase space coordinates (x, E, mu)
        proj_PparMu = [p%vpar/SPEED_OF_LIGHT, &
                      p%mu*mass_AMU*ATOMIC_MASS_UNIT/1000.d0/EL_CHG] ! p%mu is 0.5 vperp^2/B in [J/(kg·T)]
      type is (particle_kinetic_relativistic)
        ! For relativistic kinetic (full-orbit) particles with Cartesian momentum (px, py, pz)
        ! Convert to cylindrical, compute p_par and mu from p and B
        block
          use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
          real*8 :: p_cyl(3), p_perp(3), B_hat(3), B_norm, p_par
          B_vec = calculate_B(sim%fields, particle%i_elm, particle%st(1), particle%st(2), particle%x(3))
          B_norm = norm2(B_vec)
          if (B_norm > 0.d0) then
            B_hat = B_vec / B_norm
            p_cyl = vector_cartesian_to_cylindrical(particle%x(3), p%p)
            p_par = dot_product(p_cyl, B_hat)
            p_perp = p_cyl - p_par * B_hat
            proj_PparMu = [p_par / (mass_AMU * SPEED_OF_LIGHT), &
                          (p_perp(1)**2 + p_perp(2)**2 + p_perp(3)**2) / (2.d0 * B_norm * mass_AMU) &
                          * ATOMIC_MASS_UNIT / EL_CHG / 1.d3]
          else
            proj_PparMu = [0.d0, 0.d0]
          end if
        end block
    end select  
  end function proj_PparMu

  !> Project minor radius [m]
  pure function proj_rminor(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    real*8                                   :: proj_rminor(ndim)
    select type(p => particle)
      type is (particle_kinetic_leapfrog)
        proj_rminor = [((particle%x(1)-R_geo)**2 + particle%x(2)**2)**0.5]
    end select  
  end function proj_rminor

  !> Project R [m] and Z [m] coordinates
  pure function proj_RZ(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    real*8                                   :: proj_RZ(ndim)
    proj_RZ = [particle%x(1), particle%x(2)]
  end function proj_RZ
  
  !> Project normalized poloidal flux psi_N in [0,1]
  function proj_psiN(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    real*8                                   :: proj_psiN(ndim)
    real*8                                   :: psi, E(3), B(3), U
    
    ! Compute psi at particle location
    call sim%fields%calc_EBpsiU(sim%time, particle%i_elm, particle%st, particle%x(3), E, B, psi, U)
    
    ! Normalize to [0,1]: psi_N = (psi - psi_axis) / (psi_bnd - psi_axis)
    proj_psiN = [max((psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis),0.d0)]
  end function proj_psiN

  function proj_sqrtPsiN(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    real*8                                   :: proj_sqrtPsiN(ndim)
    real*8                                   :: psi, E(3), B(3), U
    
    ! Compute psi at particle location
    call sim%fields%calc_EBpsiU(sim%time, particle%i_elm, particle%st, particle%x(3), E, B, psi, U)
    
    ! Normalize to [0,1]: psi_N = (psi - psi_axis) / (psi_bnd - psi_axis)
    proj_sqrtPsiN = [sqrt(max((psi - ES%Psi_axis) / (ES%Psi_bnd - ES%Psi_axis),0.d0))]
  end function proj_sqrtPsiN

  !> Project (E, P_phi) — the two invariants of motion for the GC equations.
  !>   dim 1: kinetic energy in keV
  !>   dim 2: canonical toroidal angular momentum in [AMU·m²/s]
  !>           Pphi = Pphi_SI / ATOMIC_MASS_UNIT
  function proj_EPphi(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    real*8                                   :: proj_EPphi(ndim)
    real*8                                   :: mass_AMU, B(3), Bnorm, Bphi, psi, U
    real*8                                   :: gamma, p_total, E_kin_keV, P_phi_SI, P_phi_amu
    real*8                                   :: E_dummy(3)
    integer                                  :: q_charge
    
    mass_AMU = sim%groups(group)%mass
    
    select type(p => particle)
    type is (particle_gc_relativistic)
      q_charge = p%q  ! charge in units of e
      
      ! Compute B at particle location
      call sim%fields%calc_EBpsiU(sim%time, p%i_elm, p%st, p%x(3), E_dummy, B, psi, U)
      Bnorm = norm2(B)
      Bphi  = B(3)
      
      ! Kinetic energy in keV
      p_total = sqrt(p%p(1)**2 + 2.d0 * mass_AMU * p%p(2) * Bnorm)
      gamma = sqrt(1.d0 + (p_total / (mass_AMU * SPEED_OF_LIGHT))**2)
      E_kin_keV = (gamma - 1.d0) * mass_AMU * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 &
                  / EL_CHG / 1.d3  ! in keV
      
      ! Canonical toroidal angular momentum in SI [kg·m²/s]
      ! P_phi_SI = R * ppar_SI * Bphi/B + q * e * psi
      P_phi_SI = p%x(1) * p%p(1) * ATOMIC_MASS_UNIT * Bphi / Bnorm &
               + dble(q_charge) * EL_CHG * psi
      
      ! Convert to AMU·m²/s
      P_phi_amu = P_phi_SI / ATOMIC_MASS_UNIT
      
      proj_EPphi = [E_kin_keV, P_phi_amu]
    class default
      proj_EPphi = [0.d0, 0.d0]
    end select
  end function proj_EPphi

  !> Project (E, P_phi, mu) — the 3 invariants of motion for the GC equations.
  !>   dim 1: kinetic energy in keV
  !>   dim 2: canonical toroidal angular momentum in [AMU·m²/s]
  !>           Pphi = Pphi_SI / ATOMIC_MASS_UNIT
  !>   dim 3: magnetic moment in [keV/T]
  function proj_EPphiMu(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    real*8                                   :: proj_EPphiMu(ndim)
    real*8                                   :: mass_AMU, B(3), Bnorm, Bphi, psi, U
    real*8                                   :: gamma, p_total, E_kin_keV, P_phi_SI, P_phi_amu, mu_keV_T
    real*8                                   :: E_dummy(3)
    integer                                  :: q_charge
    
    mass_AMU = sim%groups(group)%mass
    
    select type(p => particle)
    type is (particle_gc_relativistic)
      q_charge = p%q  ! charge in units of e
      
      ! Compute B at particle location
      call sim%fields%calc_EBpsiU(sim%time, p%i_elm, p%st, p%x(3), E_dummy, B, psi, U)
      Bnorm = norm2(B)
      Bphi  = B(3)
      
      ! Kinetic energy in keV
      p_total = sqrt(p%p(1)**2 + 2.d0 * mass_AMU * p%p(2) * Bnorm)
      gamma = sqrt(1.d0 + (p_total / (mass_AMU * SPEED_OF_LIGHT))**2)
      E_kin_keV = (gamma - 1.d0) * mass_AMU * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 &
                  / EL_CHG / 1.d3  ! in keV
      
      ! Canonical toroidal angular momentum in SI [kg·m²/s]
      ! P_phi_SI = R * ppar_SI * Bphi/B + q * e * psi
      P_phi_SI = p%x(1) * p%p(1) * ATOMIC_MASS_UNIT * Bphi / Bnorm &
               + dble(q_charge) * EL_CHG * psi
      
      ! Convert to AMU·m²/s
      P_phi_amu = P_phi_SI / ATOMIC_MASS_UNIT
      
      ! Magnetic moment in [keV/T]
      mu_keV_T = p%p(2) * ATOMIC_MASS_UNIT / EL_CHG / 1.d3
      
      proj_EPphiMu = [E_kin_keV, P_phi_amu, mu_keV_T]
    class default
      proj_EPphiMu = [0.d0, 0.d0, 0.d0]
    end select
  end function proj_EPphiMu

  !> Project (E,  mu) — the 2 invariants of motion for the GC equations.
  !>   dim 1: kinetic energy in keV
  !>   dim 2: magnetic moment in [keV/T]
  function proj_EMu(ndim, sim, group, particle)
    type(particle_sim),           intent(in) :: sim
    integer,                      intent(in) :: ndim, group
    class(particle_base),         intent(in) :: particle
    real*8                                   :: proj_EMu(ndim)
    real*8                                   :: mass_AMU, B(3), Bnorm, Bphi, psi, U
    real*8                                   :: gamma, p_total, E_kin_keV, P_phi_SI, P_phi_amu, mu_keV_T
    real*8                                   :: E_dummy(3)
    integer                                  :: q_charge
    
    mass_AMU = sim%groups(group)%mass
    
    select type(p => particle)
    type is (particle_gc_relativistic)
      q_charge = p%q  ! charge in units of e
      
      ! Compute B at particle location
      call sim%fields%calc_EBpsiU(sim%time, p%i_elm, p%st, p%x(3), E_dummy, B, psi, U)
      Bnorm = norm2(B)
      Bphi  = B(3)
      
      ! Kinetic energy in keV
      p_total = sqrt(p%p(1)**2 + 2.d0 * mass_AMU * p%p(2) * Bnorm)
      gamma = sqrt(1.d0 + (p_total / (mass_AMU * SPEED_OF_LIGHT))**2)
      E_kin_keV = (gamma - 1.d0) * mass_AMU * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 &
                  / EL_CHG / 1.d3  ! in keV
      
      ! Magnetic moment in [keV/T]
      mu_keV_T = p%p(2) * ATOMIC_MASS_UNIT / EL_CHG / 1.d3
      
      proj_EMu = [E_kin_keV, mu_keV_T]
    type is (particle_kinetic_relativistic)
      ! For relativistic kinetic (full-orbit) particles with Cartesian momentum
      q_charge = p%q
      call sim%fields%calc_EBpsiU(sim%time, p%i_elm, p%st, p%x(3), E_dummy, B, psi, U)
      Bnorm = norm2(B)
      ! Kinetic energy in keV from Cartesian momentum
      p_total = sqrt(p%p(1)**2 + p%p(2)**2 + p%p(3)**2)
      gamma = sqrt(1.d0 + (p_total / (mass_AMU * SPEED_OF_LIGHT))**2)
      E_kin_keV = (gamma - 1.d0) * mass_AMU * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 &
                  / EL_CHG / 1.d3
      ! Magnetic moment from p_perp^2 / (2*mass*B)  [keV/T]
      block
        use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
        real*8 :: p_cyl(3), p_perp(3), p_par, B_hat(3)
        if (Bnorm > 0.d0) then
          B_hat = B / Bnorm
          p_cyl = vector_cartesian_to_cylindrical(p%x(3), p%p)
          p_par = dot_product(p_cyl, B_hat)
          p_perp = p_cyl - p_par * B_hat
          mu_keV_T = (p_perp(1)**2 + p_perp(2)**2 + p_perp(3)**2) / (2.d0 * Bnorm * mass_AMU) &
                     * ATOMIC_MASS_UNIT / EL_CHG / 1.d3
        else
          mu_keV_T = 0.d0
        end if
      end block
      proj_EMu = [E_kin_keV, mu_keV_T]
    class default
      proj_EMu = [0.d0, 0.d0]
    end select
  end function proj_EMu

end module mod_projection_functions_phase
