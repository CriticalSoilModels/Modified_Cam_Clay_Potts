! Functions for calculating and updateing the state parameters

module mod_state_params

   ! use statements here
   use kind_precision_module, only: dp, i32
	use mod_stress_invariants, only: calc_mean_stress

   implicit none

contains

   pure function calc_pore_water_pressure(K_water, porosity, dEps) result(pore_pressure)
      real(dp), intent(in) :: K_water, porosity, dEps(6)

      pore_pressure = K_water/porosity * sum(dEps(1:3))

   end function calc_pore_water_pressure

   pure function calc_cam_clay_bulk_modulus_2(void_ratio,  kappa, abs_mean_stress) result(K)
      real(kind = dp), intent(in) :: void_ratio, kappa, abs_mean_stress
      real(kind = dp) :: K

      ! not sure what this one is caculation
      K = (void_ratio + 1.0_dp)/kappa * abs_mean_stress
   end function calc_cam_clay_bulk_modulus_2

   pure function calc_cam_clay_bulk_modulus(Sig, v, kappa, mean_stress) result(K)
      real(kind = dp), intent(in) :: Sig(6), lambda, N, kappa

      
      K = v/kappa * mean_stress
   end function calc_cam_clay_bulk_modulus

	pure function calc_specific_volume(N, lambda, precon_stress)
		real(kind = dp), intent(in) :: N, lambda, precon_stress
		real(kind = dp) :: v

		! Calc the mean stress constrained so that it works for the log?
      ! Constrain the preconsolidation stress so that the log can evaluate it
		! TODO: Need to make sure this makes sense
		constr_precon_stress = max(precon_stress, 1.0_dp)

      v = N - lambda * log(constr_precon_stress)

	end function calc_specific_volume

end module mod_state_params
