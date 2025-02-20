module mod_elastic_stress
   use kind_precision_module, only: dp, i32

   implicit none

contains
   !> \brief Constructs the stiffness matrix DE based on shear modulus and Poisson's ratio.
   !!
   !! This function calculates and constructs a 6x6 stiffness matrix based on the given
   !! shear modulus and Poisson's ratio.
   !!
   !! \param[in] shear_modulus Shear modulus of the material.
   !! \param[in] poisson_ratio Poisson's ratio of the material.
   !! \return Stiffness matrix DE of size (6, 6).
   !!
   function construct_stiffness_matrix( shear_modulus, poisson_ratio ) result(stiff_matrix)
      ! Construct the stiffness matrix DE
      real(kind = dp), intent(in) :: shear_modulus, poisson_ratio
      real(kind = dp) :: stiff_matrix(6, 6)

      ! Local variables
      real(kind = dp) :: F1
      real(kind = dp) :: F2
      integer(kind = i32) :: i

      stiff_matrix(:, :) = 0.0_dp

      ! Calc the values needed for the stiffness matrix
      F1 = 2 * shear_modulus * ( 1 - poisson_ratio ) / ( 1- 2 * poisson_ratio )
      F2 = 2 * shear_modulus * poisson_ratio / ( 1 - 2 * poisson_ratio )

      !---- Fill the stiffness matrix ----

      ! Zero the matrix
      stiff_matrix = 0.0_dp

      ! Fill the upper block
      stiff_matrix(1:3, 1:3) = F2

      ! Loop over the first three diagonals
      do i = 1, 3
         stiff_matrix(i, i) = F1
      end do
      ! Loop over the 4th through 6th diagonal
      do i = 4, 6
         stiff_matrix(i, i) = shear_modulus
      end do
   end function construct_stiffness_matrix

   pure function calc_shear_modulus(nu, K) result(G)
      real(kind = dp), intent(in) :: nu, K
      real(kind = dp) :: G

      G = 3.0_dp * ( 1.0_dp - 2.0_dp*nu) / ( 2.0_dp * (1.0_dp+nu)) * K
   end function calc_shear_modulus

end module mod_elastic_stress
