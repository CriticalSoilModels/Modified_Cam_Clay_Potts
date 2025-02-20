

   !> \brief Calculates the increment of stress.
   !!
   !! This function calculates the increment of stress \f$ d\sigma \f$ given the stiffness
   !! matrix \f$ DE \f$ and the strain \f$ \epsilon \f$.
   !!
   !! \param[in] strain Strain vector \f$ \epsilon \f$ of size (6).
   !! \param[in] stiff_matrix Stiffness matrix \f$ DE \f$ of size (6, 6).
   !! \return Stress increment vector \f$ d\sigma \f$ of size (6).
   !!
   function calc_elastic_stress_increment(strain, stiff_matrix ) result(stress_inc)
      real(kind = dp), intent(in) :: strain(6), stiff_matrix(6,6)
      real(kind = dp) :: stress_inc(6)

      ! Calc the stress increment
      stress_inc = matmul(stiff_matrix, strain)

   end function
   
   subroutine PrincipalSig(IOpt, S, xN1, xN2, xN3, S1, S2, S3,&
      P, Q,J,theta)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions
      !            from cartesian stress vector
      !
      !  IOpt            I   I     flag to calculate principal direction (IOpt = 1)
      !  IntGlo          I   I     global ID of Gauss point or particle
      !  S               I   R()   cartesian stress
      !  xN1, xN2, xN3   O   R()   principal direction
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      !  J               O   R     sqrt J2
      !  theta           O   R     lode angle
      !
      !-------------------------------------------------------------------

      implicit none

      ! arguments
      integer, intent(in) :: IOpt
      real(kind = dp), intent(in) :: S(6)
      real(kind = dp), intent(out) :: xN1(3), xN2(3), xN3(3),&
         S1, S2, S3, P, Q, J, theta

      if (IOpt .eq. 1) then
         call Eig_3(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q,J,theta) ! Calculate principal direction
      else
         call Eig_3a(0,S,S1,S2,S3,P,Q) ! Do not calculate principal direction
      end if

   end subroutine PrincipalSig

   !-------------------------------------------------------------------
   subroutine Eig_3(iOpt, St, xN1, xN2, xN3, S1, S2, S3, P, Q, &
      J, theta)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions
      !            from cartesian stress vector
      !
      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
      !
      !  IOpt            I   I     flag for output writing (IOpt = 1)
      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
      !  xN1, xN2, xN3   O   R()   principal direction
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      !  J               O   R     sqrt J2
      !  theta           O   R     lode angle
      !
      !-------------------------------------------------------------------

      implicit none

      ! arguments
      integer, intent(in) :: IOpt
      real(kind = dp), intent(in) :: St(6)
      real(kind = dp), intent(out) :: xN1(3), xN2(3), xN3(3),&
         S1, S2, S3, P, Q, J, theta

      ! local variables
      real(kind = dp) :: A(3,3), V(3,3)
      real(kind = dp) :: abs_max_s, tol
      real(kind = dp) :: tau, sign_tau, t, c, s
      real(kind = dp) :: temp1, temp2, temp3
      integer :: i, k, it, itmax, ip, iq
      integer :: iS1, iS2, iS3

      ! Put cartesian stress vector into matrix A
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      ! Set V to unity matrix
      V(1,1) = 1
      V(2,1) = 0
      V(3,1) = 0

      V(1,2) = 0
      V(2,2) = 1
      V(3,2) = 0

      V(1,3) = 0
      V(2,3) = 0
      V(3,3) = 1

      ! get maximum value of cartesian stress vector
      abs_max_s = 0.0
      do i = 1,6
         if (abs(St(i)) .gt. abs_max_s) abs_max_s = abs(St(i))
      end do

      ! set tolerance
      tol = 1d-16 * abs_max_s

      ! get principal stresses and directions iteratively
      it = 0
      itmax = 50
      do while ( (it .lt. itmax) .and.&  
        (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) .gt. tol) )

         it = it + 1
         do k = 1,3
            if (k .eq. 1) then
               ip = 1
               iq = 2
               else if (k .eq.2) then
                  ip = 2
                  iq = 3
            else
               ip = 1
               iq = 3
            end if

            if (abs(A(ip,iq)) .gt. 1d-50) then
               tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
               if (tau .ge. 0.0) then
                  sign_tau = 1.0
               else
                  sign_tau = -1.0
               end if

               t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
               c = 1.0 / sqrt(1.0 + t**2)
               s = t * c

               temp1 = c * A(1, ip) - s * A(1, iq)
               temp2 = c * A(2, ip) - s * A(2, iq)
               temp3 = c * A(3, ip) - s * A(3, iq)
               A(1, iq) = s * A(1, ip) + c * A(1, iq)
               A(2, iq) = s * A(2, ip) + c * A(2, iq)
               A(3, iq) = s * A(3, ip) + c * A(3, iq)
               A(1, ip) = temp1
               A(2, ip) = temp2
               A(3, ip) = temp3

               temp1 = c * V(1, ip) - s * V(1, iq)
               temp2 = c * V(2, ip) - s * V(2, iq)
               temp3 = c * V(3, ip) - s * V(3, iq)
               V(1, iq) = s * V(1, ip) + c * V(1, iq)
               V(2, iq) = s * V(2, ip) + c * V(2, iq)
               V(3, iq) = s * V(3, ip) + c * V(3, iq)
               V(1, ip) = temp1
               V(2, ip) = temp2
               V(3, ip) = temp3

               temp1 = c * A(ip, 1) - s * A(iq, 1)
               temp2 = c * A(ip, 2) - s * A(iq, 2)
               temp3 = c * A(ip, 3) - s * A(iq, 3)
               A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
               A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
               A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
               A(ip, 1) = temp1
               A(ip, 2) = temp2
               A(ip, 3) = temp3

            end if 
         end do ! A(ip,iq)<>0

      end do ! k

      ! get principal stresses from diagonal of A
      S1 = A(1, 1)
      S2 = A(2, 2)
      S3 = A(3, 3)

      ! derived invariants
      P = (S1 + S2 + S3) / 3.0d0
      Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )
      J = Q /sqrt(3.0d0)
      theta = atan( (1/sqrt(3.0d0)) * ( ( (2*(S2-S3)) / (S1-S3) ) - 1. ) )
      ! if hyrostatic then theta (Load angle) is not defined. 
      ! Since J2 is also zero, make g(theta) be some dummy value - doesnt matter for F 

      if (S1 == S2 .and. S2 == S3) then
         theta = 0 
      end if 

      ! Sort eigenvalues S1 <= S2 <= S3
      iS1 = 1
      iS2 = 2
      iS3 = 3

      if (S1 .gt. S2) then
         t   = S2
         S2  = S1
         S1  = t
         it  = iS2
         iS2 = iS1
         iS1 = it
      end if

      if (S2 .gt. S3) then
         t   = S3
         S3  = S2
         S2  = t
         it  = iS3
         iS3 = iS2
         iS2 = it
      end if

      if (S1 .gt. S2) then
         t   = S2
         S2  = S1
         S1  = t
         it  = iS2
         iS2 = iS1
         iS1 = it
      end if

      ! get corresponding principal directions from V
      do i = 1,3
         xN1(i) = V(i, is1)
         xN2(i) = V(i, is2)
         xN3(i) = V(i, is3)
      end do

      ! optional output writing


   END subroutine Eig_3


   !-------------------------------------------------------------------
   subroutine Eig_3a(iOpt, St, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses from cartesian stress vector
      !
      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
      !
      !  IOpt            I   I     flag for output writing (IOpt = 1)
      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      !
      !-------------------------------------------------------------------

      implicit none

      ! arguments
      integer, intent(in) :: IOpt
      real(kind = dp), intent(in) :: St(6)
      real(kind = dp), intent(out) :: S1, S2, S3, P, Q

      ! local variables
      real(kind = dp) :: A(3,3)
      real(kind = dp) :: abs_max_s, tol
      real(kind = dp) :: tau, sign_tau, t, c, s
      real(kind = dp) :: temp1, temp2, temp3
      integer :: i, k, it, itmax, ip, iq

      ! Put cartesian stress vector into matrix A
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      ! get maximum value of cartesian stress vector
      abs_max_s = 0.0
      do i = 1,6
         if (abs(St(i)) .gt. abs_max_s) abs_max_s = abs(St(i))
      end do

      ! set tolerance
      tol = 1d-20 * abs_max_s

      ! get principal stresses and directions iteratively
      it = 0
      itmax = 50
      do while ( (it .lt. itmax) .and.&     
      (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) .gt. tol) )

         it = it + 1
         do k = 1,3
            if (k .eq. 1) then
               ip = 1
               iq = 2
            else if (k .eq.2) then
               ip = 2
               iq = 3
            else
               ip = 1
               iq = 3
            end if

            if (abs(A(ip,iq)) .gt. 1d-50) then

               tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
               if (tau .ge. 0.0) then
                  sign_tau = 1.0
               else
                  sign_tau = -1.0
               end if

               t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
               c = 1.0 / sqrt(1.0 + t**2)
               s = t * c

               temp1 = c * A(1, ip) - s * A(1, iq)
               temp2 = c * A(2, ip) - s * A(2, iq)
               temp3 = c * A(3, ip) - s * A(3, iq)
               A(1, iq) = s * A(1, ip) + c * A(1, iq)
               A(2, iq) = s * A(2, ip) + c * A(2, iq)
               A(3, iq) = s * A(3, ip) + c * A(3, iq)
               A(1, ip) = temp1
               A(2, ip) = temp2
               A(3, ip) = temp3

               temp1 = c * A(ip, 1) - s * A(iq, 1)
               temp2 = c * A(ip, 2) - s * A(iq, 2)
               temp3 = c * A(ip, 3) - s * A(iq, 3)
               A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
               A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
               A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
               A(ip, 1) = temp1
               A(ip, 2) = temp2
               A(ip, 3) = temp3

            end if ! A(ip,iq)<>0

         end do ! k

         ! optional output writing

      end do ! while

      ! get principal stresses from diagonal of A
      S1 = A(1, 1)
      S2 = A(2, 2)
      S3 = A(3, 3)

      ! derived invariants
      P = (S1 + S2 + S3) / 3.
      Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )

      ! Sort eigenvalues S1 <= S2 <= S3
      if (S1 .gt. S2) then
         t   = S2
         S2  = S1
         S1  = t
      end if

      if (S2 .gt. S3) then
         t   = S3
         S3  = S2
         S2  = t
      end if

      if (S1 .gt. S2) then
         t   = S2
         S2  = S1
         S1  = t
      end if

      ! optional output writing

   end subroutine Eig_3a

       !> \brief Calculates the new stress based on the given strain, stiffness matrix, and initial stress.
    !!
    !! This function calculates the new stress state by computing the stress increment
    !! using the stiffness matrix and strain, and then adding this increment to the initial stress.
    !!
    !! \param[in] strain Strain vector \f$ \epsilon \f$ of size (6).
    !! \param[in] stiff_matrix Stiffness matrix \f$ DE \f$ of size (6, 6).
    !! \param[in] init_stress Initial stress vector \f$ \sigma_0 \f$ of size (6).
    !! \return New stress vector \f$ \sigma \f$ of size (6).
    !!
   function calc_elastic_stress(strain, stiff_matrix, init_stress) result(new_stress)
      real(kind = dp), intent(in) :: strain(6), stiff_matrix(6,6), init_stress(6)
      real(kind = dp) :: new_stress(6)

      ! Local variables
      real(kind =dp) :: stress_inc(6) ! Stress increment

      ! Calc the stress increment
      stress_inc = calc_elastic_stress_increment(strain, stiff_matrix)

      ! Calc the new stress
      new_stress = init_stress + stress_inc
  end function


  subroutine stiffnessMCC(p,xe0,xkappa,xNu,xG,xK)
   !-----------------------------------------------------------------------
   !	input:	p,xe0,xkappa,xNu
   !	output:	xG,xK
   !	local:	r (xG-xK ratio)
   !-----------------------------------------------------------------------
   implicit none
   real(kind = dp), intent(in) :: p,xe0,xkappa,xNu
   real(kind = dp), intent(out) :: xG,xK
   real(kind = dp) :: r, abs_p

   abs_p = abs(p)
   xK = (xe0 + 1)/xkappa * abs_p  !bulk modulus at start of time step, assumed constant
   r = 3.0_dp * ( 1.0_dp - 2.0_dp*xnu) / ( 2.0_dp * (1.0_dp+xnu))
   xG = r*xK !shear modulus

end subroutine stiffnessMCC


   !	Subroutine to form elastic D matrix
   !-------------------------------------------------------------------
Subroutine FormDEMCC(Sig0, xkappa,xNu, xN, xlambda, D, Id, xG, xK)
   !-------------------------------------------------------------------
   !    Function:  To form the elastic material stiffness matrix for MCC model (Hooke)
   ! Inputs:
   !        xkappa:         Slope of swelling line (U/R line) in e-ln(p') plane
   !           xNu:         Poisson's ratio
   !           xe0:         Initial void ratio
   !            Id:         (First) dimension of D
   ! Outputs:
   !            xG:         Shear modulus
   !            xK:         Bulk modulus
   !        D(i,j):         Resulting matrix
   !                              D1  D2  D2 o  o  o
   !     Structure of             D2  D1  D2 o  o  o
   !     elastic D matrix         D2  D2  D1 o  o  o
   !                              o   o   o  G  o  o
   !                              o   o   o  o  G  o
   !                              o   o   o  o  o  G
   !-------------------------------------------------------------------

   real(kind = dp) :: Sig0(6), D, xG, xK, v
   real(kind = dp) ,intent(in) :: xNu,xN,xlambda,xkappa
   integer(kind = i32) :: Id, I,J
   Dimension D(Id,Id)
   real(kind = dp) :: p,r,FAC,D1,D2
   
   D = 0.0
   P = max(-1.0_dp * calc_mean_stress(Sig0), 1.0_dp)

   v = xN - (xlambda*(log(P)))
   !v = 0.5 +1
   xK = ((v)/xkappa) * P  !bulk modulus at start of time step, assumed constant
   r = 3. * ( 1. - 2.*xNu) / ( 2. * (1.+xNu))
   xG = r*xK
   FAC= 2*xG / (1D0 - 2*xNU)
   D1 = FAC * (1D0 - xNU)
   D2 = FAC *   xNU

   Do I=1,3
      Do J=1,3
         D(I,J)=D2
      End Do
      D(I,I)=D1
   End Do

   Do I=4,6
      D(I,I)=xG
   End Do

End Subroutine FormDEMCC
