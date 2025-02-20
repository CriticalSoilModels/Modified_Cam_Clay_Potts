module mod_implicit_stress_integ
   use kind_precision_module, only: dp, i32
   use mod_elastic_stress, only: construct_stiffness_matrix, calc_shear_modulus
   use mod_stress_invariants, only: Get_invariants
   use mod_state_params, only: calc_cam_clay_bulk_modulus

   implicit none

contains
   ! pure function calc_yield_function(mean_stress, J, lode_angle, phi_cs, preconsol_stress) result(F)
   !    real(kind = dp), intent(in) :: mean_stress, J, lode_angle, phi_cs, preconsol_stress
   !    real(kind = dp) :: F


   ! if (mean_stress > (preconsol_stress/2.0_dp)) then
   !    ! If p is greater than pp/2, then Modified Cam clay
   !    F = yield_MCC_log(p, j, pp, theta, xphics)
   ! else
   !    F = yield_Hvor_log(p, j, pp, theta, xphics,m_Hrov, xgamma,xlambda,xkappa,xN)
   ! end if

   ! end function evaluate_yield_function

   subroutine implicit_predictor_corrector_integration(xkappa,XNu,&
      dEps,xphics,xlambda,m_Hrov, xgamma,xN,FTOL,MaxIter,Sig,EpsP,dEpsP,pp,v)
      !-------------------------------------------------------------------
      !	input
      !        xkappa:         Slope of swelling line (U/R line) in e-ln(p') plane
      !           xNu:         Poisson's ratio
      !           xe0:         Initial void ratio
      !          dEps:         Strain increment
      !        xphics:         Critical state angle of shearing resistance, in radians
      !          FTOL:         Tolerence on the yield surface
      !       MaxIter:         Maximum iterations for the algorithm
      !          zeta:         (v / lambda + kappa)
      !---------------------------------------------------------------
      !	input/output:
      !            Sig:        Updated stress (Voight notation)
      !           EpsP:        Plastic strain
      !          dEpsP:        Incremental plastic strain
      !             pp:        State variable - preconsolidation pressure
      !---------------------------------------------------------------
      !	output:
      !          dEpsP:        Incremental plastic strain
      !---------------------------------------------------------------
      !	local:
      !           Sigu:        Local variables for calculation - Current stress
      !          EpsPu:        Local variables for calculation - Plastic strain
      !         dEpsPu:        Local variables for calculation - Incremental plastic strain
      !              D:        Elastic D Matrix
      !             xG:        Shear Mod
      !             xK:        Bulk Mod
      !          dSigu:        Stress increment (total)
      !           iOpt:        Integer for eigen value options
      !  xN1, xN2, xN3:        Eigen directions (not used directly)
      !     S1, S2, S3:        Principal stresses
      !    p,q,j,theta:        Invariants
      !              F:        Yield surface value
      !        counter:        Counter for the implicit algorithm (exits if higher)
      !           dgdp:        Derivative of plastic potential with pressure = plastic volumetric strain
      !         dfdsig:        Derivative of yield potential with current stress state (n vector)
      !         dgdsig:        Derivative of plastic potential with current stress state (m vector)
      !          sqrt3:        Square root of 3
      !         gtheta:        used in 3D algorithm in yield surface instead of M
      !        result1:        dummy variable for calculations
      !        result2:        dummy variable for calculations
      !        dlambda:        Lambda dot
      !              A:        Hardening/Softening parameter
      !              v:        specific volume

      !-------------------------------------------------------------------
      implicit none
      ! Input variables
      real(kind = dp),intent(in) :: xkappa,xNu,xphics,FTOL
      real(kind = dp),intent(in) :: xlambda,m_Hrov, xgamma,xN
      integer(kind = i32) :: MaxIter
      real(kind = dp), intent(in), dimension(6)  :: dEps

      ! Input/Output variables
      real(kind = dp), intent(inout), dimension(6)  :: Sig,EpsP
      real(kind = dp), intent(inout)  :: pp,v

      ! Output variables
      real(kind = dp), intent(inout), dimension(6)  :: dEpsP

      ! Local variables
      real(kind = dp)  :: xK,xG,xN1(3),xN2(3),xN3(3),S1,S2,S3
      real(kind = dp)  :: p,q,j,theta,F,dgdp,sqrt3,zeta,beta
      real(kind = dp)  :: gtheta,result1(6),result2,dlambda,A
      real(kind = dp), dimension(6)  :: Sigu,EpsPu,dSigu
      real(kind = dp), dimension(6)  :: dfdsig,dgdsig
      real(kind = dp),dimension(6,6) :: D
      integer(kind = i32) :: iOpt,counter
      logical :: check

      ! Initialization
      DEpsP = 0.0d0
      F = 0.0d0

      !Store variables for updating
      Sigu = Sig
      EpsPu = EpsP
      beta = 0.75


      call Get_invariants(Sigu, p, j, theta)


      ! Calc the stiffness matrix

      ! Calc the camclay update for the bulk modulus
      K = calc_cam_clay_bulk_modulus(Sigu, lambda, N, kappa)

      ! Calc the shear modulus
      G = calc_shear_modulus(nu, K)

      D = construct_stiffness_matrix(G, nu)

      dSigu = matmul(D, dEps)
      Sigu = Sigu + dSigu

      iOpt = 1

      call PrincipalSig(iOpt, Sigu, xN1, xN2, xN3, S1, S2, S3, p, q, &
         j, theta)

      p = max(-p, 1.)

      v = xN - (xlambda*(log(p)) )

      zeta = v / (xlambda - xkappa)

      if (p > (pp/2.0)) then
         ! If p is greater than pp/2, then Modified Cam clay
         F = yield_MCC_log(p, j, pp, theta, xphics)
      else
         F = yield_Hvor_log(p, j, pp, theta, xphics,m_Hrov, xgamma,xlambda,xkappa,xN)
      end if

      if (F <= FTOL) then
         ! Prediction of the stress and strain values are correct and the values can be updated and returned
         ! Update Sig, EpsP, dEpsP
         Sig = Sigu
         EpsP(:) = 0
         dEpsP(:) = 0

         ! Update state parameters values
         pp = pp

         ! Exit the subroutine
         return
      end if

      ! Initialize the counter keep check on the iterations
      counter = 0
      check = .true.

      do while(check .and. counter <= MaxIter)

         if (p > (pp/2.0)) then
            ! If p is greater than pp/2, then Modified Cam clay
            call derivatives_MCC_log(Sigu,p,j,xphics,theta,pp,dgdp,dfdsig,dgdsig)
         else
            call derivatives_Hvorslev_log(Sigu,p,j,xphics,theta,pp,m_Hrov,xgamma, &
               xlambda,xkappa,xN,dgdp,dfdsig,dgdsig)

         end if

         ! Calculate A (make function later)
         sqrt3 = sqrt(3.0d0)
         gtheta = cos(theta) + ((sin(theta) * sin(xphics)) / sqrt3)
         gtheta = sin(xphics) / gtheta

         if (p > (pp/2.0)) then
            ! If p is greater than pp/2, then Modified Cam clay
            !A = zeta * (pp/(p**2))* (1 - (( j / (p*gtheta) )**2) )
            A = zeta * (pp) * dgdp * (-(gtheta**2)*(p))
            !A = (((gtheta**2)*pp)* ( ((gtheta**2)*(p**2)) - (3*(j**2))) ) &
            !         / (xlambda-xkappa)
         else
            !A = ((v)/ xlambda) * exp((xgamma - xN - (xkappa*log(pp/p)))/ xlambda) &
            !         * (pp/p) * (1 - (( j / (p*gtheta) )**2) )
            !A = - A
            A = (2.0*(1.0-beta)*(((gtheta**2)*(p**2)) - &
               (q**2)))/(xlambda*xkappa*p)
            A = zeta * dgdp * (1.0 - (xkappa/xlambda)) * (1.0-beta)

         end if


         ! n * D * m = dfdsig * D * dgdsig
         Call FormDEMCC(Sigu, xkappa,xNu, xN, xlambda, D, 6, xG, xK)
         result1 = matmul(D,dgdsig)
         result2 = dot_product(dfdsig, result1)

         ! dlambda
         dlambda = F / (result2 + A)

         ! Update the stress
         Sigu = Sigu + (dlambda*result1)


         ! Acc plastic strain
         EpsPu = EpsPu + (dlambda * dgdsig)

         ! Update the state parameters (pp)
         pp = pp * exp(zeta * dlambda * dgdp)

         ! Calc the yield function value
         iOpt = 1

         call PrincipalSig(iOpt, Sigu, xN1, xN2, xN3, S1, S2, S3, p,&
            q,j,theta)

         p = max(-p, 1.)

         v = xN - (xlambda*(log(p)) )
         !v = 0.5 + 1
         zeta = v / (xlambda - xkappa)

         if (p > (pp/2.0)) then
            ! If p is greater than pp/2, then Modified Cam clay
            F = yield_MCC_log(p, j, pp, theta, xphics)
         else
            F = yield_Hvor_log(p, j, pp, theta, xphics,m_Hrov, xgamma, xlambda, xkappa, xN)
         end if

         if (F <= FTOL) then
            check = .false.
         end if

         ! Update the counter
         counter = counter + 1

      end do

      ! Return the integrated parameters
      Sig = Sigu
      dEpsP = EpsPu-EpsP
      EpsP = EpsPu

   end subroutine implicit_predictor_corrector_integration


   ! subroutine implicit_predictor_corrector_integration(xkappa,XNu,&
   !    dEps,xphics,xlambda,m_Hrov, xgamma,xN,FTOL,MaxIter,Sig,EpsP,dEpsP,pp,v)
   !    !-------------------------------------------------------------------
   !    !	input
   !    !        xkappa:         Slope of swelling line (U/R line) in e-ln(p') plane
   !    !           xNu:         Poisson's ratio
   !    !           xe0:         Initial void ratio
   !    !          dEps:         Strain increment
   !    !        xphics:         Critical state angle of shearing resistance, in radians
   !    !          FTOL:         Tolerence on the yield surface
   !    !       MaxIter:         Maximum iterations for the algorithm
   !    !          zeta:         (v / lambda + kappa)
   !    !---------------------------------------------------------------
   !    !	input/output:
   !    !            Sig:        Updated stress (Voight notation)
   !    !           EpsP:        Plastic strain
   !    !          dEpsP:        Incremental plastic strain
   !    !             pp:        State variable - preconsolidation pressure
   !    !---------------------------------------------------------------
   !    !	output:
   !    !          dEpsP:        Incremental plastic strain
   !    !---------------------------------------------------------------
   !    !	local:
   !    !           Sigu:        Local variables for calculation - Current stress
   !    !          EpsPu:        Local variables for calculation - Plastic strain
   !    !         dEpsPu:        Local variables for calculation - Incremental plastic strain
   !    !              D:        Elastic D Matrix
   !    !             xG:        Shear Mod
   !    !             xK:        Bulk Mod
   !    !          dSigu:        Stress increment (total)
   !    !           iOpt:        Integer for eigen value options
   !    !  xN1, xN2, xN3:        Eigen directions (not used directly)
   !    !     S1, S2, S3:        Principal stresses
   !    !    p,q,j,theta:        Invariants
   !    !              F:        Yield surface value
   !    !        counter:        Counter for the implicit algorithm (exits if higher)
   !    !           dgdp:        Derivative of plastic potential with pressure = plastic volumetric strain
   !    !         dfdsig:        Derivative of yield potential with current stress state (n vector)
   !    !         dgdsig:        Derivative of plastic potential with current stress state (m vector)
   !    !          sqrt3:        Square root of 3
   !    !         gtheta:        used in 3D algorithm in yield surface instead of M
   !    !        result1:        dummy variable for calculations
   !    !        result2:        dummy variable for calculations
   !    !        dlambda:        Lambda dot
   !    !              A:        Hardening/Softening parameter
   !    !              v:        specific volume

   !    !-------------------------------------------------------------------
   !    implicit none
   !    ! Input variables
   !    real(kind = dp),intent(in) :: xkappa,xNu,xphics,FTOL
   !    real(kind = dp),intent(in) :: xlambda,m_Hrov, xgamma,xN
   !    integer(kind = i32) :: MaxIter
   !    real(kind = dp), intent(in), dimension(6)  :: dEps

   !    ! Input/Output variables
   !    real(kind = dp), intent(inout), dimension(6)  :: Sig,EpsP
   !    real(kind = dp), intent(inout)  :: pp,v

   !    ! Output variables
   !    real(kind = dp), intent(inout), dimension(6)  :: dEpsP

   !    ! Local variables
   !    real(kind = dp)  :: xK,xG,xN1(3),xN2(3),xN3(3),S1,S2,S3
   !    real(kind = dp)  :: p,q,j,theta,F,dgdp,sqrt3,zeta,beta
   !    real(kind = dp)  :: gtheta,result1(6),result2,dlambda,A
   !    real(kind = dp), dimension(6)  :: Sigu,EpsPu,dSigu
   !    real(kind = dp), dimension(6)  :: dfdsig,dgdsig
   !    real(kind = dp),dimension(6,6) :: D
   !    integer(kind = i32) :: iOpt,counter
   !    logical :: check

   !    ! Initialization
   !    DEpsP = 0.0d0
   !    F = 0.0d0

   !    !Store variables for updating
   !    Sigu = Sig
   !    EpsPu = EpsP
   !    beta = 0.75

   !    ! Update G,K and evaluate D
   !    Call FormDEMCC(Sigu, xkappa,xNu, xN, xlambda, D, 6, xG, xK)
   !    write(10, '(A, F10.5)') "Value of G at start of imp func: ", xG
   !    write(10, '(A, E10.5)') "Value of K at start of imp func: ", xK

   !    call MatVec(D, 6, dEps, 6, dSigu)

   !    write(10, '(A, 6F10.5)') "Value of dSigu: ", dSigu
   !    Sigu = Sigu + dSigu


   !    ! Calculate F for the updated Sig - find variants first
   !    iOpt = 1
   !    call PrincipalSig(iOpt, Sigu, xN1, xN2, xN3, S1, S2, S3, p, q, &
   !       j, theta)
   !    p = max(-p, 1.)

   !    v = xN - (xlambda*(log(p)) )
   !    zeta = v / (xlambda - xkappa)

   !    write(10, '(A, F10.5)') "Value of p at start of imp func: ", p
   !    write(10, '(A, F10.5)') "Value of j at start of imp func: ", j

   !    if (p > (pp/2.0)) then
   !       ! If p is greater than pp/2, then Modified Cam clay
   !       F = yield_MCC_log(p, j, pp, theta, xphics)
   !    else
   !       F = yield_Hvor_log(p, j, pp, theta, xphics,m_Hrov, xgamma,xlambda,xkappa,xN)
   !    end if

   !    write(10, '(A, E20.5)') "Value of F at start of imp func: ", F

   !    if (F <= FTOL) then
   !       ! Prediction of the stress and strain values are correct and the values can be updated and returned
   !       ! Update Sig, EpsP, dEpsP
   !       Sig = Sigu
   !       EpsP(:) = 0
   !       dEpsP(:) = 0

   !       ! Update state parameters values
   !       pp = pp

   !       ! Exit the subroutine
   !       return
   !    end if

   !    ! Initialize the counter keep check on the iterations
   !    counter = 0
   !    check = .true.

   !    do while(check .and. counter <= MaxIter)

   !       if (p > (pp/2.0)) then
   !          ! If p is greater than pp/2, then Modified Cam clay
   !          call derivatives_MCC_log(Sigu,p,j,xphics,theta,pp,dgdp,dfdsig,dgdsig)
   !       else
   !          call derivatives_Hvorslev_log(Sigu,p,j,xphics,theta,pp,m_Hrov,xgamma, &
   !             xlambda,xkappa,xN,dgdp,dfdsig,dgdsig)

   !       end if

   !       ! Calculate A (make function later)
   !       sqrt3 = sqrt(3.0d0)
   !       gtheta = cos(theta) + ((sin(theta) * sin(xphics)) / sqrt3)
   !       gtheta = sin(xphics) / gtheta

   !       if (p > (pp/2.0)) then
   !          ! If p is greater than pp/2, then Modified Cam clay
   !          !A = zeta * (pp/(p**2))* (1 - (( j / (p*gtheta) )**2) )
   !          A = zeta * (pp) * dgdp * (-(gtheta**2)*(p))
   !          !A = (((gtheta**2)*pp)* ( ((gtheta**2)*(p**2)) - (3*(j**2))) ) &
   !          !         / (xlambda-xkappa)
   !       else
   !          !A = ((v)/ xlambda) * exp((xgamma - xN - (xkappa*log(pp/p)))/ xlambda) &
   !          !         * (pp/p) * (1 - (( j / (p*gtheta) )**2) )
   !          !A = - A
   !          A = (2.0*(1.0-beta)*(((gtheta**2)*(p**2)) - &
   !             (q**2)))/(xlambda*xkappa*p)
   !          A = zeta * dgdp * (1.0 - (xkappa/xlambda)) * (1.0-beta)

   !       end if


   !       ! n * D * m = dfdsig * D * dgdsig
   !       Call FormDEMCC(Sigu, xkappa,xNu, xN, xlambda, D, 6, xG, xK)
   !       result1 = matmul(D,dgdsig)
   !       result2 = dot_product(dfdsig, result1)

   !       ! dlambda
   !       dlambda = F / (result2 + A)

   !       ! Update the stress
   !       Sigu = Sigu + (dlambda*result1)


   !       ! Acc plastic strain
   !       EpsPu = EpsPu + (dlambda * dgdsig)

   !       ! Update the state parameters (pp)
   !       pp = pp * exp(zeta * dlambda * dgdp)

   !       ! Calc the yield function value
   !       iOpt = 1

   !       call PrincipalSig(iOpt, Sigu, xN1, xN2, xN3, S1, S2, S3, p,&
   !          q,j,theta)

   !       p = max(-p, 1.)

   !       v = xN - (xlambda*(log(p)) )
   !       !v = 0.5 + 1
   !       zeta = v / (xlambda - xkappa)

   !       if (p > (pp/2.0)) then
   !          ! If p is greater than pp/2, then Modified Cam clay
   !          F = yield_MCC_log(p, j, pp, theta, xphics)
   !       else
   !          F = yield_Hvor_log(p, j, pp, theta, xphics,m_Hrov, xgamma, xlambda, xkappa, xN)
   !       end if

   !       if (F <= FTOL) then
   !          check = .false.
   !       end if

   !       ! Update the counter
   !       counter = counter + 1

   !    end do

   !    ! Return the integrated parameters
   !    Sig = Sigu
   !    dEpsP = EpsPu-EpsP
   !    EpsP = EpsPu

   ! end subroutine implicit_predictor_corrector_integration

end module mod_implicit_stress_integ
