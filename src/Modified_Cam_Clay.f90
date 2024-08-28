module mod_UMAT_Cam_Clay
   use,intrinsic :: iso_fortran_env, only: dp => real64
   use,intrinsic :: iso_fortran_env, only: i32 => int32

   implicit none

   ! Modified cam clay with mohr-coulomb sufrace in the devaitoric plane
   ! Theory and equations derived from "Finite element analysis in geotechncial engineering", by David M. Potts and Lidija Zdrakovic, 1999

contains

   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------

   ! *USER SUBROUTINE
   SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
      RPL,DDSDDT,DRPLDE,DRPLDT,&
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      !---------------------------------------------------------------
      !	Input
      !         PROPS:         Material parameters
      !        DSTRAN:         Strain increment
      !---------------------------------------------------------------
      !	Input/output:
      !         STRESS:        Stress (Voight notation - S1,S2,S3,S12,S23,S13)
      !         STATEV:        State variables
      !          dEpsP:        Incremental plastic strain
      !             pp:        State variable - preconsolidation pressure
      !---------------------------------------------------------------
      !	Output:
      !         DDSDDE:        Material stiffness matrix
      !---------------------------------------------------------------
      !	Other:
      ! Other common variables not used for this UMAT. For more information please visit
      ! https://www.brown.edu/Departments/Engineering/Courses/En2340/Programming/ABAQUS/Usermat.for
      !
      !         DDSDDT:        Variation of the stress increments with respect to the temperature.
      !            RPL:        Volumetric heat generation per unit time at the end of the increment
      !                        caused by mechanical working of the material
      !         DRPLDE:        Variation of RPL with respect to the strain increments.
      !           TIME:        TIME(1): Value of step time at the beginning of the current increment.
      !                        TIME(2): Value of total time at the beginning of the current increment.
      !         PREDEF:        Array of interpolated values of predefined field variables at this point
      !                        at the start of the increment, based on the values read in at the nodes.
      !          DPRED:        Array of increments of predefined field variables.
      !         COORDS:        An array containing the coordinates of this point. These are the current
      !                        coordinates if geometric nonlinearity is accounted for during the step
      !           DROT:        Rotation increment matrix. This matrix represents the increment of rigid
      !                        body rotation of the basis system in which the components of stress
      !                        (STRESS) and strain (STRAN) are stored.
      !         DFGRD0:        DFGRD0(3,3): Array containing the deformation gradient at the beginning of the increment.
      !         DFGRD1:        DFGRD1(3,3): Array containing the deformation gradient at the end of the increment.
      !---------------------------------------------------------------
      !	Local: Material parameters
      !         xphics:        Critical state angle of shearing resistance
      !            xNu:        Poisson's ratio
      !         xkappa:        Slope of swelling line (U/R line) in e-ln(p') plane
      !        xlambda:        Slope of virgin consolidation line in e-ln(p') plane
      !            xe0:        Initial void ratio
      !           zeta:        (v / lambda + kappa)
      !             xG:        Shear Modulus
      !             xK:        Bulk Modulus
      ! Local: State parameters
      !             pp:        Preconsolidation pressure
      !           EpsP:        Accumulated plastic strain
      ! Local: Stresses and strains
      !            Sig:        Stress (Voight notation)
      !           dEps:        Strain increment
      !          dEpsP:        Plastic strain in the time increment
      ! Local: Implicit algorithm inputs for tolerence and control
      !          FTOL:         Tolerence on the yield surface
      !       MaxIter:         Maximum iterations for the algorithm
      !---------------------------------------------------------------

      ! DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
      ! INCLUDE 'ABA_PARAM.INC'
      implicit none
      integer(kind = i32) :: NTENS,NPROPS,NSTATEV, NDI, NSHR, NOEL, NPT, LAYER, KSPT,KSTEP,KINC
      real(kind = dp) :: STRESS,STATEV,DDSDDE,SSE,SPD,SCD,        &
                         RPL,DDSDDT,DRPLDE,DRPLDT, STRAN,DSTRAN,  &
                         TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED, &
                         PROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,  &
                         DFGRD1

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),&
         DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),&
         STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),&
         PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      ! Userdefined parameters
      real(kind = dp) :: xphics, xNu, xkappa, xlambda, xG, xK
      real(kind = dp) :: m_Hrov,xgamma,xN
      real(kind = dp) :: pp,v, EpsP(6), FTOL
      real(kind = dp), dimension(6) :: Sig, dEps, dEpsP
      real(kind = dp) :: pwp,porosity,Bulk_Water
      integer(kind = i32) :: MaxIter, i, status
      character(len=100) :: output_file
      logical :: file_exists
      
      ! Make sure all of the variables are used
      if (.False.) then
         ! This condiiton is purposefully set to zero so that the compiler is 
         ! is tricked into thinking that all the variables are used.
         print *, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, STRAN,    &
                  TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, &
                  NSHR, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,   &
                  NOEL, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      end if

      !-------------------------------------------------------------------
      ! SUBROUTINE START
      !-------------------------------------------------------------------

      ! Initialization: Get parameters from Props and STATEV
      xphics = Props(1)
      xNu = Props(2)
      xkappa = Props(3)
      xlambda = Props(4)
      xN = Props(5)
      m_Hrov = Props(6)
      xgamma = xN - ((xlambda-xkappa)*log(2.0))

      pp = STATEV(1)
      v = STATEV(2)
      pwp = STATEV(3)
      ! Accumulated plastic strains
      do i = 1,NTENS
         EpsP(i) = STATEV(3+i)
      end do
      ! helping parameter


      ! Initialize stress and strain
      Sig=stress
      dEps=dstran
      do i = 1, 3
         dEps(i) = -dEps(i)
      end do
      
      porosity = 0.6
      Bulk_Water = 2100000
      pwp =  pwp + ((Bulk_Water/porosity)*(dEps(1)+dEps(2)+dEps(3)))
   
      ! Diagnostics file for debugging (optional)
      ! Make a file named diagnostics_output.txt in the folder or else you will have error 
      ! Erase data in file before every run
      output_file = 'diagnostics_output.txt' 
      open(unit=10, file=output_file, status='old', action='write', position='append', iostat=status)
      

      if (status /= 0) then
         print *, "Error opening file for writing."
         stop
      end if

      ! Write header information
      write(10, '(A)') "Diagnostics Output"
      write(10, '(A)') "=================="

      write(10, '(A, 6F10.5)') "Value of sig: ", Sig
      write(10, '(A, 6F10.5)') "Value of axial strain: ", dEps


      !-------------------------------------------------------------------
      ! Do the predictor corrector scheme.
      ! The subroutine calculate the plastic and elastic parts and returns the updated stress and state variables
      ! Set tolerance for yield surface and the maximum iterations the algorithm can do
      ! Reccommended tolerance error (10-6 to 10-9)
      FTOL = 1e-3
      MaxIter = 10
      call implicit_predictor_corrector_integration(xkappa,XNu,&
      dEps,xphics,xlambda,m_Hrov, xgamma,xN,FTOL,MaxIter,Sig,EpsP,dEpsP,pp,v)

      !-------------------------------------------------------------------
      ! update state variables
      STATEV(1) = pp
      STATEV(2) = v
      STATEV(3) = pwp
      do i = 1,NTENS
         STATEV(3+i) = EpsP(i)
      end do
      !-------------------------------------------------------------------
      ! if (isundr  ==  1) then        Calculation of pore pressure not needed because done outside the subroutine
      !   Swp = Swp0 - BulkW*(dEpsV)
      ! end if

      !-------------------------------------------------------------------
      ! update stress and stiffness matrix
      stress=Sig

      write(10, '(A, 6F10.5)') "Value of integrated stress at end: ", stress

      ! Calculate effective/elastic D-matrix
      Call FormDEMCC(stress, xkappa,xNu, xN, xlambda, DDSDDE, 6, xG, xK) ! also updates K and G
      write(10, '(A, F10.5)') "Value of G at end: ", xG
      write(10, '(A, F10.5)') "Value of K at end: ", xK

      write(10, '(A)') "End of step"
      write(10, '(A)') "=================="
      close(10)

   End SUBROUTINE UMAT


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------


   !	Subroutine for implicit stress estimate
   !-------------------------------------------------------------------
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

      ! Update G,K and evaluate D
      Call FormDEMCC(Sigu, xkappa,xNu, xN, xlambda, D, 6, xG, xK)
      write(10, '(A, F10.5)') "Value of G at start of imp func: ", xG
      write(10, '(A, E10.5)') "Value of K at start of imp func: ", xK

      call MatVec(D, 6, dEps, 6, dSigu)

      write(10, '(A, 6F10.5)') "Value of dSigu: ", dSigu
      Sigu = Sigu + dSigu


      ! Calculate F for the updated Sig - find variants first
      iOpt = 1
      call PrincipalSig(iOpt, Sigu, xN1, xN2, xN3, S1, S2, S3, p, q, &
         j, theta)
         p = max(-p, 1.)

      v = xN - (xlambda*(log(p)) )
      zeta = v / (xlambda - xkappa) 

      write(10, '(A, F10.5)') "Value of p at start of imp func: ", p
      write(10, '(A, F10.5)') "Value of j at start of imp func: ", j

      if (p > (pp/2.0)) then
         ! If p is greater than pp/2, then Modified Cam clay 
         F = yield_MCC_log(p, j, pp, theta, xphics)
      else 
         F = yield_Hvor_log(p, j, pp, theta, xphics,m_Hrov, xgamma,xlambda,xkappa,xN)
      end if
      
      write(10, '(A, E20.5)') "Value of F at start of imp func: ", F

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

         write(10, '(A, I5)') "Value of counter : ", counter
         write(10, '(A, F10.5)') "Value of F before derivative: ", F
         write(10, '(A, 6F10.5)') "Value of Sigu before derivative: ", Sigu
         write(10, '(A, F10.5)') "Value of p before derivative:: ", p
         write(10, '(A, F10.5)') "Value of j before derivative:: ", j
         write(10, '(A, F10.5)') "Value of theta before derivative:: ", theta
         write(10, '(A, F10.5)') "Value of pp before derivative:: ", pp

      
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

         write(10, '(A, F10.5)') "Value of A:: ", A

         ! n * D * m = dfdsig * D * dgdsig
         Call FormDEMCC(Sigu, xkappa,xNu, xN, xlambda, D, 6, xG, xK)
         result1 = matmul(D,dgdsig)
         result2 = dot_product(dfdsig, result1)

         write(10, '(A, 6E15.5)') "Value of result1 before lambda:: ", result1
         write(10, '(A, E15.5)') "Value of result2 before lambda:: ", result2

         ! dlambda
         dlambda = F / (result2 + A)

         write(10, '(A, F10.5)') "Value of dlambda:: ", dlambda

         write(10, '(A, 6F10.5)') "Value of Sigu:: ", Sigu

         write(10, '(A, 6F10.5)') "Value of dlambda*result1:: ", dlambda*result1

         ! Update the stress
         Sigu = Sigu + (dlambda*result1)

         write(10, '(A, 6F10.5)') "Value of Sigu:: ", Sigu

         ! Acc plastic strain
         EpsPu = EpsPu + (dlambda * dgdsig)

         ! Update the state parameters (pp)
         pp = pp * exp(zeta * dlambda * dgdp)

         write(10, '(A, 6F10.5)') "Value of pp:: ", pp

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


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------


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
      P = MAX(-(Sig0(1)+Sig0(2)+Sig0(3))/3., 1d0)
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


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------



   !-----------------------------------------------------------------------
   ! Function computing the yield function for MOdified cam clay 
   !-----------------------------------------------------------------------
   function yield_MCC(p, j, pp, theta, xphics)
      !-----------------------------------------------------------------------
      !  function: computing the yield function
      !
      !  input:	p, q, pp, g_theta
      !  output:	yield
      !-----------------------------------------------------------------------
      implicit none

      real(kind = dp) :: p, j, pp, theta, xphics, g_theta
      real(kind = dp) :: yield_MCC

      !write(10, '(A, F10.5)') "Value of theta for F calculation: ", theta

      g_theta = cos(theta) + ( (sin(theta) * sin(xphics)) / sqrt(3.0d0) )

      !write(10, '(A, F10.5)') "Value of g_theta denominator for F calculation: ", g_theta

      g_theta = sin(xphics) / g_theta

      !write(10, '(A, F10.5)') "Value of g_theta for F calculation: ", g_theta

      yield_MCC = ( j / (p*g_theta) )**2 - ((pp/p) - 1)

      !write(10, '(A, F10.5)') "Value of yield for F calculation: ", yield

   end function yield_MCC


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------

   !-----------------------------------------------------------------------
   ! Function computing the yield function for MOdified cam clay 
   !-----------------------------------------------------------------------
   function yield_MCC_log(p, j, pp, theta, xphics)
      !-----------------------------------------------------------------------
      !  function: computing the yield function
      !
      !  input:	p, q, pp, g_theta
      !  output:	yield
      !-----------------------------------------------------------------------
      implicit none

      real(kind = dp) :: p, j, pp, theta, xphics, g_theta
      real(kind = dp) :: yield_MCC_log

      !write(10, '(A, F10.5)') "Value of theta for F calculation: ", theta

      g_theta = cos(theta) + ( (sin(theta) * sin(xphics)) / sqrt(3.0d0) )

      !write(10, '(A, F10.5)') "Value of g_theta denominator for F calculation: ", g_theta

      g_theta = sin(xphics) / g_theta

      !write(10, '(A, F10.5)') "Value of g_theta for F calculation: ", g_theta

      yield_MCC_log = ((j**2) / 3.0d0) - ((g_theta**2)*p*pp) + ((g_theta**2)*(p**2))

      !write(10, '(A, F10.5)') "Value of yield for F calculation: ", yield

   end function yield_MCC_log

 


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------


   !-----------------------------------------------------------------------
   ! Function computing the yield function for Hvorslev surface
   !-----------------------------------------------------------------------
   function yield_Hvor(p, j, pp, theta, xphics,m_Hrov, xgamma, xlambda, xkappa, xN)
      !-----------------------------------------------------------------------
      !  function: computing the yield function
      !
      !  input:	p, q, pp, g_theta, m_Hrov, xgamma, xlambda, xkappa, xN 
      !  output:	yield
      !-----------------------------------------------------------------------
      implicit none

      real(kind = dp) :: p, j, pp, theta, xphics, g_theta
      real(kind = dp) :: m_Hrov, xgamma, xlambda, xkappa, xN 
      real(kind = dp) :: yield_Hvor

      !write(10, '(A, F10.5)') "Value of theta for F calculation: ", theta

      g_theta = cos(theta) + ( (sin(theta) * sin(xphics)) / sqrt(3.0d0) )

      !write(10, '(A, F10.5)') "Value of g_theta denominator for F calculation: ", g_theta

      g_theta = sin(xphics) / g_theta

      !write(10, '(A, F10.5)') "Value of g_theta for F calculation: ", g_theta

      yield_Hvor = ( j / (g_theta - m_Hrov) ) - (( m_Hrov / (g_theta - m_Hrov) )* p ) &
                - (pp*(exp((xgamma - xN)/xlambda))*(exp((-xkappa/xlambda)*log(pp/p))))

      !write(10, '(A, F10.5)') "Value of yield for F calculation: ", yield

   end function yield_Hvor


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------

   !-----------------------------------------------------------------------
   ! Function computing the yield function for Hvorslev surface
   !-----------------------------------------------------------------------
   function yield_Hvor_log(p, j, pp, theta, xphics,m_Hrov, xgamma, xlambda, xkappa, xN)
      !-----------------------------------------------------------------------
      !  function: computing the yield function
      !
      !  input:	p, q, pp, g_theta, m_Hrov, xgamma, xlambda, xkappa, xN 
      !  output:	yield
      !-----------------------------------------------------------------------
      implicit none

      real(kind = dp) :: p, j, pp, theta, xphics, g_theta
      real(kind = dp) :: m_Hrov, xgamma, xlambda, xkappa, xN, beta
      real(kind = dp) :: yield_Hvor_log
      beta = 0.75 

      g_theta = cos(theta) + ( (sin(theta) * sin(xphics)) / sqrt(3.0d0) )
      g_theta = sin(xphics) / g_theta
      yield_Hvor_log = log( (j*sqrt(3.0)) / (g_theta) ) - &
                     ((beta + ((xkappa/xlambda)*(1.0-beta)))*log(p)) - &
                     ((1.0-(xkappa/xlambda))*(1.0-beta)*log(pp/2.0))

      !write(10, '(A, F10.5)') "Value of yield for F calculation: ", yield

   end function yield_Hvor_log


   !-----------------------------------------------------------------------
   ! Subroutine computing the derivatives
   !-----------------------------------------------------------------------
   subroutine derivatives_MCC(sig,p,j,xphics,theta,pp,dgdp,dfdsig,dgdsig)
      !-----------------------------------------------------------------------
      !	input:	sig,p,j,xphics,theta
      !	output:	dfdsig,dgdsig
      !-----------------------------------------------------------------------
      implicit none
      real(kind = dp), intent(in)  :: p, j, xphics, theta, pp
      real(kind = dp), intent(out), dimension(6) :: dfdsig,dgdsig
      real(kind = dp), intent(out) :: dgdp

      ! Local Variables
      real(kind = dp), dimension(6)  :: sig
      real(kind = dp) :: g_theta,dfdp,dfdj,dfdtheta
      real(kind = dp) :: dgdj,dgdtheta,dets,theta2
      real(kind = dp) :: dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
      real(kind = dp), dimension(6)  :: dpdsig,djdsig,d_dets_dsig,dthetadsig
      integer(kind = i32) :: i

      
      do i = 1, 3
         sig(i) = -sig(i)
      end do 

      g_theta = cos(theta) + ((sin(theta) * sin(xphics)) / sqrt(3.0d0))
      g_theta = sin(xphics) / g_theta
      write(10, '(A, F10.5)') "Value of g_theta in derivative fn:: ", g_theta

      dfdp = (1/p) * ( 1 - ( ( j / (p*g_theta) )**2) )
      !dfdp = (1/p) * ( (pp/p) - (2.0d0 * ( j / (p*g_theta) )**2) )
      !dfdp = -(((g_theta)**2)*pp) + (2.0d0*((g_theta)**2)*p)
      dfdj = (2*j) / ((p*g_theta)**2)
      !dfdj = (2*j)
      dfdtheta = (2*(j**2)) / (sqrt(3.0d0)*(p**2)*g_theta*(sin(xphics)))
      dfdtheta = dfdtheta * ( (cos(theta)*sin(xphics)) - sin(theta))

      write(10, '(A, F10.5)') "Value of dfdp in derivative fn:: ", dfdp
      write(10, '(A, F10.5)') "Value of dfdj in derivative fn:: ", dfdj
      write(10, '(A, E15.5)') "Value of dfdtheta in derivative fn:: ", dfdtheta

      dgdp = (1/p) * ( 1 - ( ( j / (p*g_theta) )**2) )
      !dgdp = (1/p) * ( (pp/p) - (2.0d0 * ( j / (p*g_theta) )**2) )
      !dgdp = -(((g_theta)**2)*pp) + (2.0d0*((g_theta)**2)*p)
      dgdj = (2*j) / ((p*g_theta)**2)
      !dgdj = (2*j)
      dgdtheta = 0

      write(10, '(A, F10.5)') "Value of dgdp in derivative fn:: ", dgdp
      write(10, '(A, F10.5)') "Value of dgdj in derivative fn:: ", dgdj
      write(10, '(A, F10.5)') "Value of dgdtheta in derivative fn:: ", dgdtheta


      dpdsig = 0.33333333d0*[1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0]
      djdsig = (1/(2*j))*[sig(1)-p, sig(2)-p, sig(3)-p, 2*sig(4), &
         2*sig(5), 2*sig(6)]

      write(10, '(A, 6F10.5)') "Value of djdsig in derivative fn:: ", djdsig   

      dummy1 = (2*sig(1))-sig(2)-sig(3)   
      dummy2 = (2*sig(2))-sig(1)-sig(3)  
      dummy3 = (2*sig(3))-sig(1)-sig(2)  
      dummy4 = (2*(sig(4)**2) -(sig(5)**2) -(sig(6)**2))
      dummy5 = (2*(sig(5)**2) -(sig(4)**2) -(sig(6)**2))
      dummy6 = (2*(sig(6)**2) -(sig(5)**2) -(sig(4)**2))

      write(10, '(A, E15.5)') "Value of dummy1 in derivative fn:: ", dummy1
      write(10, '(A, E15.5)') "Value of dummy2 in derivative fn:: ", dummy2
      write(10, '(A, E15.5)') "Value of dummy3 in derivative fn:: ", dummy3

      d_dets_dsig = [ ((2.0/27.0)*dummy2*dummy3) + ((1.0/27.0)*(dummy1**2)) - (1.0/3.0)*(dummy5), &
      ((2.0/27.0)*dummy3*dummy1) + ((1.0/27.0)*(dummy2**2)) - (1.0/3.0)*(dummy6), &
      ((2.0/27.0)*dummy1*dummy2) + ((1.0/27.0)*(dummy3**2)) - (1.0/3.0)*(dummy4), &
      (-(2.0/3.0)*sig(4)*dummy3) + (sig(5)*sig(6)), &
      (-(2.0/3.0)*sig(5)*dummy1) + (sig(4)*sig(6)),&
      (-(2.0/3.0)*sig(6)*dummy2) + (sig(5)*sig(4)) ]


      write(10, '(A, 6E15.5)') "Value of d_dets_dsig in derivative fn:: ", d_dets_dsig

      dets = ((sig(1)-p)*(sig(2)-p)*(sig(3)-p)) - ((sig(1)-p)*(sig(5)**2)) 
      dets =  dets  - ((sig(2)-p)*(sig(6)**2)) - ((sig(3)-p)*(sig(4)**2)) 
      dets = dets   + (2.0*sig(4)*sig(5)*sig(6))
      write(10, '(A, E15.5)') "Value of dets in derivative fn:: ", dets
      
      dthetadsig = (((dets/j)*djdsig) ) 
      dthetadsig = 4*(((dets/j)*djdsig) ) 

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      ! dthetadsig = (dthetadsig - d_dets_dsig) 
      dthetadsig = ( d_dets_dsig - dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      !dthetadsig = dthetadsig * (sqrt(3.0d0)/2.0)
      dthetadsig = dthetadsig * (2.0/sqrt(3.0d0)) * (j**4)

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      !write(10, '(A, E15.5)') "Value of dthetadsig in derivative fn:: ", cos(3.0*theta2)
      !write(10, '(A, E15.5)') "Value of dthetadsig in derivative fn:: ", theta2
      !write(10, '(A, E15.5)') "Value of dthetadsig in derivative fn:: ", 3.0*theta2
      !write(10, '(A, E15.5)') "Value of dthetadsig in derivative fn:: ", cos(3.0*theta2)*(j**3)


      !dthetadsig = dthetadsig / (cos(3.0*theta2)*(j**3))
      dthetadsig = dthetadsig / ((3.0*(j**8)) + (2.0*(dets**2)))
      

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      dfdsig = (dfdp * dpdsig) + (dfdj * djdsig) + (dfdtheta * dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dfdsig in derivative fn:: ", dfdsig

      dgdsig = (dgdp * dpdsig) + (dgdj * djdsig) + (dgdtheta * dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dgdsig in derivative fn:: ", dgdsig

      do i = 1, 3
         sig(i) = -sig(i)
      end do 

   end subroutine derivatives_MCC


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------

   !-----------------------------------------------------------------------
   ! Subroutine computing the derivatives
   !-----------------------------------------------------------------------
   subroutine derivatives_MCC_log(sig,p,j,xphics,theta,pp,dgdp,dfdsig,dgdsig)
      !-----------------------------------------------------------------------
      !	input:	sig,p,j,xphics,theta
      !	output:	dfdsig,dgdsig
      !-----------------------------------------------------------------------
      implicit none
      real(kind = dp), intent(in)  :: p, j, xphics, theta, pp
      real(kind = dp), intent(out), dimension(6) :: dfdsig,dgdsig
      real(kind = dp), intent(out) :: dgdp

      ! Local Variables
      real(kind = dp), dimension(6)  :: sig
      real(kind = dp) :: g_theta,dfdp,dfdj,dfdtheta
      real(kind = dp) :: dgdj,dgdtheta,dets,theta2
      real(kind = dp) :: dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
      real(kind = dp), dimension(6)  :: dpdsig,djdsig,d_dets_dsig,dthetadsig
      integer(kind = i32) :: i

      
      do i = 1, 3
         sig(i) = -sig(i)
      end do 

      g_theta = cos(theta) + ((sin(theta) * sin(xphics)) / sqrt(3.0d0))
      g_theta = sin(xphics) / g_theta
      write(10, '(A, F10.5)') "Value of g_theta in derivative fn:: ", g_theta

      dfdp = ((g_theta)**2)*((2.0*p)-pp)
      dfdj = (2*j)/3
      dfdtheta = (2*g_theta*p)*(-pp+p)*(g_theta**2)
      dfdtheta = dfdtheta * ( ((cos(theta)*sin(xphics))/sqrt(3.0d0)) - sin(theta))

      write(10, '(A, F10.5)') "Value of dfdp in derivative fn:: ", dfdp
      write(10, '(A, F10.5)') "Value of dfdj in derivative fn:: ", dfdj
      write(10, '(A, E15.5)') "Value of dfdtheta in derivative fn:: ", dfdtheta

      dgdp = ((g_theta)**2)*((2.0*p)-pp)
      dgdj = (2*j)/3
      dgdtheta = 0

      write(10, '(A, F10.5)') "Value of dgdp in derivative fn:: ", dgdp
      write(10, '(A, F10.5)') "Value of dgdj in derivative fn:: ", dgdj
      write(10, '(A, F10.5)') "Value of dgdtheta in derivative fn:: ", dgdtheta


      dpdsig = 0.33333333d0*[1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0]
      djdsig = (1/(2*j))*[sig(1)-p, sig(2)-p, sig(3)-p, 2*sig(4), &
         2*sig(5), 2*sig(6)]

      write(10, '(A, 6F10.5)') "Value of djdsig in derivative fn:: ", djdsig   

      dummy1 = (2*sig(1))-sig(2)-sig(3)   
      dummy2 = (2*sig(2))-sig(1)-sig(3)  
      dummy3 = (2*sig(3))-sig(1)-sig(2)  
      dummy4 = (2*(sig(4)**2) -(sig(5)**2) -(sig(6)**2))
      dummy5 = (2*(sig(5)**2) -(sig(4)**2) -(sig(6)**2))
      dummy6 = (2*(sig(6)**2) -(sig(5)**2) -(sig(4)**2))

      write(10, '(A, E15.5)') "Value of dummy1 in derivative fn:: ", dummy1
      write(10, '(A, E15.5)') "Value of dummy2 in derivative fn:: ", dummy2
      write(10, '(A, E15.5)') "Value of dummy3 in derivative fn:: ", dummy3

      d_dets_dsig = [ ((2.0/27.0)*dummy2*dummy3) + ((1.0/27.0)*(dummy1**2)) - (1.0/3.0)*(dummy5), &
      ((2.0/27.0)*dummy3*dummy1) + ((1.0/27.0)*(dummy2**2)) - (1.0/3.0)*(dummy6), &
      ((2.0/27.0)*dummy1*dummy2) + ((1.0/27.0)*(dummy3**2)) - (1.0/3.0)*(dummy4), &
      (-(2.0/3.0)*sig(4)*dummy3) + (sig(5)*sig(6)), &
      (-(2.0/3.0)*sig(5)*dummy1) + (sig(4)*sig(6)),&
      (-(2.0/3.0)*sig(6)*dummy2) + (sig(5)*sig(4)) ]


      write(10, '(A, 6E15.5)') "Value of d_dets_dsig in derivative fn:: ", d_dets_dsig

      dets = ((sig(1)-p)*(sig(2)-p)*(sig(3)-p)) - ((sig(1)-p)*(sig(5)**2)) 
      dets =  dets  - ((sig(2)-p)*(sig(6)**2)) - ((sig(3)-p)*(sig(4)**2)) 
      dets = dets   + (2.0*sig(4)*sig(5)*sig(6))
      write(10, '(A, E15.5)') "Value of dets in derivative fn:: ", dets
      
      dthetadsig = (((dets/j)*djdsig) ) 
      dthetadsig = 4*(((dets/j)*djdsig) ) 

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      ! dthetadsig = (dthetadsig - d_dets_dsig) 
      dthetadsig = ( d_dets_dsig - dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      !dthetadsig = dthetadsig * (sqrt(3.0d0)/2.0)
      dthetadsig = dthetadsig * (2.0/sqrt(3.0d0)) * (j**4)

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      !write(10, '(A, E15.5)') "Value of dthetadsig in derivative fn:: ", cos(3.0*theta2)
      !write(10, '(A, E15.5)') "Value of dthetadsig in derivative fn:: ", theta2
      !write(10, '(A, E15.5)') "Value of dthetadsig in derivative fn:: ", 3.0*theta2
      !write(10, '(A, E15.5)') "Value of dthetadsig in derivative fn:: ", cos(3.0*theta2)*(j**3)


      !dthetadsig = dthetadsig / (cos(3.0*theta2)*(j**3))
      dthetadsig = dthetadsig / ((3.0*(j**8)) + (2.0*(dets**2)))
      

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      dfdsig = (dfdp * dpdsig) + (dfdj * djdsig) + (dfdtheta * dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dfdsig in derivative fn:: ", dfdsig

      dgdsig = (dgdp * dpdsig) + (dgdj * djdsig) + (dgdtheta * dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dgdsig in derivative fn:: ", dgdsig

      do i = 1, 3
         sig(i) = -sig(i)
      end do 

   end subroutine derivatives_MCC_log


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------

   !-----------------------------------------------------------------------
   ! Subroutine computing the derivatives
   !-----------------------------------------------------------------------
   subroutine derivatives_Hvorslev(sig,p,j,xphics,theta,pp,m_Hvor,xgamma, &
      xlambda,xkappa,xN,dgdp,dfdsig,dgdsig)
      !-----------------------------------------------------------------------
      !	input:	sig,p,j,xphics,theta,pp,m_Hvor,xgamma,xlambda,xkappa,xN
      !	output:	dgdp,dfdsig,dgdsig
      !-----------------------------------------------------------------------
      implicit none
      real(kind = dp), intent(in)  :: p, j, xphics, theta, pp
      real(kind = dp), intent(in)  :: m_Hvor,xgamma,xlambda,xkappa,xN
      real(kind = dp), intent(out), dimension(6) :: dfdsig,dgdsig
      real(kind = dp), intent(out) :: dgdp

      ! Local Variables
      real(kind = dp), dimension(6)  :: sig
      real(kind = dp) :: g_theta,dfdp,dfdj,dfdtheta
      real(kind = dp) :: dgdj,dgdtheta,dets,theta2
      real(kind = dp) :: dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
      real(kind = dp), dimension(6)  :: dpdsig,djdsig,d_dets_dsig,dthetadsig
      integer(kind = i32) :: i

      
      do i = 1, 3
         sig(i) = -sig(i)
      end do 

      g_theta = cos(theta) + ((sin(theta) * sin(xphics)) / sqrt(3.0d0))
      g_theta = sin(xphics) / g_theta
      write(10, '(A, F10.5)') "Value of g_theta in derivative fn:: ", g_theta

      dfdp = -(m_Hvor/(g_theta-m_Hvor)) - &
      ((xkappa/xlambda)*(1/p)*(exp((xgamma - xN)/xlambda))*(exp((-xkappa/xlambda)*(log(pp/p)))))

      dfdj = (1/(g_theta-m_Hvor))

      dfdtheta = ((j-m_Hvor)/((g_theta-m_Hvor)**2))*(g_theta**2)*(1/sin(xphics))
      dfdtheta = dfdtheta * ( ((cos(theta)*sin(xphics))/sqrt(3.0)) - sin(theta))

      write(10, '(A, F10.5)') "Value of dfdp in derivative fn:: ", dfdp
      write(10, '(A, F10.5)') "Value of dfdj in derivative fn:: ", dfdj
      write(10, '(A, E15.5)') "Value of dfdtheta in derivative fn:: ", dfdtheta

      dgdp = (1/p) * ( 1 - ( ( j / (p*g_theta) )**2) )
      dgdj = (2*j) / ((p*g_theta)**2)
      dgdtheta = 0

      write(10, '(A, F10.5)') "Value of dgdp in derivative fn:: ", dgdp
      write(10, '(A, F10.5)') "Value of dgdj in derivative fn:: ", dgdj
      write(10, '(A, F10.5)') "Value of dgdtheta in derivative fn:: ", dgdtheta


      dpdsig = 0.33333333d0*[1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0]
      djdsig = (1/(2*j))*[sig(1)-p, sig(2)-p, sig(3)-p, 2*sig(4), &
         2*sig(5), 2*sig(6)]

      write(10, '(A, 6F10.5)') "Value of djdsig in derivative fn:: ", djdsig   

      dummy1 = (2*sig(1))-sig(2)-sig(3)   
      dummy2 = (2*sig(2))-sig(1)-sig(3)  
      dummy3 = (2*sig(3))-sig(1)-sig(2)  
      dummy4 = (2*(sig(4)**2) -(sig(5)**2) -(sig(6)**2))
      dummy5 = (2*(sig(5)**2) -(sig(4)**2) -(sig(6)**2))
      dummy6 = (2*(sig(6)**2) -(sig(5)**2) -(sig(4)**2))

      write(10, '(A, E15.5)') "Value of dummy1 in derivative fn:: ", dummy1
      write(10, '(A, E15.5)') "Value of dummy2 in derivative fn:: ", dummy2
      write(10, '(A, E15.5)') "Value of dummy3 in derivative fn:: ", dummy3

      d_dets_dsig = [ ((2.0/27.0)*dummy2*dummy3) + ((1.0/27.0)*(dummy1**2)) - (1.0/3.0)*(dummy5), &
      ((2.0/27.0)*dummy3*dummy1) + ((1.0/27.0)*(dummy2**2)) - (1.0/3.0)*(dummy6), &
      ((2.0/27.0)*dummy1*dummy2) + ((1.0/27.0)*(dummy3**2)) - (1.0/3.0)*(dummy4), &
      (-(2.0/3.0)*sig(4)*dummy3) + (sig(5)*sig(6)), &
      (-(2.0/3.0)*sig(5)*dummy1) + (sig(4)*sig(6)),&
      (-(2.0/3.0)*sig(6)*dummy2) + (sig(5)*sig(4)) ]


      write(10, '(A, 6E15.5)') "Value of d_dets_dsig in derivative fn:: ", d_dets_dsig

      dets = ((sig(1)-p)*(sig(2)-p)*(sig(3)-p)) - ((sig(1)-p)*(sig(5)**2)) 
      dets =  dets  - ((sig(2)-p)*(sig(6)**2)) - ((sig(3)-p)*(sig(4)**2)) 
      dets = dets   + (2.0*sig(4)*sig(5)*sig(6))
      write(10, '(A, E15.5)') "Value of dets in derivative fn:: ", dets
      
      dthetadsig = (((dets/j)*djdsig) ) 
      dthetadsig = 4*(((dets/j)*djdsig) ) 

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      ! dthetadsig = (dthetadsig - d_dets_dsig) 
      dthetadsig = ( d_dets_dsig - dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      !dthetadsig = dthetadsig * (sqrt(3.0d0)/2.0)
      dthetadsig = dthetadsig * (2.0/sqrt(3.0d0)) * (j**4)

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      !dthetadsig = dthetadsig / (cos(3.0*theta2)*(j**3))
      dthetadsig = dthetadsig / ((3.0*(j**8)) + (2.0*(dets**2)))
      

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      dfdsig = (dfdp * dpdsig) + (dfdj * djdsig) + (dfdtheta * dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dfdsig in derivative fn:: ", dfdsig

      dgdsig = (dgdp * dpdsig) + (dgdj * djdsig) + (dgdtheta * dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dgdsig in derivative fn:: ", dgdsig

      do i = 1, 3
         sig(i) = -sig(i)
      end do 

   end subroutine derivatives_Hvorslev
   

   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------

   !-----------------------------------------------------------------------
   ! Subroutine computing the derivatives
   !-----------------------------------------------------------------------
   subroutine derivatives_Hvorslev_log(sig,p,j,xphics,theta,pp,m_Hvor,xgamma, &
      xlambda,xkappa,xN,dgdp,dfdsig,dgdsig)
      !-----------------------------------------------------------------------
      !	input:	sig,p,j,xphics,theta,pp,m_Hvor,xgamma,xlambda,xkappa,xN
      !	output:	dgdp,dfdsig,dgdsig
      !-----------------------------------------------------------------------
      implicit none
      real(kind = dp), intent(in)  :: p, j, xphics, theta, pp
      real(kind = dp), intent(in)  :: m_Hvor,xgamma,xlambda,xkappa,xN
      real(kind = dp), intent(out), dimension(6) :: dfdsig,dgdsig
      real(kind = dp), intent(out) :: dgdp

      ! Local Variables
      real(kind = dp), dimension(6)  :: sig
      real(kind = dp) :: g_theta,dfdp,dfdj,dfdtheta,beta
      real(kind = dp) :: dgdj,dgdtheta,dets,theta2
      real(kind = dp) :: dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
      real(kind = dp), dimension(6)  :: dpdsig,djdsig,d_dets_dsig,dthetadsig
      integer(kind = i32) :: i

      
      do i = 1, 3
         sig(i) = -sig(i)
      end do 

      g_theta = cos(theta) + ((sin(theta) * sin(xphics)) / sqrt(3.0d0))
      g_theta = sin(xphics) / g_theta
      write(10, '(A, F10.5)') "Value of g_theta in derivative fn:: ", g_theta

      beta = 0.75 

      dfdp = -(beta + ((xkappa/xlambda)*(1-beta))  ) * (1/p)   

      dfdj = (1/(j))

      dfdtheta = -g_theta*(1/sin(xphics))
      dfdtheta = dfdtheta * ( ((cos(theta)*sin(xphics))/sqrt(3.0)) - sin(theta))

      write(10, '(A, F10.5)') "Value of dfdp in derivative fn:: ", dfdp
      write(10, '(A, F10.5)') "Value of dfdj in derivative fn:: ", dfdj
      write(10, '(A, E15.5)') "Value of dfdtheta in derivative fn:: ", dfdtheta

      
      dgdp = ((g_theta)**2)*((2.0*p)-pp)
      dgdj = (2*j)/3
      dgdtheta = 0

      write(10, '(A, F10.5)') "Value of dgdp in derivative fn:: ", dgdp
      write(10, '(A, F10.5)') "Value of dgdj in derivative fn:: ", dgdj
      write(10, '(A, F10.5)') "Value of dgdtheta in derivative fn:: ", dgdtheta


      dpdsig = 0.33333333d0*[1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0]
      djdsig = (1/(2*j))*[sig(1)-p, sig(2)-p, sig(3)-p, 2*sig(4), &
         2*sig(5), 2*sig(6)]

      write(10, '(A, 6F10.5)') "Value of djdsig in derivative fn:: ", djdsig   

      dummy1 = (2*sig(1))-sig(2)-sig(3)   
      dummy2 = (2*sig(2))-sig(1)-sig(3)  
      dummy3 = (2*sig(3))-sig(1)-sig(2)  
      dummy4 = (2*(sig(4)**2) -(sig(5)**2) -(sig(6)**2))
      dummy5 = (2*(sig(5)**2) -(sig(4)**2) -(sig(6)**2))
      dummy6 = (2*(sig(6)**2) -(sig(5)**2) -(sig(4)**2))

      write(10, '(A, E15.5)') "Value of dummy1 in derivative fn:: ", dummy1
      write(10, '(A, E15.5)') "Value of dummy2 in derivative fn:: ", dummy2
      write(10, '(A, E15.5)') "Value of dummy3 in derivative fn:: ", dummy3

      d_dets_dsig = [ ((2.0/27.0)*dummy2*dummy3) + ((1.0/27.0)*(dummy1**2)) - (1.0/3.0)*(dummy5), &
      ((2.0/27.0)*dummy3*dummy1) + ((1.0/27.0)*(dummy2**2)) - (1.0/3.0)*(dummy6), &
      ((2.0/27.0)*dummy1*dummy2) + ((1.0/27.0)*(dummy3**2)) - (1.0/3.0)*(dummy4), &
      (-(2.0/3.0)*sig(4)*dummy3) + (sig(5)*sig(6)), &
      (-(2.0/3.0)*sig(5)*dummy1) + (sig(4)*sig(6)),&
      (-(2.0/3.0)*sig(6)*dummy2) + (sig(5)*sig(4)) ]


      write(10, '(A, 6E15.5)') "Value of d_dets_dsig in derivative fn:: ", d_dets_dsig

      dets = ((sig(1)-p)*(sig(2)-p)*(sig(3)-p)) - ((sig(1)-p)*(sig(5)**2)) 
      dets =  dets  - ((sig(2)-p)*(sig(6)**2)) - ((sig(3)-p)*(sig(4)**2)) 
      dets = dets   + (2.0*sig(4)*sig(5)*sig(6))
      write(10, '(A, E15.5)') "Value of dets in derivative fn:: ", dets
      
      dthetadsig = (((dets/j)*djdsig) ) 
      dthetadsig = 4*(((dets/j)*djdsig) ) 

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      ! dthetadsig = (dthetadsig - d_dets_dsig) 
      dthetadsig = ( d_dets_dsig - dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      !dthetadsig = dthetadsig * (sqrt(3.0d0)/2.0)
      dthetadsig = dthetadsig * (2.0/sqrt(3.0d0)) * (j**4)

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      !dthetadsig = dthetadsig / (cos(3.0*theta2)*(j**3))
      dthetadsig = dthetadsig / ((3.0*(j**8)) + (2.0*(dets**2)))
      

      write(10, '(A, 6E15.5)') "Value of dthetadsig in derivative fn:: ", dthetadsig

      dfdsig = (dfdp * dpdsig) + (dfdj * djdsig) + (dfdtheta * dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dfdsig in derivative fn:: ", dfdsig

      dgdsig = (dgdp * dpdsig) + (dgdj * djdsig) + (dgdtheta * dthetadsig)

      write(10, '(A, 6E15.5)') "Value of dgdsig in derivative fn:: ", dgdsig

      do i = 1, 3
         sig(i) = -sig(i)
      end do 

   end subroutine derivatives_Hvorslev_log

   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   

   !-----------------------------------------------------------------------
   ! Subroutine computing the material stiffness parameters
   !-----------------------------------------------------------------------
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
      r = 3. * ( 1. - 2.*xNu) / ( 2. * (1.+xNu))
      xG = r*xK !shear modulus

   end subroutine stiffnessMCC


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------



   !-------------------------------------------------------------------
   Subroutine MatVec(xMat,IM,Vec,N,VecR)
      !-------------------------------------------------------------------
      !     Calculate VecR = xMat*Vec
      !
      ! I   xMat  : (Square) Matrix (IM,*)
      ! I   Vec   : Vector
      ! I   N     : Number of rows/colums
      ! O   VecR  : Resulting vector
      !
      !-------------------------------------------------------------------
      Implicit none
      real(kind = dp) :: xMat(:, :), Vec(:), VecR(:),X
      integer(kind = i32) :: I,J,N,IM
      ! Dimension xMat(IM,*),Vec(*),VecR(*)
      !-------------------------------------------------------------------
      Do I=1,N
         X=0
         Do J=1,N
            X=X+xMat(I,J)*Vec(J)
         End Do
         VecR(I)=X
      End Do
      Return
   End Subroutine MatVec


   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------


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

end module mod_UMAT_Cam_Clay




