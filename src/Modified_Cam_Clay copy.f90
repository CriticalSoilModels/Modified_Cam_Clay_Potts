module mod_UMAT_Cam_Clay
   use kind_precision_module, only: dp, i32
   use mod_implicit_stress_integ, only: implicit_predictor_corrector_integration


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
      open(unit=10, file=output_file, status='replace', action='write', iostat=status)
      

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



end module mod_UMAT_Cam_Clay




