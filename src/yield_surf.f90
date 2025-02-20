module mod_yield_surface

   use kind_precision_module, only: dp, i32

contains

   ! Function to calculate the yield surface
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


   ! Derivatives of the yield surface w/ respect to the invariants

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

   
end module mod_yield_surface
