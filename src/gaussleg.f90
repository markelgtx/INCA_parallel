     





!******************************************************************************

!Gauss-Legendreren subroutine. t eta w balioak ematen diguna

!Modification of 
!https://github.com/NREL/OpenWARP/blob/master/source/NemohImproved/Nemoh/Solver/Core/Gaussm3.f90
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays t(1:n) and  w(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula. 
!* For detailed explanations finding weights & abscissas, see  "Numerical Recipes in Fortran */

      SUBROUTINE sub_GauLeg(X1,X2,t,w,n)
	!beste modu bat X1 eta X2 kendu, eta XL eta XM
	!bukaera t(i) = -z eta t(n+1-i)= +z izango dira

      IMPLICIT NONE
      
      INTEGER :: m,i,j
      INTEGER,intent(in) :: n   !Number of Gaussian points
      REAL*8, dimension(n),intent(out) :: W,t
      REAL*8 :: XM,XL,X1,X2,EPS,P1,P2,P3,pi,Z1,Z,PP

	!Relative precision
   	EPS = 1.D-14

  	!double precision arccosine. Pi value=3.14159
 	pi = DACOS(-1.D0)

	!N = number of Gauss Points
	!Roots are symmetric in the interval - so only need to find half of them  
   	m = int((dble(n) + 1.0d0) / 2.0d0)
	
	!The coats are going to be X1 = -1 and X2 = 1, Gauss-Legendre 
      	XM=0.5D0*(X1+X2)
      	XL=0.5D0*(X2-X1)


	!Loop over the desired roots
      	DO i = 1,m
         Z = DCOS (pi * (dble(i) - 0.25D0)/(dble(n) + 0.5D0))
	!Starting with the above approximation to the i-th root,
	!we enter the main loop of refinement by NEWTON'S method   
 10      P1 = 1.D0
         P2 = 0.D0

	!Loop up the recurrence relation to get the Legendre
	!polynomial evaluated at z                
         DO j = 1,n
            P3 = P2
            P2 = P1
            P1 = ((2.D0 * dble(j) - 1.D0) * Z * P2 - (dble(j)-1.D0)*P3)&
                /dble(j)
         END DO
!p1 is now the desired Legendre polynomial.
!We next compute pp, its derivative, by a standard relation involving also p2, 
!the polynomial of one lower order. 
         PP = dble(n) * (Z * P1 - P2)/(Z * Z - 1.D0)
         Z1 = Z
         Z = Z1 - P1/PP	      ! Newton's Method  */

         IF (DABS(Z-Z1) .GT. EPS) GO TO 10

	! Roots will be symmetric about the origin  
         t(i) = XM - XL * Z
         t(n + 1 - i) = XM + XL * Z
	!Compute the weight and its symmetric counterpart 
         W(i) = 2.D0 * XL/((1.D0 - Z * Z) * PP * PP)
         W(n + 1 - i) = W(i)
      END DO  

	RETURN   !not neccesary
      END SUBROUTINE Sub_GauLeg


