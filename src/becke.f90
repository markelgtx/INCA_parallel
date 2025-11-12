!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine becke(rrg,sumq,rgrid,nquad,cent,w_beck,Ps)                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Based on A. D. Becke's paper from 1987, A multicenter numerical integration scheme for   !
! polyatomic molecules.                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes Becke's weights for each (sym reduced) point of the grid, with a weight for each!
! center. Then normalizes this weight so that the sum of all the weights of a single point !
! is 1. Finally stores the weight of the corresponding grid point in the array w_beck,     !
! since we will only use the weight of the quadrature centre from where the curren pòints  !
! have been originated.                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! As input: grid points, number of points per quadrature center, number of quad. centers   !
! As output: becke weights for each grid point.                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 implicit none
 integer, intent(in) :: rgrid, nquad !total grid, number of quadratures
 double precision, intent(in), dimension(3,nquad) :: cent !center of the quadratures
 double precision, intent(in), dimension(3,rgrid) :: rrg  !(reduced) grid points
 double precision, intent(in), dimension(nquad) :: Ps !weight of each center
 integer, dimension(nquad) :: sumq            !number of grid points per quad.
 double precision, intent(out), dimension(rgrid) :: w_beck !becke weights for each point
 !local variables
 double precision, dimension(rgrid,nquad) :: w_becke !weights of different quad. centers at same point
 integer :: i1,i,j,sm
 double precision, dimension(nquad,nquad) :: Rij, Xi !distance between centers, equation A4 xi
 double precision, dimension(nquad) :: P             !equation 13
 double precision :: ri, rj, s_ij, Ptot
 double precision ::  mu_ij, v_ij, a_ij, u_ij !eq.11, eq.a2, eq.a5, eq.a6 

 !compute R_ij
 do i=1,nquad
  do j=1,nquad
     Rij(i,j)=dsqrt(sum((cent(:,i)-cent(:,j))**2.d0))  !distance between centers
     Xi(i,j)=Ps(i)/Ps(j)                              !relacio entre centres
  end do
 end do     
 !compute Becke weights from all centres for all the grid points
 w_becke=0.d0     
 sm=0     
 do i1=1,rgrid !for all the grid points
   sm=sm+1     !count grid points
   do i=1,nquad   
        P(i)=1.d0
        do j=1,nquad
           if (i.ne.j) then        
             ri=dsqrt(sum((rrg(:,i1)-cent(:,i))**2.d0)) !distance to center i (from point r)
             rj=dsqrt(sum((rrg(:,i1)-cent(:,j))**2.d0)) !distance to center j (from point r)
             mu_ij=(ri-rj)*(Rij(i,j)**(-1.d0))
             u_ij=(xi(i,j)-1.d0)/(xi(i,j)+1.d0)  
             a_ij=u_ij*((u_ij**2.d0)-1.d0)**(-1.d0)
             v_ij=mu_ij+a_ij*(1.d0-mu_ij**2.d0) 
             s_ij=0.5d0*(1.d0-f_k(v_ij))
             P(i)=P(i)*s_ij
           end if              
        end do
       w_becke(i1,i)=P(i) !store weight of quadrature i at point i1     
   end do 
   Ptot=sum(P)
   w_becke(i1,:)=w_becke(i1,:)*(Ptot**(-1.d0)) !normalize weight to fulfill equation 3   
 end do
 
 !store single becke weight for each grid point in a single array
 sm=0
 w_beck=0.d0
 do i=1,nquad
  do j=1,sumq(i)
    sm=sm+1   
    w_beck(sm)=w_becke(sm,i)
  end do       
 end do

 contains        
     function f_k(val)   !equation 20, with value k=3
     implicit none
     double precision :: f_k
     double precision, intent(in) :: val
     double precision :: vl
     integer :: i
     integer, parameter :: k=3
     vl=val
         do i=1,k         
            f_k=pf(vl)
            vl=f_k
         end do       
     end function
     function pf(vl)  !equation 19
     implicit none
     double precision :: pf
     double precision, intent(in) :: vl
         pf=1.5d0*vl-0.5d0*vl**(3.d0)     
     end function
end subroutine becke         

