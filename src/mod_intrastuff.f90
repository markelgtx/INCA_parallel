module intrastuff 
use geninfo
use quadratures
implicit none
integer :: i,j,k,l !primitive quartets
double precision :: sqe !dsqrt(e_ijkl) save time and compute it once
!Grid independent variables for the intracule:
double precision :: a_ik, a_jl, a_ijkl 
double precision :: e_ik, e_jl, e_ijkl 
double precision, dimension(3) :: R_i, R_k, R_j, R_l !X_i, ...
double precision, dimension(3) :: R_ik, R_jl, R_ijkl !X_ik, ...
double precision :: Alf_ijkl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
!!!!!!!!!!!!!Functions for the first integral screening!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        function J_ik(i,k,R_i_k_2) !J upper bound for the first integral screening
        implicit none
        double precision :: J_ik       
        integer, intent(in) :: i,k
        double precision, intent(in) :: R_i_k_2
        double precision, dimension(3) :: x_max
        double precision :: r_part
        integer :: n_nodes, ii
    
        do ii=1,3 !loop for J_ik(x_i(*), (y) and (z)
                n_nodes=TMN(i,ii)+TMN(k,ii) + 1          !obtain x_max
                if (n_nodes.eq.1) then
                        x_max(ii)=0.d0
                else if (n_nodes.eq.2) then
                        x_max(ii)=0.707106781186548d0 
                else if (n_nodes.eq.3) then 
                        x_max(ii)=1.224744871391589d0 
                else if (n_nodes.eq.4) then
                        x_max(ii)=1.650680123885785d0
                else if (n_nodes.eq.5) then 
                        x_max(ii)=2.020182870456086d0 
                else if (n_nodes.eq.6) then 
                        x_max(ii)=2.350604973674492d0
                else if (n_nodes.eq.7) then 
                        x_max(ii)=2.651961356835233d0 
                else if (n_nodes.eq.8) then 
                        x_max(ii)=2.930637420257244d0
                else if (n_nodes.eq.9) then 
                        x_max(ii)=3.190993201781528d0  
                else if (n_nodes.eq.10) then 
                        x_max(ii)=3.436159118837738d0
                end if
        end do
                
        J_ik=pi**(1.5d0)*(2.d0*a_ik)**(-(dble(TMN(i,1)+TMN(k,1)+TMN(i,2)+TMN(k,2)+TMN(i,3)+TMN(k,3))+1.5d0))&
        *dexp(-2.d0*e_ik*(R_i_k_2))
        r_part=1.d0
        do ii=1,3
              r_part=r_part &
              *((x_max(ii)+Alpha(k)*dsqrt(2.d0/a_ik)*dabs(Cartes(Ra(k),ii)-Cartes(Ra(i),ii)))**(2.d0*dble(TMN(i,ii))))&
              *((x_max(ii)+Alpha(i)*dsqrt(2.d0/a_ik)*dabs(Cartes(Ra(k),ii)-Cartes(Ra(i),ii)))**(2.d0*dble(TMN(k,ii))))
        end do
        J_ik=J_ik*r_part
        end function
 
        function J_jl(j,l,R_j_l_2) !J upper bound for the first integral screening
        double precision :: J_jl
        integer, intent(in) :: j,l
        double precision, intent(in) :: R_j_l_2
        double precision, dimension(3) :: x_max 
        integer :: n_nodes, ii
        n_nodes=0  
        do ii=1,3 !loop for J_ik(x_i(*), (y) and (z)
                n_nodes=TMN(j,ii)+TMN(l,ii) + 1
                if (n_nodes.eq.1) then 
                        x_max(ii)=0.d0
                else if (n_nodes.eq.2) then
                        x_max(ii)=0.707106781186548d0 
                else if (n_nodes.eq.3) then 
                        x_max(ii)=1.224744871391589d0 
                else if (n_nodes.eq.4) then 
                        x_max(ii)=1.650680123885785d0
                else if (n_nodes.eq.5) then 
                        x_max(ii)=2.020182870456086d0 
                else if (n_nodes.eq.6) then 
                        x_max(ii)=2.350604973674492d0
                else if (n_nodes.eq.7) then 
                        x_max(ii)=2.651961356835233d0 
                else if (n_nodes.eq.8) then 
                        x_max(ii)=2.930637420257244d0
                else if (n_nodes.eq.9) then 
                        x_max(ii)=3.190993201781528d0  
                else if (n_nodes.eq.10) then 
                        x_max(ii)=3.436159118837738d0
                end if
        end do
  
        J_jl=pi**(1.5d0)*(2.d0*a_jl)**(-(dble(TMN(j,1)+TMN(l,1)+TMN(j,2)+TMN(l,2)+TMN(j,3)+TMN(l,3))+1.5d0))&
        *dexp(-2.d0*e_jl*(R_j_l_2))&
        *((x_max(1)+Alpha(l)*dsqrt(2.d0/a_jl)*dabs(Cartes(Ra(l),1)-Cartes(Ra(j),1)))**(2.d0*dble(TMN(j,1))))&
        *((x_max(1)+Alpha(j)*dsqrt(2.d0/a_jl)*dabs(Cartes(Ra(l),1)-Cartes(Ra(j),1)))**(2.d0*dble(TMN(l,1))))&
        *((x_max(2)+Alpha(l)*dsqrt(2.d0/a_jl)*dabs(Cartes(Ra(l),2)-Cartes(Ra(j),2)))**(2.d0*dble(TMN(j,2))))&
        *((x_max(2)+Alpha(j)*dsqrt(2.d0/a_jl)*dabs(Cartes(Ra(l),2)-Cartes(Ra(j),2)))**(2.d0*dble(TMN(l,2))))&
        *((x_max(3)+Alpha(l)*dsqrt(2.d0/a_jl)*dabs(Cartes(Ra(l),3)-Cartes(Ra(j),3)))**(2.d0*dble(TMN(j,3))))&
        *((x_max(3)+Alpha(j)*dsqrt(2.d0/a_jl)*dabs(Cartes(Ra(l),3)-Cartes(Ra(j),3)))**(2.d0*dble(TMN(l,3))))
        end function
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

   !     function A_ind(R_i_k_2,R_j_l_2)  !grid independent part of the intracule
   !     implicit none
   !     double precision, intent(in) :: R_i_k_2, R_j_l_2
   !     double precision :: A_ind
   !     A_ind=(a_ijkl)**(-1.5d0)* dexp(-e_ik*R_i_k_2-e_jl*R_j_l_2) !eq.18
 ! 
 !       end function A_ind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!eq 15!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!


  function Wr(rhat,r,ax)  !for 1st integral screening
        double precision :: Wr
        double precision, intent(in) :: rhat,r
        integer, intent(in) :: ax
        Wr=((dsqrt(a_ijkl)**(-1.d0)*rhat)+(alf_ijkl-0.5d0)*r+(r_ijkl(ax)-Cartes(Ra(i),ax)))**dble(TMN(i,ax))*&
        ((dsqrt(a_ijkl)**(-1.d0)*rhat)+(alf_ijkl+0.5d0)*r+(r_ijkl(ax)-Cartes(Ra(j),ax)))**dble(TMN(j,ax))*&
        ((dsqrt(a_ijkl)**(-1.d0)*rhat)+(alf_ijkl-0.5d0)*r+(r_ijkl(ax)-Cartes(Ra(k),ax)))**dble(TMN(k,ax))*&
        ((dsqrt(a_ijkl)**(-1.d0)*rhat)+(alf_ijkl+0.5d0)*r+(r_ijkl(ax)-Cartes(Ra(l),ax)))**dble(TMN(l,ax))
  end function Wr
      

 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine polycoef(C_r,np,n,ax,rh,w_r,ipiv)
!Evaluates W function (eqn.15) at xx(np=L+1) points and sets up a system of linear equations (see
!matrix M). 
!Solves this linear equations using Lapack subroutine dgesv.
!As output, we get the polynomial coefficients(C_r(np)) of V_r. This process is repeated for 
!x (ax=1), y (ax=2) and z (ax=3).

!n is the number of nodes of the quadrature
!np is the number of unknowns
!w_n is a matrix with the n weights of the G-Hermite quadrature (in module quadratures)
!rh is a matrix with the n nodes of the G-Hermite quadrature (in module quadratures)
IMPLICIT NONE
!global variables
INTEGER, intent(in) :: np, n, ax !number of points, number of nodes, axis
double precision, intent(out) :: C_r(:)
double precision, intent(in) :: w_r(:), rh(:)
double precision, intent(inout) :: ipiv(:) !for dgbsv
!local variables
double precision, dimension(np,np) :: M  
double precision, dimension(np) :: xx    !Points where we evaluate W
integer :: Ltot
integer :: i1, j1, l1

!variables for dgbsv!!!!!
integer :: INFO

if (np.eq.1) then !only one coefficient
   C_r=w_r(1)     !see equation 15 (in this case we only have 1 node)
else
   Ltot=np-1 !degree of the polynomial    
   !generate 'np' points to evaluate the polynomial (xx(np)) !linearly independent!! 
   do i1=1,np
     xx(i1)=-0.5d0+dble(i1)
   end do   
   !build M matrix (M and C_r form the augmented matrix)
   do i1=1,np
       l1=Ltot !start from the maximum degree
       do j1=1,np
          M(i1,j1)=(sqe*(xx(i1)+r_ik(ax)-r_jl(ax)))**(dble(l1)) !r'^l               
          l1=l1-1
       end do
   end do
   C_r=0.d0
   do i1=1,np
       do j1=1,n        
         C_r(i1)=C_r(i1)+w_r(j1)*Wr(rh(j1),xx(i1),ax) !V_ijkl(X') (~eq.16)
       end do 
   end do 
   !!!Solve linear system with lapack!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call dgesv(np,1,M,np,ipiv,C_r,np,INFO)
   if (INFO.ne.0) write(*,*) "ERROR, cannot solve linear equations"
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if
END SUBROUTINE Polycoef

subroutine gauherm(Lrtot,n,rh,w_r) 
        !performs gauss-hermite quadrature                          
        !we will obtain the nodes and weights depending on the degree of the polinomial (2n-1) --> (n) 
        !(n is nx ny or nz in intracule.f90)
         implicit none
         !global variable
         integer, intent(in) :: Lrtot !degree of polynomial
         integer, intent(out) :: n    !number of gauss-hermite nodes
         double precision,intent(out) :: rh(:), w_r(:) !nodes
         !local variables
         integer :: factn, i2
         double precision :: Cons
         double precision, parameter :: pi=4.d0*datan(1.d0)
         if (Lrtot.gt.0) then
         !compute the number of nodes for exact Hermite quadrature
             if(MOD(Lrtot,2).eq.0) then !even 
                   n=int((dble(Lrtot)*0.5d0) +1.0d0) !nx->number of nodes
             else 
                   n=int((dble(Lrtot)+1.0d0)*0.5d0) 
             end if
         else
          !only one node
            n=1
         end if  
         !store rh values depending on pol. degree
          if (n.eq.1) then !1 node-->degree of the pol is 0 or 1 
           rh(1)=0.0000000000000d0  
          else if (n.eq.2) then !2 nodes-->degree of the pol is 3 or 2                
           rh(2)=0.7071067811865475d0 
           rh(1)=-0.7071067811865475d0 
          else if (n.eq.3) then !3 nodes-->degree of pol is 5 or 4
           rh(2)=0.d0
           rh(3)= 1.224744871391589d0 
           rh(1)=-1.224744871391589d0  
          else if (n.eq.4) then   !7 or 6
           rh(3)=0.5246476232752903d0 
           rh(4)=1.650680123885785d0 
           rh(2)=-0.5246476232752903d0 
           rh(1)=-1.650680123885785d0 
          else if (n.eq.5) then    !9 or 8
           rh(1)=-2.020182870456086d0
           rh(2)=-0.9585724646138185d0
           rh(3)=0.000000000000000d0
           rh(4)=2.020182870456086d0
           rh(5)=0.9585724646138185d0
          else if (n.eq.6) then    !11 or 10
           rh(4)=0.4360774119276165d0 
           rh(5)=1.335849074013697d0 
           rh(6)=2.350604973674492d0   
           rh(3)=-0.4360774119276165d0 
           rh(2)=-1.335849074013697d0 
           rh(1)=-2.350604973674492d0
          else if (n.eq.7) then   !13 or 12
           rh(4)=0.d0
           rh(5)=0.8162878828589647d0 
           rh(6)=1.673551628767471d0 
           rh(7)=2.651961356835233d0
           rh(3)=-0.8162878828589647d0 
           rh(2)=-1.673551628767471d0 
           rh(1)=-2.651961356835233d0
          else if (n.eq.8) then    !15 or 14
           rh(5)=0.3811869902073221d0 
           rh(6)=1.157193712446780d0 
           rh(7)=1.981656756695843d0 
           rh(8)=2.930637420257244d0
           rh(4)=-0.3811869902073221d0 
           rh(3)=-1.157193712446780d0 
           rh(2)=-1.981656756695843d0 
           rh(1)=-2.930637420257244d0
          else if (n.eq.9) then     !17 or 16
           rh(6)=0.7235510187528376d0 
           rh(7)=1.468553289216668d0 
           rh(8)=2.266580584531843d0 
           rh(9)=3.190993201781528d0 
           rh(5)=0.d0
           rh(4)=-0.7235510187528376d0 
           rh(3)=-1.468553289216668d0 
           rh(2)=-2.266580584531843d0 
           rh(1)=-3.190993201781528d0 
          else if (n.eq.10) then     !19 or 18
            rh(6)=0.3429013272237046d0 
            rh(7)=1.036610829789514d0 
            rh(8)=1.756683649299882d0 
            rh(9)=2.532731674232790d0 
            rh(10)=3.436159118837738d0 
            rh(5)=-0.3429013272237046d0 
            rh(4)=-1.036610829789514d0 
            rh(3)=-1.756683649299882d0 
            rh(2)=-2.532731674232790d0 
            rh(1)=-3.436159118837738d0 
          else if (n.eq.11) then     !degree of pol is 21 or 20
            rh(6)=0.d0
            rh(7)=0.6568095668820998d0 
            rh(8)=1.326557084494933d0 
            rh(9)=2.025948015825755d0 
            rh(10)=2.783290099781652d0 
            rh(11)=3.668470846559583d0 
            rh(5)=-0.6568095668820998d0 
            rh(4)=-1.326557084494933d0 
            rh(3)=-2.025948015825755d0 
            rh(2)=-2.783290099781652d0 
            rh(1)=-3.668470846559583d0 
          else 
           write(*,*) "Error, total ang. momenta greater than 3"
          end if
          
          factn=1     !compute weights using the equation
          do i2=1,n   
           factn=factn*i2  
          end do
          Cons=((2.d0)**(dble(n-1)) * dble(factn)* dsqrt(pi))/dble(n*n)
          do i2=1,n
           w_r(i2)=cons*(1.d0/(Hermite(rh(i2),n))**2.d0)  
          end do  

        
          contains
            function Hermite(x,n) !
            double precision :: Hermite
            integer, intent(in) :: n !numer of nodes of pol
            double precision, intent(in) :: x !roots of Hermite pol
           
           if (n.eq.1) then !0 deg
             Hermite=1.d0
           else if (n.eq.2) then
             Hermite=2.d0*x 
           else if (n.eq.3) then
             Hermite=4.d0*(x**2.d0) -2.d0
           else if (n.eq.4) then
             Hermite=8.d0*(x**3.d0) -12.d0*x
           else if (n.eq.5) then
             Hermite=16.d0*(x**4.d0)-48.d0*(x**2.d0)+12.d0  
           else if (n.eq.6) then
             Hermite=32.d0*(x**5.d0)-160.d0*(x**3.d0)+120.d0*x  
           else if (n.eq.7) then
             Hermite=64.d0*(x**6.d0)-480.d0*(x**4.d0)+720.d0*(x**2.d0)-120.d0
           else if (n.eq.8) then
             Hermite=128.d0*(x**7.d0)-1344.d0*(x**5.d0)+3360.d0*(x**3.d0)-1680.d0*x    
           else if (n.eq.9) then
             Hermite=256.d0*(x**8.d0)-3584.d0*(x**6.d0)+13440.d0*(x**4.d0)-13440.d0*(x**2.d0)+1680.d0  
           else if (n.eq.10) then
             Hermite=512.d0*(x**9.d0)-9216.d0*(x**7.d0)+48384.d0*(x**5.d0)-80640.d0*(x**3.d0)+30240.d0*x
           else if (n.eq.11) then
           Hermite=1024.d0*(x**10.d0)-23040.d0*(x**8.d0)+161280.d0*(x**6.d0)-403200.d0*(x**4.d0)+302400.d0*(x**2.d0)-30240.d0  
           end if       
           end function
        
        end subroutine gauherm 
      
end module intrastuff
