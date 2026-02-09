!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine becke(rrg,sumq,rgrid,nquad,cent,w_beck,Ps)                                      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! UPDATED: Now implements Salvador's TFVC scheme (JCP 139, 071103 2013)                    !
  ! - Uses robust boundary shift for large size differences                                  !
  ! - Uses k=4 stiffness parameter                                                           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   implicit none
   integer, intent(in) :: rgrid, nquad !total grid, number of quadratures
   double precision, intent(in), dimension(3,nquad) :: cent !center of the quadratures
   double precision, intent(in), dimension(3,rgrid) :: rrg  !(reduced) grid points
   double precision, intent(in), dimension(nquad) :: Ps !weight (POPULATION) of each center
   integer, dimension(nquad) :: sumq             !number of grid points per quad.
   double precision, intent(out), dimension(rgrid) :: w_beck !becke weights for each point
   
   !local variables
   double precision, dimension(rgrid,nquad) :: w_becke !weights of different quad. centers at same point
   integer :: i1,i,j,sm
   double precision, dimension(nquad,nquad) :: Rij, Xi !distance between centers, Xi = Chi ratio
   double precision, dimension(nquad) :: P             !equation 13
   double precision :: ri, rj, s_ij, Ptot
   double precision :: mu_ij, v_ij 
   double precision :: term  ! New local variable for Salvador's formula
  
   !compute R_ij and Xi (Chi)
   do i=1,nquad
    do j=1,nquad
       Rij(i,j)=dsqrt(sum((cent(:,i)-cent(:,j))**2.d0))  !distance between centers
       
       ! TFVC CHANGE 1: Use Square Root of Population ratio
       ! If Ps contains "pair counts" (Volume), sqrt gives "Radius" ratio
       if (Ps(j) .gt. 1.d-12) then
          !Xi(i,j) = dsqrt(Ps(i)/Ps(j))
              Xi(i,j) = (Ps(i)/Ps(j))**(1.d0/3.d0) ! Volume ratio to the 1/3 gives Radius ratio
              !cube root is more robust for large size differences, as it reduces the magnitude of Xi and thus the shift term in Salvador's formula    
       else
          Xi(i,j) = 1.d0 ! Fallback if Pop is zero
       end if
    end do
   end do      
  
   !compute Becke/TFVC weights from all centres for all the grid points
   w_becke=0.d0      
   sm=0      
   do i1=1,rgrid !for all the grid points
      sm=sm+1      !count grid points
      do i=1,nquad   
          P(i)=1.d0
          do j=1,nquad
              if (i.ne.j) then        
                 ri=dsqrt(sum((rrg(:,i1)-cent(:,i))**2.d0)) !distance to center i
                 rj=dsqrt(sum((rrg(:,i1)-cent(:,j))**2.d0)) !distance to center j                
                 ! 1. Calculate Standard Confocal Coordinate mu [-1, 1]
                 if (Rij(i,j) .gt. 1.d-12) then
                     mu_ij=(ri-rj)*(Rij(i,j)**(-1.d0))       !equation 11
                 else
                     mu_ij=0.d0 ! Should not happen if centers are merged
                 end if
  
                 ! 2. TFVC CHANGE 2: Salvador's Transformation (Eq 9 in Source 2)
                 ! This replaces Becke's Eq A2, A5, A6.
                 ! Formula: v = (1 + mu - chi*(1-mu)) / (1 + mu + chi*(1-mu))
                 
                 term = Xi(i,j) * (1.d0 - mu_ij)
                 v_ij = (1.d0 + mu_ij - term) / (1.d0 + mu_ij + term)
                 
                 ! 3. Cutoff Profile
                 s_ij=0.5d0*(1.d0-f_k(v_ij))             !equation 21
                 P(i)=P(i)*s_ij                          !equation 13
              end if                          
          end do
         w_becke(i1,i)=P(i) !store weight of quadrature i at point i1      
      end do 
      Ptot=sum(P)
      if (Ptot .gt. 1.d-15) then
          w_becke(i1,:)=w_becke(i1,:)*(Ptot**(-1.d0)) !normalize weight
      end if
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
        function f_k(val)   !equation 20
        implicit none
        double precision :: f_k
        double precision, intent(in) :: val
        double precision :: vl
        integer :: i
        ! TFVC CHANGE 3: Salvador recommends k=4 for this scheme
        integer, parameter :: k=4 
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