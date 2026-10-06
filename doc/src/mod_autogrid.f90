!! Automatically calculates optimal multicenter quadrature points and weights
!! based on molecular geometry and Bragg-Slater radii.
module autogrid
    use geninfo     ! contains cartes, natoms, an (atomic numbers)
    use intrainfo   ! contains nquad, cent, Ps, sfalpha, etc.
    implicit none
 
 contains
 
    !****************************************************************************
    !! Computes the coordinates, populations, and scaling factors for integration centers.
    subroutine centercalc()
       integer :: i,j,k, nqmax, sm, sm1, smp
       double precision, allocatable, dimension(:,:) :: c1
       double precision, allocatable, dimension(:) :: dist
       double precision, allocatable, dimension(:) :: nelec_c1, nelec_cent
       double precision, allocatable, dimension(:) :: alpha_c1 
       
       double precision, parameter :: trsh=1.d-1, zero=1.d-15, trs=5.d-1
       double precision :: z, k_val, r_j, r_k, pairs_intra, w_rad_sum
       double precision :: total_intra_pairs
       logical :: equal
 
       nqmax=(natoms*(natoms-1))+1
       allocate(c1(3,nqmax))
       allocate(nelec_c1(nqmax))
       allocate(alpha_c1(nqmax)) 
 
       ! --- 1. SETUP ORIGIN CENTER (Index 1) ---
       c1(:,1) = 0.d0 
       
       total_intra_pairs = 0.d0
       w_rad_sum = 0.d0
       
       do i = 1, natoms
           r_j = get_bs_radius(int(an(i)))
           pairs_intra = an(i) * (an(i) - 1.d0) / 2.d0
           total_intra_pairs = total_intra_pairs + pairs_intra
           w_rad_sum = w_rad_sum + (pairs_intra * r_j)
       end do
       
       nelec_c1(1) = total_intra_pairs
       if (total_intra_pairs > zero) then
           alpha_c1(1) = w_rad_sum / total_intra_pairs
       else
           alpha_c1(1) = get_bs_radius(1) 
       end if
 
       ! --- 2. GENERATE BOND CENTERS ---
       sm=1
       write(*,*) "Calculating centers for intracule function..."
       
       do j=1,natoms-1
          do k=j+1,natoms
             if ((.not.nohydro) .or. ((an(j).gt.1).and.(an(k).gt.1))) then
                if (j.ne.k) then
                   ! Positive Vector Center
                   sm=sm+1       
                   c1(:,sm)=cartes(j,:)-cartes(k,:)
                   nelec_c1(sm) = an(j)*an(k)*0.5d0 
                   r_j = get_bs_radius(int(an(j)))
                   r_k = get_bs_radius(int(an(k)))
                   alpha_c1(sm) = (r_j + r_k) / 2.d0
 
                   z=dabs(c1(3,sm))
                   if (z.lt.trs) c1(3,sm)=0.d0 
                   smp=sm  
                   
                   ! Negative Vector Center
                   sm=sm+1 
                   c1(:,sm) = -c1(:,smp) 
                   nelec_c1(sm) = nelec_c1(smp) 
                   alpha_c1(sm) = alpha_c1(smp) 
                end if  
             end if
          end do
       end do
 
       ! --- 3. MERGE EQUIVALENT CENTERS ---
       nqmax=sm
       sm=0 
       
       do i=2,nqmax-1
          z = sum(abs(c1(:,i))) 
          if (z.ne.0.d0) then                     
             do j=i+1,nqmax
                equal=all(abs(c1(:,i)-c1(:,j)).lt.trsh) 
                if (equal) then                        
                   alpha_c1(i) = (alpha_c1(i)*nelec_c1(i) + alpha_c1(j)*nelec_c1(j)) &
                               / (nelec_c1(i) + nelec_c1(j))
                   nelec_c1(i) = nelec_c1(i) + nelec_c1(j)
                   
                   sm=sm+1
                   c1(:,j)=0.d0   
                   nelec_c1(j)=0.d0                   
                end if   
             end do
          end if
       end do
 
       ! --- 4. PACK FINAL ARRAYS ---
       nquad=nqmax-sm
       allocate(cent(3,nquad)); allocate(nelec_cent(nquad))
       allocate(nradc(nquad)); allocate(nangc(nquad))
       allocate(sfalpha(nquad)); allocate(Ps(nquad))
       allocate(dist(nquad))
 
       sm1=0 
       do i=1,nqmax
          z = sum(abs(c1(:,i)))
          if (i == 1 .or. z.ne.0.d0) then  
             sm1=sm1+1
             cent(:,sm1) = c1(:,i)
             Ps(sm1)      = nelec_c1(i)  
             sfalpha(sm1) = alpha_c1(i)  
             write(*,*) "Center ", sm1, " stored. Weight(Pop)=", Ps(sm1), " Alpha=", sfalpha(sm1)
          end if
       end do 
 
       if (betaone) then 
            do i=1,nquad
                 sfalpha(i)=1.d0
                 Ps(i)=1.d0
           end do
       end if
 
       do i=1,nquad 
          nradc(i)=nrad
          nangc(i)=nang
          dist(i)=sqrt(sum(cent(:,i)**2))
       end do
 
       deallocate(c1, nelec_c1, alpha_c1)
    end subroutine centercalc
 
    !****************************************************************************
    !! Returns the Bragg-Slater radius (in Bohr) for a given atomic number.
    double precision function get_bs_radius(iz)
       implicit none
       integer, intent(in) :: iz
       
       select case (iz)
       case (1)  ; get_bs_radius = 0.35d0  
       case (2)  ; get_bs_radius = 0.35d0
       case (3:5); get_bs_radius = 1.45d0 
       case (6)  ; get_bs_radius = 0.70d0
       case (7)  ; get_bs_radius = 0.65d0
       case (8)  ; get_bs_radius = 0.60d0
       case (9)  ; get_bs_radius = 0.50d0
       case (10) ; get_bs_radius = 0.38d0 
       case (15) ; get_bs_radius = 1.00d0
       case (16) ; get_bs_radius = 1.00d0
       case (17) ; get_bs_radius = 1.00d0
       case default; get_bs_radius = 1.00d0
       end select
       
       get_bs_radius = get_bs_radius * 1.8897d0
    end function get_bs_radius
 
 end module autogrid