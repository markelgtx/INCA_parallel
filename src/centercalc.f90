subroutine centercalc()
   use geninfo !contains cartes and natoms, an (atomic numbers)
   use intrainfo !contains nquad and cent
   implicit none

   integer :: i,j,k, nqmax, sm, sm1, smp
   double precision, allocatable, dimension(:,:) :: c1
   double precision, allocatable, dimension(:) :: dist
   double precision, allocatable, dimension(:) :: nelec_c1, nelec_cent
   ! NEW: Array to store alpha candidates before packing into sfalpha
   double precision, allocatable, dimension(:) :: alpha_c1 
   
   double precision, parameter :: trsh=1.d-1, zero=1.d-15, trs=5.d-1
   double precision :: z, k_val, r_j, r_k, pairs_intra, w_rad_sum
   double precision :: total_intra_pairs
   logical :: equal

   ! Function to get radii (defined below or in module)
   double precision :: get_bs_radius

   nqmax=(natoms*(natoms-1))+1
   allocate(c1(3,nqmax))
   allocate(nelec_c1(nqmax))
   allocate(alpha_c1(nqmax)) ! Allocate temp alpha array

   ! --- 1. SETUP ORIGIN CENTER (Index 1) ---
   c1(:,1) = 0.d0 
   
   ! Calculate parameters for the Origin (Super-Center for intra-pairs)
   total_intra_pairs = 0.d0
   w_rad_sum = 0.d0
   
   do i = 1, natoms
       r_j = get_bs_radius(int(an(i)))
       ! Intra-atomic pairs: N*(N-1)/2
       pairs_intra = an(i) * (an(i) - 1.d0) / 2.d0
       total_intra_pairs = total_intra_pairs + pairs_intra
       w_rad_sum = w_rad_sum + (pairs_intra * r_j)
   end do
   
   nelec_c1(1) = total_intra_pairs
   ! Alpha is weighted average radius of all atoms
   if (total_intra_pairs > zero) then
       alpha_c1(1) = w_rad_sum / total_intra_pairs
   else
       alpha_c1(1) = get_bs_radius(1) ! Fallback (e.g. for pure H system? rare)
   end if

   !write(*,*) "Origin Center: Pop=", nelec_c1(1), " Alpha=", alpha_c1(1)

   ! --- 2. GENERATE BOND CENTERS ---
   sm=1
   write(*,*) "Calculating centers for intracule function..."
   
   ! (Your existing loop structure, slightly modified for alpha/pop)
   do j=1,natoms-1
      do k=j+1,natoms
         ! Check conditions (hydro/non-hydro logic preserved)
         if ((.not.nohydro) .or. ((an(j).gt.1).and.(an(k).gt.1))) then
            if (j.ne.k) then
               !write(*,*) "centers between atoms ", j, " and ", k 
               
               ! -- Positive Vector Center --
               sm=sm+1       
               c1(:,sm)=cartes(j,:)-cartes(k,:)
               
               ! POPULATION: N_j * N_k * 0.5 (Half pairs at +s, Half at -s)
               nelec_c1(sm) = an(j)*an(k)*0.5d0 
               
               ! ALPHA: Average of Bragg-Slater radii
               r_j = get_bs_radius(int(an(j)))
               r_k = get_bs_radius(int(an(k)))
               alpha_c1(sm) = (r_j + r_k) / 2.d0

               ! Symmetry check (preserve your logic)
               z=dabs(c1(3,sm))
               if (z.lt.trs) c1(3,sm)=0.d0 

               smp=sm  
               
               ! -- Negative Vector Center --
               sm=sm+1 
               c1(:,sm) = -c1(:,smp) ! Vector inversion
               nelec_c1(sm) = nelec_c1(smp) ! Same population
               alpha_c1(sm) = alpha_c1(smp) ! Same alpha
            end if  
         end if
      end do
   end do

   ! --- 3. MERGE EQUIVALENT CENTERS (Super-Centers) ---
   ! This is crucial: If vectors are identical, we must SUM their populations
   nqmax=sm
   sm=0 ! Reset counter for "zeros" (reusing your logic)
   
   do i=2,nqmax-1
      z = sum(abs(c1(:,i))) ! check if center i was already zeroed
      if (z.ne.0.d0) then                     
         do j=i+1,nqmax
            equal=all(abs(c1(:,i)-c1(:,j)).lt.trsh) ! Use trsh (0.1) or strict tolerance?
            if (equal) then                        
               ! MERGE J INTO I
               !write(*,*) "Merging duplicate center", j, "into", i
               
               ! Update Alpha (Weighted Average)
               alpha_c1(i) = (alpha_c1(i)*nelec_c1(i) + alpha_c1(j)*nelec_c1(j)) &
                           / (nelec_c1(i) + nelec_c1(j))
               
               ! Update Population (Sum)
               nelec_c1(i) = nelec_c1(i) + nelec_c1(j)
               
               ! Zero out center J
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
      ! We keep index 1 (Origin) even if it's 0.0, so handle i=1 explicitly
      if (i == 1 .or. z.ne.0.d0) then  
         sm1=sm1+1
         cent(:,sm1) = c1(:,i)
         ! Store the computed automation parameters
         Ps(sm1)      = nelec_c1(i)  ! Correct Pair Count
         sfalpha(sm1) = alpha_c1(i)  ! Correct Radius Scale
         
         write(*,*) "Center ", sm1, " stored. Weight(Pop)=", Ps(sm1), " Alpha=", sfalpha(sm1)
      end if
   end do 

   if (betaone) then !set alpha and beta parameters to one
        do i=1,nquad
             sfalpha(i)=1.d0
             Ps(i)=1.d0
       end do
    end if

   ! Final parameter settings
   do i=1,nquad 
      nradc(i)=nrad
      nangc(i)=nang
      dist(i)=sqrt(sum(cent(:,i)**2))
   end do

   ! Cleanup
   deallocate(c1, nelec_c1, alpha_c1)

end subroutine centercalc

! --- HELPER FUNCTION ---
double precision function get_bs_radius(iz)
   implicit none
   integer, intent(in) :: iz
   ! Bragg-Slater radii in Angstroms (convert to Bohr if needed!)
   ! Source: Slater 1964
   select case (iz)
   case (1)  ! H
       get_bs_radius = 0.35d0  ! Becke's value (Slater is 0.25)
   case (2)  ! He
       get_bs_radius = 0.35d0
   case (3:5) ! Li, Be, B
       get_bs_radius = 1.45d0 ! Approximate average for simplicity
   case (6)  ! C
       get_bs_radius = 0.70d0
   case (7)  ! N
       get_bs_radius = 0.65d0
   case (8)  ! O
       get_bs_radius = 0.60d0
   case (9)  ! F
       get_bs_radius = 0.50d0
   case (10) ! Ne
       get_bs_radius = 0.38d0 
   case (15) ! P
       get_bs_radius = 1.00d0
   case (16) ! S
       get_bs_radius = 1.00d0
   case (17) ! Cl
       get_bs_radius = 1.00d0
   case default
       get_bs_radius = 1.00d0
   end select
   ! NOTE: If your code uses Bohr, multiply result by 1.8897
   get_bs_radius = get_bs_radius * 1.8897d0
end function get_bs_radius