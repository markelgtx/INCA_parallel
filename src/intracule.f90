!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine intracule(normalize_dm2p)  !Computes vector or radial intracule, or radial integration
                        !need .wfx and .dm2p as input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use geninfo !information about the primitive functions
use intrastuff !subroutines and functions to compute the intracule 
use wfxinfo
use intrainfo !information from the input
use quadratures !nodes and weights of the quadratures, +total grid points
implicit none
logical, intent(in) :: normalize_dm2p
!!!!!!!!!!!!!!!local variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: kk1, kk2 !integers we do not want to read from .dm2
integer :: dfact    !double factorial
integer :: ii,iii,sum, ig, ir, summ, np, sm, smm !some integers
double precision :: R_i_k_2, R_j_l_2  !R_ik=(R_i-R_k)^2
double precision :: Aa_ijkl !grid independent part of the intracule
double precision :: n_prim_t !total normalization factor
double precision :: DMval !value of DM2
double precision :: A_ind !eq.18 (grid independent part)
double precision :: screen1, screen2, lim !integral screenings
double precision :: U_x, U_y, U_z !coef. for 2nd integral screening
double precision, allocatable, dimension(:) :: C_x, C_y, C_z !coefficients of V 
double precision :: V_x, V_y, V_z !eq 16 
double precision, allocatable, dimension(:) :: w_m_array
!Gauss-hermite quadrature
double precision, allocatable, dimension(:) :: rh, w_r !nodes and weights for gauss hermite(coef.)
integer, parameter :: Lmax=21 !maximum angular momentum for primitives
integer, parameter :: nmax = (Lmax+1)/2 + 1
integer, parameter :: npmax = Lmax + 1
integer :: nn !number of gauss hermite nodes (for x y and z)
integer, dimension(3) :: Lrtot  !total angular momentum and number of nodes
double precision, allocatable, dimension(:) :: ipiv !for dgesv in intrastuff
!GRID POINTS
integer :: ngrid !number of grid points
double precision, allocatable, dimension(:,:) :: r !grid points
double precision, dimension(3) :: rp  !prima points(X', Y', Z')
!intracules
double precision :: r_integral                         !radial integral
double precision, allocatable, dimension(:) :: r_intra !radial intracule
double precision, allocatable, dimension(:) :: I_vec   !intracule at a point
!check accuracy of calculations
double precision :: trace_DM2prim, trDM2 !normalized and not normalized DM2prim
integer :: npairs      !number of electron pairs
double precision :: T1, T2, T3, T4, TT1, TT2, TT3, TT4, TT5 !time check
double precision :: Tread, T1screen, T2screen, Tgrid
double precision :: vee
double precision :: intracule_zero
!primitive normalization
double precision, allocatable, dimension(:) :: N_prim
!set counter for primitive quartets
integer :: quartet_count, refval
quartet_count=0

lim=thresh*(dble(nprim)*(dble(nprim)+1.d0)*0.5d0)**(-1.d0) !limit for the 1st integral screening
 
call cpu_time(T1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!Obtain grid points!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (radial_integral) then
    call gridpoints(nradc,nAngc,sfalpha,nquad,cent,Ps) !obtain the grid points
    allocate(r(3,rrgrid))
    r=rrrg
    ngrid=rrgrid
    deallocate(rrrg)
end if

if (radial_plot .or. vee_flag) then
    call gridpoints2(nblock,tart,stp,n_an_per_part)
    allocate(r(3,rgrid))
    r=rpg
    ngrid=rgrid
end if

if (cubeintra) then !vectorial plot
    call gridpoints3(center_i,step_i,np_i)
    !all gridpoints3_corner(center_i,step_i,np_i)    
    allocate(r(3,rgrid))
    r=rg
    ngrid=rgrid
end if

if (intracule_at_zero) then
    ngrid=1
    allocate(r(3,ngrid))
    r(:,1)=0.d0    
end if    
if (intracule_two_points) then
    ngrid=2
    allocate(r(3,ngrid))
    r(1,1)=x_point1
    r(2,1)=y_point1
    r(3,1)=z_point1
    r(1,2)=x_point2
    r(2,2)=y_point2
    r(3,2)=z_point2
end if    
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
allocate(I_vec(ngrid))
I_vec=0.d0
trace_DM2prim=0.d0 
trDM2=0.d0
call cpu_time(T2)
!call intracalc(r,ngrid)
!subroutine intracalc(Computes intracule function with the given points)
open(unit=5,file=dm2name, form='unformatted',access='stream') !open binary file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!Start loop over primitive quartets!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!Following Cioslowski and Liu algorithm!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sum=0     !quartets skipped in the 1st screening
summ=0    !quartets skipped in the 2nd screeening         
rewind(5)    !start reading from the begining of the dm2 file
smm=0
trDM2=0.d0  !sum of all the DM2 terms.
trace_DM2prim=0.d0 !sum of all the normalized DM2 terms.
!Time check
Tread=0.d0
T1screen=0.d0
T2screen=0.d0
Tgrid=0.d0
intracule_zero=0.d0
allocate(rh(nmax)) 
rh=0.d0
allocate(w_r(nmax))
allocate(C_x(npmax)); allocate(C_y(npmax)); allocate(C_z(npmax))
allocate(ipiv(npmax))
! Precompute w_m values for 2nd integral screening
allocate(w_m_array(npmax))
do ii = 1, npmax
    if (ii == 1) then
        w_m_array(ii) = 1.d0
    else
        w_m_array(ii) = (dble(ii) * 0.5d0)**(dble(ii) * 0.5d0) * dexp(-dble(ii) * 0.5d0)
    end if
end do
!perform primitive normalization for DM2 if requested before primitive loop
allocate(N_prim(nprim))
if (normalize_dm2p) then
    do i=1,nprim
        N_prim(i)=(2.d0*Alpha(i)/pi)**(0.75d0)&
        *dsqrt(((4.d0*Alpha(i))**(dble(TMN(i,1)+TMN(i,2)+TMN(i,3))))*&
        dble(dfact(2*TMN(i,1)-1)*dfact(2*TMN(i,2)-1)*dfact(2*TMN(i,3)-1))**(-1.d0))  
    end do
end if
do while (.true.)  !loop for primitive quartets.
    call cpu_time(TT1)    
    read(5,end=100, err=200) kk1,i,j,k,l,DMval,kk2!read a line from binary file .dm2
    !read(5,end=100, err=200) i,j,k,l,DMval
    !write(*,*) i,j,k,l
    if (i.lt.1 .or. i.gt.nprim) goto 200
    if (j.lt.1 .or. k.lt.1 .or. l.lt.1) goto 200
    quartet_count=quartet_count+1
    !if (i.eq.0) goto 100 !file is finished            
    call cpu_time(TT2) 
    Tread=Tread+(TT2-TT1)
    trDM2=trDM2+DMval
    smm=smm+1
    if (normalize_dm2p) then
        !normalization of DM2--> Normalize the primitives
    !    N_prim_i=(2.d0*Alpha(i)/pi)**(0.75d0)&
    !    *dsqrt(((4.d0*Alpha(i))**(dble(TMN(i,1)+TMN(i,2)+TMN(i,3))))*&
    !    dble(dfact(2*TMN(i,1)-1)*dfact(2*TMN(i,2)-1)*dfact(2*TMN(i,3)-1))**(-1.d0))  
   
    !    N_prim_j=(2.d0*Alpha(j)/pi)**(0.75d0)&
    !    *dsqrt(((4.d0*Alpha(j))**(dble(TMN(j,1)+TMN(j,2)+TMN(j,3))))*&
    !    dble(dfact(2*TMN(j,1)-1)*dfact(2*TMN(j,2)-1)*dfact(2*TMN(j,3)-1))**(-1.d0)) 
   
    !    N_prim_k=(2.d0*Alpha(k)/pi)**(0.75d0)&
    !    *dsqrt(((4.d0*Alpha(k))**(dble(TMN(k,1)+TMN(k,2)+TMN(k,3))))*&
    !    dble(dfact(2*TMN(k,1)-1)*dfact(2*TMN(k,2)-1)*dfact(2*TMN(k,3)-1))**(-1.d0)) 
         
    !    N_prim_l=(2.d0*Alpha(l)/pi)**(0.75d0)&
    !    *dsqrt(((4.d0*Alpha(l))**(dble(TMN(l,1)+TMN(l,2)+TMN(l,3))))*&
    !    dble(dfact(2*TMN(l,1)-1)*dfact(2*TMN(l,2)-1)*dfact(2*TMN(l,3)-1))**(-1.d0)) 
   
        !compute DMval with the normalization of primitives
        n_prim_t=N_prim(i)*N_prim(j)*N_prim(k)*N_prim(l)
        DMval=n_prim_t*DMval 
    end if     
    trace_DM2prim=trace_DM2prim+DMval !sum all the DM2 quartets to check accuracy           
    !compute the first variables
    a_ik=Alpha(i)+Alpha(k)  
    a_jl=Alpha(j)+Alpha(l)                  !eqn. 10                         
    e_ik=Alpha(i)*Alpha(k)*a_ik**(-1.d0)
    e_jl=Alpha(j)*Alpha(l)*a_jl**(-1.d0)   
                
    R_i_k_2=(Cartes(Ra(i),1)-Cartes(Ra(k),1))**2.d0+& !this is needed to compute A_ijkl (eqn. 18)
            (Cartes(Ra(i),2)-Cartes(Ra(k),2))**2.d0+&      !R_i_k_2=(R_i-R_k)²
            (Cartes(Ra(i),3)-Cartes(Ra(k),3))**2.d0    
                
    R_j_l_2=(Cartes(Ra(j),1)-Cartes(Ra(l),1))**2.d0+&
            (Cartes(Ra(j),2)-Cartes(Ra(l),2))**2.d0+&
            (Cartes(Ra(j),3)-Cartes(Ra(l),3))**2.d0                   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1st integral screening!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                   
    screen1=dabs(DMval)*dsqrt(J_ik(i,k,R_i_K_2)*J_jl(j,l,R_j_l_2))            
    call cpu_time(TT3)
    T1screen=T1screen-(TT3-TT2)         
    if (screen1.ge.lim) then    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
        !compute the other variables (eqn. 12)
        a_ijkl=a_ik+a_jl
        !and also careful with this!
        e_ijkl=(a_ik*a_jl)*a_ijkl**(-1.d0)
        sqe=dsqrt(e_ijkl)
        !careful with these expressions!
        do ii=1,3                   !X_ik,Y_ik,Z_ik,...
            R_ik(ii)=(Alpha(i)*Cartes(Ra(i),ii)+Alpha(k)*Cartes(Ra(k),ii))*a_ik**(-1.d0)
            R_jl(ii)=(Alpha(j)*Cartes(Ra(j),ii)+Alpha(l)*Cartes(Ra(l),ii))*a_jl**(-1.d0)
            R_ijkl(ii)=(a_ik*R_ik(ii)+a_jl*R_jl(ii))*a_ijkl**(-1.d0)
        end do                      
        Alf_ijkl=0.5d0*(a_ik-a_jl)*a_ijkl**(-1.d0)
        !compute Aa_ijkl (the grid independent part)  !eq.18
        A_ind=(a_ijkl)**(-1.5d0)* dexp(-e_ik*R_i_k_2-e_jl*R_j_l_2)
        Aa_ijkl=DMval*A_ind 
        !!!!!!!!!!!!Calculate coeficients of V!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        sm=0
        !calculate w_m
        ! --- x direction ---
        Lrtot(1)=TMN(i,1)+TMN(j,1)+TMN(k,1)+TMN(l,1) !Lrtot=degree of eqn 15 (I don't use Lmax)    
        call gauherm(Lrtot(1), nn, rh, w_r) !obtain Gauss-Hermite nodes(rh) and weights(w_r) (2L+1)
        np=Lrtot(1)+1              !np is the number of coefficients to represent the polynomial V
        call polycoef(C_x, np, nn, 1, rh, w_r, ipiv) !obtain coefficients of V in x direction
        U_x=dot_product(C_x(1:np),w_m_array(1:np)) !compute U coefficients for 2nd integral screening (eq.43 of the paper)
        ! --- y direction ---
        Lrtot(2)=TMN(i,2)+TMN(j,2)+TMN(k,2)+TMN(l,2) !Lrtot=degree of eqn 15 (I don't use Lmax)
        call gauherm(Lrtot(2), nn, rh, w_r) !obtain Gauss-Hermite nodes(rh) and weights(w_r) (2L+1)
        np=Lrtot(2)+1              !np is the number of coefficients to represent the polynomial V
        call polycoef(C_y, np, nn, 2, rh, w_r, ipiv) !obtain coefficients of V in y direction
        U_y=dot_product(C_y(1:np),w_m_array(1:np)) !compute U coefficients for 2nd integral screening (eq.43 of the paper)  
        ! --- z direction ---
        Lrtot(3)=TMN(i,3)+TMN(j,3)+TMN(k,3)+TMN(l,3) !Lrtot=degree of eqn 15 (I don't use Lmax)
        call gauherm(Lrtot(3), nn, rh, w_r) !obtain Gauss-Hermite nodes(rh) and weights(w_r) (2L+1)
        np=Lrtot(3)+1              !np is the number of coefficients to represent the polynomial V
        call polycoef(C_z, np, nn, 3, rh, w_r, ipiv) !obtain coefficients of V in z direction
        U_z=dot_product(C_z(1:np),w_m_array(1:np)) !compute U coefficients for 2nd integral screening (eq.43 of the paper)             
        !!!!!!!!!!2nd Integral Screening!!!!!!!!!!!!!!!!!!!!!                      
        screen2=dabs(Aa_ijkl*U_x*U_y*U_z) !eqn.40
        call CPU_time(TT4)
        T2screen=T2screen+(TT4-TT3)                
        if (screen2.ge.lim) then 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                             
            !loop over grid points
            do ig=1,nGrid                                             
                V_x=0.d0; V_y=0.d0; V_z=0.d0  
                rp(1)=sqe*(r(1,ig)+r_ik(1)-r_jl(1)) !compute R'(eq. 19)
                rp(2)=sqe*(r(2,ig)+r_ik(2)-r_jl(2))
                rp(3)=sqe*(r(3,ig)+r_ik(3)-r_jl(3)) 
                ! ----Horner's method for x ----
                do iii=1,Lrtot(1)+1   !Ltot+1 is the number of coefficients we have
                    V_x=V_x*rp(1)+C_x(iii)
                end do  !end loop over the polynomial coefficients
                ! ----Horner's method for y ----
                do iii=1,Lrtot(2)+1   !Ltot+1 is the number of coefficients we have
                    V_y=V_y*rp(2)+C_y(iii)
                end do  !end loop over the polynomial coefficients
                ! ----Horner's method for z ----
                do iii=1,Lrtot(3)+1   !Ltot+1 is the number of coefficients we have
                    V_z=V_z*rp(3)+C_z(iii)
                end do  !end loop over the polynomial coefficients
                !!!!!!!!!!!!!Calculate intracule at the grid point!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                I_vec(ig)=I_vec(ig)+Aa_ijkl*dexp(-(rp(1)**2.d0+rp(2)**2.d0+rp(3)**2.d0))*V_x*V_y*V_z
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  
            end do       !End loop over grid points
            call CPU_time(TT5)
            Tgrid=Tgrid+(TT5-TT4)
        else                     
            summ=summ+1  !count quartets that do not pass 2nd screening 
        end if                                                                
    else               
        sum=sum+1    !count quartets that do not pass 1st screening                                 
    end if                
end do !end loop over quartets 
!it comes here when .dm2 file is finished
100 continue
write(*,*) "Reached end of file."
write(*,*) "Total number of quartets read:", quartet_count
goto 300
200 continue
write(*,*) "REad error or corrupted record near quartet", quartet_count+1
write(*,*) "Ended loop for primitives"
write(*,*) "Stopping safely"
goto 300

300 continue
write(*,*) "Loop over primtives completed"
deallocate(C_x)
deallocate(C_y)   !deallocate polynomial coefficients
deallocate(C_z) 
!check if we read a symmetry reduced dmp2prim or the full dm2prim
refval=(((nprim*(nprim+1)/2)+1)*(nprim*(nprim+1)/2))/2
write(*,*) "last quartet", i,j,k,l
write(*,*) "nprim", nprim, "refval", refval, "quartet_count", quartet_count
if (quartet_count.gt.refval) then
    write(*,*) "Full DM2prim read"
    write(*,*) "I(r)=I(-r) intracule symmetry is respected"
else
    write(*,*) "CAUTION!!!!!!!!!!!!!!"
    write(*,*) "Symmetry reduced DM2prim. I(r)=/ I(-r) intracule symmetry is NOT respected"
    if (.not.nosym) then
        write(*,*) "You should run with the option NOSYM to get correct intracule values"
        write(*,*) "Exiting..."
       !stop
    end if
end if
!end subroutine intracalc 
!as output gives a vector with the intracule at the given points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!Compute the integrals using the quadrature !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!weights and the vectorial intracule        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cpu_time(T3)
if (radial_plot) then        !compute radial intracule
    open(unit=3, file=r_plot_name)   
    allocate(r_intra(nradi)) 
    ig=0
    sm=0
    r_intra=0.d0
    ir=0
    do i=1,nradi   
        write(*,*) i, radi(i), smn(i)
        do k=1,smn(i) !sum reduced angular points per radi  
            sm=sm+1            
            r_intra(i)=r_intra(i)+w_ang(sm)*I_vec(sm)!Perform the angular quadrature 
        end do
        write(3,*) radi(i), r_intra(i), r_intra(i)*2.0d0*pi*radi(i)**2, r_intra(i)*2.0d0*pi*radi(i)
        !'(F8.4,1X,D20.16,1X,E15.10,1X,E15.10)'
        !write(3,*) i, radi(i), r_intra(i), r_intra(i)*2.d0*pi*radi(i)**2, r_intra*2*pi*radi(i)
    end do 
    deallocate(r_intra)
    close(3)
end if     
if (radial_integral) then !integral of the intracule
    r_integral=0.d0 
    vee=0.d0
    do i=1,rrgrid
        r_integral=r_integral+rweight(i)*I_vec(i)
        vee=vee+rweight_vee(i)*I_vec(i)
    end do 
    deallocate(rweight)
end if
if (cubeintra) then
    open(unit=3, file=cubeintraname)
    write(3,*) "CUBE FILE"
    write(3,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
    write(3,*) natoms, center_i(1), center_i(2), center_i(3)
    write(3,*) np_i(1), step_i(1), 0.d0, 0.d0
    write(3,*) np_i(2), 0.d0, step_i(2), 0.d0
    write(3,*) np_i(3), 0.d0, 0.d0, step_i(3)
    do i=1,natoms
        write(3,*) an(i), chrg(i), cartes(i,1), cartes(i,2), cartes(i,3)
    end do
    do i=1,ngrid
        !write(3,'(D25.16)') I_vec(i)/2.0d0
        write(3,*) I_vec(i)
    end do  
    close(3)
end if       
if (intracule_at_zero) then
    intracule_zero=I_vec(1)/2.d0
    write(*, '(A, ES25.16)') 'INTRACULE AT ZERO = ', intracule_zero
end if 
if (intracule_two_points) then
    write(*, *) x_point1, y_point1, z_point1, I_vec(1)
    write(*, *) x_point2, y_point2, z_point2, I_vec(2)
end if
! 40 format(6(E16.6E3))
call cpu_time(T4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "-----------Computational time----------------"
write(*,*) "Total CPU time", T4-T1 
write(*,*) "Grid points computing time", T2-T1
write(*,*) "Primitive quartet loop time", T3-T2
write(*,*) "CPU time for intracule integrations", T4-T3
write(*,*) "---Primitive quartet time analysis-----------"
write(*,*) "Reading .dm2p file-->", Tread
write(*,*) "1st primitive screening-->", T1screen
write(*,*) "2nd primitive screening-->", T2screen
write(*,*) "Grid points loop-->", Tgrid
write(*,*) "---------------------------------------------"
write(*,*) "------Grid point + primitive quartet info----"
write(*,*) "Original grid points", maxgrid
write(*,*) "Symmetry reduced grid points", rgrid
write(*,*) "Total number of reduced grid points", rrgrid  
write(*,*) "---------------------------------------------"
write(*,*) "--------------Accuracy check-----------------"
write(*,*) "Sum of all DM2prim terms=", trDM2, trace_DM2prim
write(*,*) "Thresholds used for DM2prim", trsh1, trsh2
write(*,*) "Threshold used for screenings(Tau)=", thresh
write(*,*) "Computed limit for the screenings", lim
write(*,*) "---------------------------------------------"
!!!!!!!!!!!!!!!!!!!!!Print output!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=3,file=outname) !open output file 
write(*,*) "*******Job Finished************"
if (radial_integral) then   
    write(3,*) "--------------RADIAL INTEGRAL----------------" 
    write(3,*) "Number of centres", nquad
    if (nquad.gt.1) then
        write(3,*) "Centres for radial integration"
    end if
    do i=1,nquad
        write(3,*) i, ":::", cent(:,i)
    end do
    write(3,*) "Weights for centres", Ps(:)
    write(3,*) "Alpha parameter", sfalpha(:)
    write(3,*) "Gauss-Legendre nodes", nradc(:)
    write(3,*) "Gauss-Lebedev nodes", nangc(:)
    write(3,*) "---------------------------------------------"
    write(3, '(A, ES25.16)') 'TOTAL VALUE OF INTRACULE = ', r_integral
    write(3, '(A, ES25.16)') 'TOTAL VALUE OF Vee = ', vee
    !write(3, '(A, ES25.16)') 'TOTAL Energy = ', toteng
    npairs=(nelec)*(nelec-1)/2
    write(3,*) "Radial_integral error=", r_integral-dble(npairs)
else if (radial_plot) then
    write(3,*) "RADIAL PLOT, I(S) vs s"    
    write(3,*) "Gauss-Lebedev nodes", n_an_per_part(:)
    write(3,*) "From", radi(1), "to", radi(nradi)
else if (cubeintra) then
    write(3,*) "CUBEFILE GENERATED"
else if (intracule_at_zero) then
    write(3, '(A, ES25.16)') 'INTRACULE AT ZERO = ', intracule_zero
end if
close(3)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function dfact(a) !computes double factorial
integer :: dfact
integer, intent(in) :: a
dfact=1
if (mod(a,2).eq.0) then !n is even
  do i=1,a/2
      dfact=dfact*(2*i)
  end do
else 
  do i=1,(a+1)/2
      dfact=dfact*(2*i-1)
  end do
end if
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
