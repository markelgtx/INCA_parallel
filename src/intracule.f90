!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine intracule(normalize_dm2p)  !Computes vector or radial intracule, or radial integration
                        !need .wfx and .dm2p as input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use geninfo !information about the primitive functions
use numbers
use intrastuff !subroutines and functions to compute the intracule 
use wfxinfo
use intrainfo !information from the input
use quadratures !nodes and weights of the quadratures, +total grid points
use omp_lib
implicit none
logical, intent(in) :: normalize_dm2p
!!!!!!!!!!!!!!!local variables from the paper!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer :: kk1, kk2 !integers we do not want to read from .dm2
integer :: ii,iii,sum, ig, ir, summ, np, sm, smm !some integers
! Primitive Data 
double precision :: Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,Xl,Yl,Zl !primitive centers
double precision :: alfi,alfj,alfk,alfl !primitive exponents
integer :: ti,mi,ni,tj,mj,nj,tk,mk,nk,tl,ml,nl !angular momenta of primitives

!Calculation Variables 
double precision :: aik, ajl, eik, ejl            !variables from eqn 10
double precision ::  Xik, Yik, Zik, Xjl, Yjl, Zjl !variables from eqn 10
double precision :: Rik2, Rjl2  !R_ik=(R_i-R_k)^2 !for equation A8
double precision :: Jik, Jjl !J of equation A8 for screening
double precision :: Aijkl !grid independent part of the intracule (eqn 18)
double precision :: screen1, screen2, lim !integral screenings
double precision :: zeta, eijkl, alfijkl, sqe, invz !variables from eqn 12 in module intrastuff
double precision :: Xijkl, Yijkl, Zijkl          !variables from eqn 12

! Screening and Polynomial Coefficients
double precision :: Ux, Uy, Uz !coef. for 2nd integral screening (eqn. 43)
double precision :: Vx, Vy, Vz !eq 16 
integer, dimension(3) :: Lrtot  !total angular momentum and number of nodes
integer :: nn !number of nodes for gauss hermite (eqn 16) output from gauherm, input to polycoef

! Grid and Output
!GRID POINTS
integer :: ngrid !number of grid points
double precision, allocatable, dimension(:,:) :: r !grid points
double precision, dimension(3) :: rp  !prima points(X', Y', Z')
double precision, allocatable, dimension(:) :: rintra !radial intracule
double precision, allocatable, dimension(:) :: Ivec   !intracule at a point
double precision :: rintegral, vee, intraculezero !integral of the intracule and intracule at zero

! Constants and Arrays
integer, parameter :: Lmax=21 !maximum angular momentum for primitives (eq 16)
integer, parameter :: nmax = (Lmax+1)/2 + 1
integer, parameter :: npmax = Lmax + 1
double precision, allocatable, dimension(:) :: wmarray !eq A10, for 2nd integral screening (eq 43)
double precision, allocatable, dimension(:) :: normprim !primitive normalization factors

!primitive counters and total normalization
integer :: i, j, k, l !indices for primitive functions
double precision :: nprimt, DMval !total normalization factor
double precision :: traceDM2prim, trDM2 !normalized and not normalized DM2prim
integer :: npairs      !number of electron pairs
integer :: quartetcount, refval
integer :: quartetafter1stscreening, quartetafter2ndscreening

! Timing variables
double precision :: T1, T2, T3, T4, TT1, TT2, TT3, TT4, TT5 !time check
double precision :: Tread, T1screen, T2screen, Tgrid

!Paralelization variables
integer, parameter :: BATCH_SIZE=40960
integer :: buf_i(BATCH_SIZE),buf_j(BATCH_SIZE),buf_k(BATCH_SIZE),buf_l(BATCH_SIZE)
double precision :: buf_DMval(BATCH_SIZE)
integer :: batch_count, q_idx
logical :: finished_reading

!Thread private pointers/allocatables
double precision, allocatable, dimension(:) :: Ivec_private
double precision, allocatable, dimension(:) :: rh_p, wh_p
double precision, allocatable, dimension(:) :: Cx_p, Cy_p, Cz_p
integer, allocatable, dimension(:) :: ipiv_p
!thead private accumulators
double precision :: Tread_p, T1screen_p, T2screen_p, Tgrid_p
double precision :: traceDM2prim_p, trDM2_p
integer :: q_cnt_p, q_scr1_p, q_scr2_p
! Variables for verified reading
integer :: rec_len, rec_tail, ios
integer :: format_type ! 0=Flat, 1=Standard(24), 2=Bugged(0), 3=Intel64


!start global timer
!call cpu_time(T1)
T1 = omp_get_wtime()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!Obtain grid points!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "Obtaining grid points..."
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
    r(:,1)=ZERO    
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
write(*,*) "Number of grid points to compute intracule:", ngrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!allocate GLOBAL Ivec
allocate(Ivec(ngrid))
Ivec=ZERO
!precompute w_m values
allocate(wmarray(npmax))
do ii = 1, npmax
    if (ii == 1) then
        wmarray(ii) = ONE
    else
        wmarray(ii) = (dble(ii) * HALF)**(dble(ii) * HALF) * dexp(-dble(ii) * HALF)
    end if
end do
!Precompute Primitive Normalization
allocate(normprim(nprim))
if (normalize_dm2p) then
    do ii=1,nprim
        normprim(ii)=(TWO*Alpha(ii)/pi)**(THREEQUARTER)&
        *dsqrt(((FOUR*Alpha(ii))**(dble(TMN(ii,1)+TMN(ii,2)+TMN(ii,3))))*&
        dble(dfact(2*TMN(ii,1)-1)*dfact(2*TMN(ii,2)-1)*dfact(2*TMN(ii,3)-1))**(-ONE))  
    end do
end if
T2=omp_get_wtime()
! Initialize Global Counters
Tread=ZERO; T1screen=ZERO; T2screen=ZERO; Tgrid=ZERO
quartetcount=0; quartetafter1stscreening=0; quartetafter2ndscreening=0
lim=thresh*(dble(nprim)*(dble(nprim)+ONE)*HALF)**(-ONE) !limit for the 1st integral screening
write(*,*) "LIM=", lim
open(unit=5,file=dm2name, form='unformatted',access='stream') ! Existing line
! --- VERIFIED FORMAT DETECTION ---
read(5) rec_len, rec_tail
rewind(5) ! Reset to start

format_type = 0 ! Default to Flat/Stream

if (rec_len == 24 .and. rec_tail == 0) then
    format_type = 3
    write(*,*) "INFO: File Format is 64-bit Fortran (8-byte markers)."
else if (rec_len == 24) then
    format_type = 1
    write(*,*) "INFO: File Format is Standard Fortran (24-byte markers)."
else if (rec_len == 0) then
    format_type = 2
    write(*,*) "INFO: File Format is 'Zero-Marker' Fortran (4-byte markers)."
else
    ! If it's not 24 or 0, it might be flat, OR it might be dangerous garbage.
    ! We assume flat, but we will watch out.
    format_type = 0
    write(*,*) "INFO: File Format appears to be Flat Stream (No markers)."
end if
! ---------------------------------
close(5)
if (format_type.eq.1) then
    open(unit=5,file=dm2name,form='unformatted')
else if (format_type.eq.0) then  
    open(unit=5,file=dm2name, form='unformatted',access='stream')
else 
    write(*,*) "Error: Unsupported file format with markers. Please convert to stream format."
    stop
end if          
finished_reading=.false.
!!!!!!!!!!!!!!!Following Cioslowski and Liu algorithm!!!!!!!!!!!!!!!
write(*,*) "Starting PARALLEL loop over primitive quartets..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2. PARALLEL REGION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2. DEBUGGING PARALLEL REGION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i, j, k, l, DMval, q_idx, ii, iii, ig) &
!$OMP PRIVATE(Xi, Yi, Zi, Xj, Yj, Zj, Xk, Yk, Zk, Xl, Yl, Zl) &
!$OMP PRIVATE(alfi, alfj, alfk, alfl, ti, mi, ni, tj, mj, nj, tk, mk, nk, tl, ml, nl) &
!$OMP PRIVATE(aik, ajl, eik, ejl, Rik2, Rjl2, Jik, Jjl, screen1, screen2) &
!$OMP PRIVATE(zeta, eijkl, sqe, Xik, Yik, Zik, Xjl, Yjl, Zjl, invz, alfijkl) &
!$OMP PRIVATE(Xijkl, Yijkl, Zijkl, Aijkl, Lrtot, nn, np, sm) &
!$OMP PRIVATE(Ux, Uy, Uz, Vx, Vy, Vz, rp) &
!$OMP PRIVATE(rh_p, wh_p, Cx_p, Cy_p, Cz_p, ipiv_p, Ivec_private) &
!$OMP PRIVATE(TT1, TT2, TT3, TT4, TT5, nprimt) &
!$OMP PRIVATE(Tread_p, T1screen_p, T2screen_p, Tgrid_p) &
!$OMP PRIVATE(traceDM2prim_p, trDM2_p, q_cnt_p, q_scr1_p, q_scr2_p)
    
    ! Check which thread is running
    ! (requires 'use omp_lib' at top of subroutine)
!write(*,*) "Thread", omp_get_thread_num(), "started."
!call flush(6) 

!A. Initialize Private Accumulators
Tread_p=ZERO; T1screen_p=ZERO; T2screen_p=ZERO; Tgrid_p=ZERO
q_cnt_p=0; q_scr1_p=0; q_scr2_p=0

!B. Allocate Private Arrays
allocate(Ivec_private(ngrid))
Ivec_private=ZERO
allocate(rh_p(Lmax+1)); rh_p=ZERO
allocate(wh_p(Lmax+1)); wh_p=ZERO
allocate(Cx_p(npmax)); allocate(Cy_p(npmax)); allocate(Cz_p(npmax))
allocate(ipiv_p(npmax))

!C. Loop over batches
do while (.true.)
    
    !---SERIAL: MASTER READS---
    !$OMP SINGLE
        !write(*,*) "Master: Reading Batch..."
        !call flush(6)

        batch_count=0
        do ii=1,BATCH_SIZE
            if (finished_reading) exit
            read(5,end=999,err=200) i, j, k, l, DMval
            
            ! VALIDATION
            if (i < 1 .or. i > nprim) goto 200
            if (j < 1 .or. k < 1 .or. l < 1) goto 200

            batch_count = batch_count + 1
            buf_i(batch_count) = i
            buf_j(batch_count) = j
            buf_k(batch_count) = k
            buf_l(batch_count) = l
            buf_DMval(batch_count) = DMval
            cycle

999             finished_reading = .true.
            !write(*,*) "Master: EOF Reached."
            exit
200         if (batch_count == 0 .and. .not. finished_reading) then
                write(*,*) "Master: End of Valid Data (Padding/Error)."
            end if
            finished_reading=.true.
            exit
        end do
        !write(*,*) "Master: Read", batch_count, "quartets. Finished?", finished_reading
        !call flush(6)
    !$OMP END SINGLE
    
    ! --- BARRIER 1 (Wait for read) ---
    !$OMP BARRIER
    
    if (batch_count == 0) then
         !write(*,*) "Thread", omp_get_thread_num(), "exiting loop (empty batch)."
         !call flush(6)
         exit
    end if

    !---PARALLEL CALCULATION---
    !$OMP DO SCHEDULE(DYNAMIC)
    do q_idx = 1, batch_count
        
        ! Print only every 1000th quartet to avoid spam, but show progress
        !if (mod(q_idx, 1000) == 0) then
        !    write(*,*) "TID", omp_get_thread_num(), "processing q_idx", q_idx
        !    call flush(6)
        !end if

        i = buf_i(q_idx)
        j = buf_j(q_idx)
        k = buf_k(q_idx)
        l = buf_l(q_idx)
        DMval = buf_DMval(q_idx)
        q_cnt_p = q_cnt_p + 1

        if (normalize_dm2p) then
            nprimt=normprim(i)*normprim(j)*normprim(k)*normprim(l)
            DMval=nprimt*DMval
        end if

        ! Load info
        alfi=Alpha(i); alfj=Alpha(j); alfk=Alpha(k); alfl=Alpha(l)   
        ti=TMN(i,1); mi=TMN(i,2); ni=TMN(i,3)
        tj=TMN(j,1); mj=TMN(j,2); nj=TMN(j,3)
        tk=TMN(k,1); mk=TMN(k,2); nk=TMN(k,3)
        tl=TMN(l,1); ml=TMN(l,2); nl=TMN(l,3)
        Xi=Xn(i); Yi=Yn(i); Zi=Zn(i)
        Xj=Xn(j); Yj=Yn(j); Zj=Zn(j)
        Xk=Xn(k); Yk=Yn(k); Zk=Zn(k)    
        Xl=Xn(l); Yl=Yn(l); Zl=Zn(l)

        ! Compute vars
        aik=alfi+alfk; ajl=alfj+alfl
        eik=alfi*alfk*aik**(-ONE); ejl=alfj*alfl*ajl**(-ONE)   
        Rik2=(Xi-Xk)**TWO+(Yi-Yk)**TWO+(Zi-Zk)**TWO    
        Rjl2=(Xj-Xl)**TWO+(Yj-Yl)**TWO+(Zj-Zl)**TWO

        ! 1st Screening
        Jik=JA8(ti,tk,mi,mk,ni,nk,alfi,alfk,xi,xk,yi,yk,zi,zk,Rik2,aik,eik)
        Jjl=JA8(tj,tl,mj,ml,nj,nl,alfj,alfl,xj,xl,yj,yl,zj,zl,Rjl2,ajl,ejl)                                   
        screen1=dabs(DMval)*dsqrt(Jik*Jjl)

        if (screen1.ge.lim) then    
            q_scr1_p=q_scr1_p+1
            zeta=aik+ajl
            eijkl=(aik*ajl)*zeta**(-ONE)
            sqe=dsqrt(eijkl)
            Xik=(alfi*Xi+alfk*Xk)*(aik**(-ONE)); Yik=(alfi*Yi+alfk*Yk)*(aik**(-ONE)); Zik=(alfi*Zi+alfk*Zk)*(aik**(-ONE))
            Xjl=(alfj*Xj+alfl*Xl)*(ajl**(-ONE)); Yjl=(alfj*Yj+alfl*Yl)*(ajl**(-ONE)); Zjl=(alfj*Zj+alfl*Zl)*(ajl**(-ONE))
            Xijkl=(aik*Xik+ajl*Xjl)*zeta**(-ONE); Yijkl=(aik*Yik+ajl*Yjl)*zeta**(-ONE); Zijkl=(aik*Zik+ajl*Zjl)*zeta**(-ONE)
            invz=dsqrt(zeta)**(-ONE); alfijkl=HALF*(aik-ajl)*zeta**(-ONE)
            Aijkl=DMval*(zeta)**(-ONEANDHALF)*dexp(-eik*Rik2-ejl*Rjl2)

            ! Polycoef X
            Lrtot(1)=ti+tj+tk+tl 
            call gauherm(Lrtot(1),nn,rh_p,wh_p) 
            np=Lrtot(1)+1              
            call polycoef(Xi,Xj,Xk,Xl,Xik,Xjl,Xijkl,ti,tj,tk,tl,np,nn,rh_p,wh_p,Cx_p, sqe, invz, alfijkl, ipiv_p) 
            Ux=dot_product(Cx_p(1:np),wmarray(1:np))
            
            ! Polycoef Y
            Lrtot(2)=mi+mj+mk+ml
            call gauherm(Lrtot(2), nn, rh_p, wh_p) 
            np=Lrtot(2)+1              
            call polycoef(Yi,Yj,Yk,Yl,Yik,Yjl,Yijkl,mi,mj,mk,ml,np,nn,rh_p,wh_p,Cy_p, sqe, invz, alfijkl, ipiv_p)
            Uy=dot_product(Cy_p(1:np),wmarray(1:np))
            
            ! Polycoef Z
            Lrtot(3)=ni+nj+nk+nl
            call gauherm(Lrtot(3), nn, rh_p, wh_p) 
            np=Lrtot(3)+1              
            call polycoef(Zi,Zj,Zk,Zl,Zik,Zjl,Zijkl,ni,nj,nk,nl,np,nn,rh_p,wh_p,Cz_p, sqe, invz, alfijkl, ipiv_p)
            Uz=dot_product(Cz_p(1:np),wmarray(1:np))
            
            screen2=dabs(Aijkl*Ux*Uy*Uz)
            if (screen2.ge.lim) then 
                q_scr2_p=q_scr2_p+1
                do ig=1,nGrid                                             
                    Vx=ZERO; Vy=ZERO; Vz=ZERO  
                    rp(1)=sqe*(r(1,ig)+Xik-Xjl); rp(2)=sqe*(r(2,ig)+Yik-Yjl); rp(3)=sqe*(r(3,ig)+Zik-Zjl)
                    do iii=1,Lrtot(1)+1
                        Vx=Vx*rp(1)+Cx_p(iii)
                    end do
                    do iii=1,Lrtot(2)+1
                        Vy=Vy*rp(2)+Cy_p(iii)
                    end do
                    do iii=1,Lrtot(3)+1
                        Vz=Vz*rp(3)+Cz_p(iii)
                    end do
                    Ivec_private(ig)=Ivec_private(ig)+Aijkl*dexp(-(rp(1)**TWO+rp(2)**TWO+rp(3)**TWO))*Vx*Vy*Vz
                end do       
            end if    
        end if   
    end do                                            
    !$OMP END DO 

    ! --- BARRIER 2 (Wait for calc) ---
    !$OMP BARRIER
    
    if (finished_reading .and. batch_count < BATCH_SIZE) then
         !write(*,*) "Thread", omp_get_thread_num(), "finished reading and exiting."
         !call flush(6)
         exit
    end if

end do 

!write(*,*) "Thread", omp_get_thread_num(), "reached CRITICAL."
!call flush(6)

!$OMP CRITICAL
    Ivec(:)=Ivec(:)+Ivec_private(:)
    quartetcount=quartetcount+q_cnt_p
    quartetafter1stscreening=quartetafter1stscreening+q_scr1_p
    quartetafter2ndscreening=quartetafter2ndscreening+q_scr2_p

    ! ... add other reductions here ...
!$OMP END CRITICAL

!write(*,*) "Thread", omp_get_thread_num(), "reached END PARALLEL."
!call flush(6)

deallocate(Ivec_private); deallocate(rh_p); deallocate(wh_p)
deallocate(Cx_p); deallocate(Cy_p); deallocate(Cz_p); deallocate(ipiv_p)
!$OMP END PARALLEL
rewind(5)
close(5)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3. POST-PROCESSING (SERIAL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "Loop over primtives completed"
!check if we read a symmetry reduced dmp2prim or the full dm2prim
refval=(((nprim*(nprim+1)/2)+1)*(nprim*(nprim+1)/2))/2
if (quartetcount.gt.refval) then
    write(*,*) "Full DM2prim read"
else
    write(*,*) "Symmetry reduced DM2prim."
    if (.not.nosym) then
        write(*,*) "Run with NOSYM option."
    end if
end if

!call cpu_time(T3)
T3=omp_get_wtime()

if (radial_plot) then        !compute radial intracule
    open(unit=3, file=r_plot_name)   
    allocate(rintra(nradi)) 
    ig=0; sm=0; rintra=ZERO; ir=0
    do i=1,nradi   
        do k=1,smn(i) !sum reduced angular points per radi  
            sm=sm+1            
            rintra(i)=rintra(i)+w_ang(sm)*Ivec(sm)!Perform the angular quadrature 
        end do
        write(3,*) radi(i), rintra(i), rintra(i)*TWO*pi*radi(i)**TWO, rintra(i)*TWO*pi*radi(i)
        !'(F8.4,1X,D20.16,1X,E15.10,1X,E15.10)'
        !write(3,*) i, radi(i), rintra(i), rintra(i)*TWO*pi*radi(i)**2, rintra*2*pi*radi(i)
    end do 
    deallocate(rintra)
    close(3)
end if    

if (radial_integral) then !integral of the intracule
    rintegral=ZERO 
    vee=ZERO
    do i=1,rrgrid
        rintegral=rintegral+rweight(i)*Ivec(i)
        vee=vee+rweight_vee(i)*Ivec(i)
        !write(*,*) i, r(1,i),r(2,i),r(3,i), Ivec(i)
    end do 
    deallocate(rweight)
end if

if (cubeintra) then
    open(unit=3, file=cubeintraname)
    write(3,*) "CUBE FILE"
    write(3,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
    write(3,*) natoms, center_i(1), center_i(2), center_i(3)
    write(3,*) np_i(1), step_i(1), ZERO, ZERO
    write(3,*) np_i(2), ZERO, step_i(2), ZERO
    write(3,*) np_i(3), ZERO, ZERO, step_i(3)
    do i=1,natoms
        write(3,*) an(i), chrg(i), cartes(i,1), cartes(i,2), cartes(i,3)
    end do
    do i=1,ngrid
        !write(3,'(D25.16)') Ivec(i)/2.0d0
        write(3,*) Ivec(i)
    end do  
    close(3)
end if       
if (intracule_at_zero) then
    intraculezero=Ivec(1)/TWO
    write(*, '(A, ES25.16)') 'INTRACULE AT ZERO = ', intraculezero
end if 
if (intracule_two_points) then
    write(*, *) x_point1, y_point1, z_point1, Ivec(1)
    write(*, *) x_point2, y_point2, z_point2, Ivec(2)
end if
! 40 format(6(E16.6E3))
!call cpu_time(T4)
T4=omp_get_wtime()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "-----------Computational time----------------"
write(*,*) "Total CPU time", T4-T1 
write(*,*) "Grid points and preparameter computing time", T2-T1
write(*,*) "Primitive quartet loop time", T3-T2
write(*,*) "CPU time for intracule integrations", T4-T3
write(*,*) "---------------------------------------------"
write(*,*) "------Grid point + primitive quartet info----"
write(*,*) "Original grid points", maxgrid
write(*,*) "Symmetry reduced grid points", rgrid
write(*,*) "Total number of reduced grid points", rrgrid  
write(*,*) "---------------------------------------------"
write(*,*) "--------------Accuracy check-----------------"
write(*,*) "Quartets in the dm2p file", quartetcount
write(*,*) "Quartets after 1st screening", quartetafter1stscreening
write(*,*) "Quartets after 2nd screening (computed quartets)", quartetafter2ndscreening
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
    write(3, '(A, ES25.16)') 'TOTAL VALUE OF INTRACULE = ', rintegral
    write(3, '(A, ES25.16)') 'TOTAL VALUE OF Vee = ', vee
    !write(3, '(A, ES25.16)') 'TOTAL Energy = ', toteng
    npairs=(nelec)*(nelec-1)/2
    write(3,*) "Radial_integral error=", rintegral-dble(npairs)
    write(3,*) nelec, npairs
else if (radial_plot) then
    write(3,*) "RADIAL PLOT, I(S) vs s"    
    write(3,*) "Gauss-Lebedev nodes", n_an_per_part(:)
    write(3,*) "From", radi(1), "to", radi(nradi)
else if (cubeintra) then
    write(3,*) "CUBEFILE GENERATED"
else if (intracule_at_zero) then
    write(3, '(A, ES25.16)') 'INTRACULE AT ZERO = ', intraculezero
end if
close(3)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
