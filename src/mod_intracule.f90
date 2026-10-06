!! Handles the core execution loop for evaluating the electron intracule.
!! Includes OpenMP parallelization and Cioslowski-Liu integral screening.
module intracule
    use geninfo 
    use numbers
    use intrastuff 
    use wfxinfo
    use intrainfo 
    use build_grid  ! Updated from build_grid
    use quadratures 
    use omp_lib
    implicit none

contains

!********************************************************************************
    !! Computes vector or radial intracule, or radial integration over grid points.
    subroutine calculate_intracule(normalize_dm2p)  
        logical, intent(in) :: normalize_dm2p
        integer :: ii,iii, ig, ir, summ, np, sm, smm !some integers
        integer :: c_idx !index for point_center_map to get the center index for a grid point
        double precision :: u_max
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
        !output
        integer, parameter :: iout=15
        open(unit=iout,file=outname)

        ! start global timer
        T1 = omp_get_wtime()
        
        ! Obtain grid points
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

        !Precompute Gauss-Hermite quadrature for the maximum possible angular momentum (Lmax) 
        call init_gauherm()

        !ADD HERE OTHER PRECOMPUTATIONS (e.g., screening thresholds, etc.)
        T2=omp_get_wtime()

        ! Initialize Global Counters
        Tread=ZERO; T1screen=ZERO; T2screen=ZERO; Tgrid=ZERO
        quartetcount=0; quartetafter1stscreening=0; quartetafter2ndscreening=0
        lim=thresh*(dble(nprim)*(dble(nprim)+ONE)*HALF)**(-ONE) !limit for the 1st integral screening
        open(unit=5,file=dm2name, form='unformatted',access='stream') ! Existing line
        
        ! --- DM2P FORMAT DETECTION ---
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
        close(5) 
        ! --------------------------------
        ! Reopen the file with the correct format
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
        !                        PARALLEL REGION
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
        !call flush(6) 

        ! Initialize Private Accumulators
        Tread_p=ZERO; T1screen_p=ZERO; T2screen_p=ZERO; Tgrid_p=ZERO
        q_cnt_p=0; q_scr1_p=0; q_scr2_p=0

        ! Allocate Private Arrays
        allocate(Ivec_private(ngrid))
        Ivec_private=ZERO
        allocate(rh_p(Lmax+1)); rh_p=ZERO
        allocate(wh_p(Lmax+1)); wh_p=ZERO
        allocate(Cx_p(npmax)); allocate(Cy_p(npmax)); allocate(Cz_p(npmax))
        allocate(ipiv_p(npmax))

        do while (.true.)
            !---SERIAL: MASTER READS---
            !$OMP SINGLE
            !call flush(6)
            batch_count=0
            !serially read a batch of data for the next parallel section. 
            do ii=1,BATCH_SIZE
                if (finished_reading) exit
                read(5,end=999,err=200) i, j, k, l, DMval                            
                ! end reading the .dm2p file if i,j,k,l are out of bounds or if we hit EOF. 
                if (i < 1 .or. i > nprim) goto 200
                if (j < 1 .or. k < 1 .or. l < 1) goto 200
                
                !batch the data for parallel processing
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

200             if (batch_count == 0 .and. .not. finished_reading) then
                    write(*,*) "Master: End of Valid Data (Error)."
                end if
                finished_reading=.true.
                exit
            end do
            !call flush(6)
            !$OMP END SINGLE
    
            ! --- BARRIER 1 (Wait for read) ---
            !$OMP BARRIER
    
            if (batch_count == 0) then
                write(*,*) "Thread", omp_get_thread_num(), "exiting loop (empty batch)."
                !call flush(6)
                exit
            end if

            !---PARALLEL CALCULATION---
            !$OMP DO SCHEDULE(DYNAMIC)
            do q_idx = 1, batch_count        
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
                ! Load primitive info
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
                    !call gauherm(Lrtot(1),nn,rh_p,wh_p)
                    nn = pre_n(Lrtot(1))
                    rh_p(1:nn) = pre_rh(1:nn, Lrtot(1))
                    wh_p(1:nn) = pre_wh(1:nn, Lrtot(1))
                    np=Lrtot(1)+1          
                    call polycoef(Xi,Xj,Xk,Xl,Xik,Xjl,Xijkl,ti,tj,tk,tl,np,nn,rh_p,wh_p,Cx_p, sqe, invz, alfijkl, ipiv_p) 
                    Ux=dot_product(Cx_p(1:np),wmarray(1:np))            
                    ! Polycoef Y
                    Lrtot(2)=mi+mj+mk+ml
                    !call gauherm(Lrtot(2), nn, rh_p, wh_p) 
                    nn = pre_n(Lrtot(2))
                    rh_p(1:nn) = pre_rh(1:nn, Lrtot(2))
                    wh_p(1:nn) = pre_wh(1:nn, Lrtot(2))
                    np=Lrtot(2)+1              
                    call polycoef(Yi,Yj,Yk,Yl,Yik,Yjl,Yijkl,mi,mj,mk,ml,np,nn,rh_p,wh_p,Cy_p, sqe, invz, alfijkl, ipiv_p)
                    Uy=dot_product(Cy_p(1:np),wmarray(1:np))            
                    ! Polycoef Z
                    Lrtot(3)=ni+nj+nk+nl
                    !call gauherm(Lrtot(3), nn, rh_p, wh_p) 
                    nn = pre_n(Lrtot(3))
                    rh_p(1:nn) = pre_rh(1:nn, Lrtot(3))
                    wh_p(1:nn) = pre_wh(1:nn, Lrtot(3))                    
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
                write(*,*) "Thread", omp_get_thread_num(), "finished reading and exiting."
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
        !POST-PROCESSING: write out results, perform radial integration if needed, etc.
        write(*,*) ">> Intracule integrations completed. Writing output..."

        write(iout,*) "=========================================================="
        write(iout,*) "               INCA INTRACULE CALCULATION                 "
        write(iout,*) "=========================================================="
        
        ! 1. Accuracy and Grid Information
        write(iout,*) "--- GRID AND ACCURACY SETTINGS ---"
        write(iout,*) "Original grid points:             ", maxgrid
        write(iout,*) "Symmetry reduced grid points:     ", rgrid
        write(iout,*) "Total surviving grid points:      ", rrgrid  
        write(iout,*) "Quartets in the dm2p file:        ", quartetcount
        write(iout,*) "Quartets after 1st screening:     ", quartetafter1stscreening
        write(iout,*) "Quartets actually computed:       ", quartetafter2ndscreening
        write(iout,'(A, ES12.4)') "Screening limit (lim):            ", lim
        write(iout,*) "Thresholds used for DM2prim:      ", trsh1, trsh2
        
        ! 2. Calculation Specific Results
        write(iout,*) ""
        write(iout,*) "--- RESULTS ---"
        
        if (radial_integral) then
            write(iout,*) "Calculation Type: Radial Integral" 
            write(iout,*) "Number of centres: ", nquad
            if (nquad.gt.1) write(iout,*) "Centres for radial integration:"
            do i=1,nquad
                write(iout,*) i, ":::", cent(:,i)
            end do
            write(iout,*) "Weights for centres: ", Ps(:)
            write(iout,*) "Alpha parameter:     ", sfalpha(:)
            write(iout,*) "Gauss-Legendre nodes:", nradc(:)
            write(iout,*) "Gauss-Lebedev nodes: ", nangc(:)
        
            rintegral=ZERO 
            vee=ZERO
            allocate(rint_local(nquad))
            allocate(vee_local(nquad))
            rint_local(:)=ZERO
            vee_local(:)=ZERO    
            
            ! Compute the integral values
            do i=1,rrgrid
                c_idx=point_center_map(i)
                rint_local(c_idx)=rint_local(c_idx)+rweight(i)*Ivec(i)
                vee_local(c_idx)=vee_local(c_idx)+rweight_vee(i)*Ivec(i)
                rintegral=rintegral+rweight(i)*Ivec(i)
                vee=vee+rweight_vee(i)*Ivec(i)
            end do 
        
            write(iout,*) "---------------------------------------------"
            write(iout, '(A, ES25.16)') 'TOTAL VALUE OF INTRACULE = ', rintegral
            write(iout, '(A, ES25.16)') 'TOTAL VALUE OF Vee       = ', vee
            npairs=(nelec)*(nelec-1)/2
            write(iout, '(A, ES15.6)')  'Radial_integral error    = ', rintegral-dble(npairs)
            write(iout,*) "---------------------------------------------"
            
            write(iout,*) "-------LOCAL DECOMPOSITION (TFVC)-------"
            write(iout,'(A6,1X,A20,1X,A15,1X,A15)') "Center", "Coordinates", "Pairs(N_pairs)", "Vee(Hartree)"
            do i=1,nquad
                if (i == 1) then
                    write(iout,'(I6,1X,3F7.3,1X,F15.8,1X,F15.8,1X,A)') &
                    i, cent(:,i), rint_local(i), vee_local(i), "(INTRA-ATOMIC)"
                else
                    write(iout,'(I6,1X,3F7.3,1X,F15.8,1X,F15.8,1X,A)') &
                    i, cent(:,i), rint_local(i), vee_local(i), "(INTER/DISP)"
                end if
            end do
            write(iout,*) "---------------------------------------------"
        
        else if (radial_plot) then
            write(iout,*) "Calculation Type: Radial Plot, I(s) vs s"    
            write(iout,*) "Radius evaluated from ", radi(1), " to ", radi(nradi), " bohr"
            write(iout,*) "Data written to: ", trim(r_plot_name)
        
        else if (cubeintra) then
            write(iout,*) "Calculation Type: Vectorial Cubefile"
            write(iout,*) "Data written to: ", trim(cubeintraname)
        end if 
        
        if (intracule_at_zero) then
            intraculezero=Ivec(1)/TWO
            write(iout, '(A, ES25.16)') 'INTRACULE AT ZERO = ', intraculezero
        end if
        
        if (intracule_two_points) then
            write(iout, *) "Point 1: ", x_point1, y_point1, z_point1, " Intracule: ", Ivec(1)
            write(iout, *) "Point 2: ", x_point2, y_point2, z_point2, " Intracule: ", Ivec(2)
        end if
        
        ! 3. Timings
        T4=omp_get_wtime()
        write(iout,*) ""
        write(iout,*) "--- COMPUTATIONAL TIME ---"
        write(iout,'(A, F10.3, A)') "Total CPU time:             ", T4-T1, " seconds"
        write(iout,'(A, F10.3, A)') "Grid generation time:       ", T2-T1, " seconds"
        write(iout,'(A, F10.3, A)') "Integration loop time:      ", T4-T3, " seconds"
        write(iout,*) "=========================================================="
        close(iout)
        
        
        ! =================================================================
        ! B. WRITE SECONDARY DATA FILES
        ! =================================================================
        
        if (radial_plot) then        
            ! The .rad file (Unit 16)
            open(unit=16, file=r_plot_name, status='replace')   
            write(16,*) "# --------------------------------------------------------"
            write(16,*) "# RADIAL INTRACULE DATA"
            write(16,*) "# Col 1: Distance s (bohr)"
            write(16,*) "# Col 2: I(s) (Radial Intracule)"
            write(16,*) "# Col 3: P(s) = 2 * pi * s^2 * I(s) (Probability Density)"
            write(16,*) "# Col 4: s * I(s)"
            write(16,*) "# --------------------------------------------------------"
            allocate(rintra(nradi)) 
            ig=0; sm=0; rintra=ZERO; ir=0
            do i=1,nradi   
                do k=1,smn(i) 
                    sm=sm+1            
                    rintra(i)=rintra(i)+w_ang(sm)*Ivec(sm)
                end do
                write(16,'(F12.6, 3ES20.8)') radi(i), rintra(i), rintra(i)*TWO*pi*radi(i)**TWO, rintra(i)*TWO*pi*radi(i)
            end do 
            deallocate(rintra)
            close(16)
        end if
        
        if (radial_integral) then
        ! --- Smart Filename Generation ---
            block
                character(len=200) :: smear_file
                integer :: dot_idx
                dot_idx = index(trim(outname), '.', back=.true.)
                if (dot_idx > 0) then
                    smear_file = outname(1:dot_idx-1) // "_smearing.dat"
                else
                    smear_file = trim(outname) // "_smearing.dat"
                end if
                !set u_max as the longest distance integration center + 5 bohr buffer
                u_max=maxval(dsqrt(sum(cent(:,1:nquad)**TWO, dim=1))) + 5.d0   
                call radial_smearing(r, rweight, Ivec, rrgrid, u_max, 1000, 0.15d0, smear_file)
            end block   
        end if
        
        if (cubeintra) then
            ! The .cube file (Unit 17)
            open(unit=17, file=cubeintraname, status='replace')
            write(17,*) "CUBE FILE"
            write(17,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
            write(17,*) natoms, center_i(1), center_i(2), center_i(3)
            write(17,*) np_i(1), step_i(1), ZERO, ZERO
            write(17,*) np_i(2), ZERO, step_i(2), ZERO
            write(17,*) np_i(3), ZERO, ZERO, step_i(3)
            do i=1,natoms
                write(17,*) an(i), chrg(i), cartes(i,1), cartes(i,2), cartes(i,3)
            end do
            do i=1,ngrid
                write(17,*) Ivec(i)
            end do  
            close(17)
        end if 
        
        write(*,*) ">> INCA job finished successfully. Output written to: ", trim(outname)
        
    end subroutine calculate_intracule

    !********************************************************************************
    !! Reconstructs Radial Intracule P(u) from 3D Multicenter Grid using Gaussian Smearing
    subroutine radial_smearing(coords, weights, density, n_grid, u_max, n_bins,sigma,filename)
        implicit none
        ! Arguments
        integer, intent(in) :: n_grid, n_bins
        double precision, intent(in) :: u_max, sigma
        double precision, dimension(3, n_grid), intent(in) :: coords   ! rrrg (3, rrgrid)
        double precision, dimension(n_grid), intent(in) :: weights    ! rweight (rrgrid)
        double precision, dimension(n_grid), intent(in) :: density    ! Ivec (rrgrid)    
        double precision, allocatable :: u_axis(:), P_u(:), I_avg(:)
        double precision :: u_i, val_i, dist, gauss, prefactor, r2_inv
        double precision :: bin_width
        integer :: i, j, j_min, j_max
        character(len=*), intent(in) :: filename
        allocate(u_axis(n_bins))
        allocate(P_u(n_bins))     ! The Probability Density P(u)
        allocate(I_avg(n_bins))   ! The Spherically Averaged Density I(u)
        P_u = 0.d0
        write(*,*) "Allocated arrays for smearing. Starting setup...", n_bins
        ! 1. Setup Radial Axis
        bin_width = u_max / dble(n_bins - 1)
        do j = 1, n_bins
            u_axis(j) = dble(j-1) * bin_width
        end do

        ! Precompute Gaussian constant
        prefactor = 1.d0 / (sigma * sqrt(2.d0 * pi))
        write(*,*) "We approximate the radial intracule using Gaussian smearing. The obtained I(s) vs s"
        write(*,*) "plot is valid for a qualitative analysis of the radial intracule peaks but not valid"
        write(*,*) "for a qualitative calculation of moments of the visualization of dispersion holes"
        write(*,*) "Starting Radial Smearing..."
        write(*,*) "Grid Points:", n_grid, " Sigma:", sigma
 
        ! 2. Smearing Loop (Project 3D volume to 1D axis)
        ! We loop over every grid point in the multicenter grid
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, u_i, val_i, j, dist, gauss, j_min, j_max) REDUCTION(+:P_u)
        do i = 1, n_grid
            ! Calculate distance of this grid point from origin
            u_i = sqrt(coords(1,i)**2 + coords(2,i)**2 + coords(3,i)**2)
            ! The "Mass" this point contributes = Volume_Weight * Density
            val_i = weights(i) * density(i)
        
            ! Optimization: Only update bins within +/- 4 sigma
            j_min = max(1, int((u_i - 4.d0*sigma)/bin_width) + 1)
            j_max = min(n_bins, int((u_i + 4.d0*sigma)/bin_width) + 1)
        
            do j = j_min, j_max
                dist = u_axis(j) - u_i
                ! Gaussian Kernel
                gauss = prefactor * exp(-(dist**2) / (2.d0 * sigma**2))
                ! Add to probability density P(u)
                ! We multiply by bin_width later to normalize the integral, 
                ! but for the density value at a point, we add the kernel value.
                P_u(j) = P_u(j) + val_i * gauss
            end do
        end do
        !$OMP END PARALLEL DO

        ! 3. Post-Processing & Output
        open(unit=25, file=filename, status='replace')
        write(25,*) "#   u (bohr)       P(u) [Prob]      I_avg(u) [Density]"
    
        do j = 1, n_bins
            if (u_axis(j) > 1.d-6) then
                r2_inv = 1.d0 / (4.d0 * pi * u_axis(j)**2)
                I_avg(j) = P_u(j) * r2_inv
            else
                I_avg(j) = 0.d0 ! Avoid singularity at u=0
            end if        
            write(25, '(3ES16.8)') u_axis(j), P_u(j), I_avg(j)
        end do
        close(25)
        write(*,*) "Radial smearing completed. Output in radial_intracule_smearing.dat"
        deallocate(u_axis, P_u, I_avg)
    end subroutine radial_smearing

end module intracule    
