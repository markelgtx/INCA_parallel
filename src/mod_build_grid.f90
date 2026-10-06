module build_grid
    use quadratures
    implicit none
    !hide internal subroutines from the rest of the program
    private :: becke

contains

!********************************************************************************
    !! Computes grid points for a multicenter radial integral from 0 to infinity.
    !! Radial: Gauss-Legendre. Angular: Gauss-Lebedev. Uses Becke/TFVC weights.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine gridpoints(nradc,nAngc,sfalpha,nquad,cent,Ps) 
        ! Compute grid points to perform a multicenter radial integral from 0 to infinity    
        ! or from a to b (specify in the input file).
        ! Radial quadrature:Gauss-Legendre ; Angular quadrature: Gauss-Levedeb.
        ! Using Becke's weights to compute multicenter integrals.
        ! Neglects grid points with negative z using intracule symmetry, I(r)=I(-r)
        ! Neglects grid points with 0 weight.
        ! The total (reduced) grid ponints are stored in a 3*rrgrid matrix
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        !global variables
        integer, intent(in) :: nquad  !number of quad.
        integer, dimension(nquad), intent(in) :: nradc, nAngc !number of radial points, number of anglular points
        double precision, dimension(nquad),intent(in) :: sfalpha  !scaling factor for the radial points
        double precision, dimension(nquad), intent(in) :: Ps
        double precision, dimension(3,nquad) :: cent !center of the quadratures
        !local variables
        !1-levedeb quadrature
        integer :: nang
        double precision, allocatable, dimension(:) :: x_lb, y_lb, z_lb
        double precision, allocatable, dimension(:) :: Wlb
        !double precision, allocatable, dimension(:) :: theta, phi
        !2-legendre quadrature
        integer :: nrad
        double precision, allocatable, dimension(:) :: xl_i, wl_i, wl2_i
        double precision, allocatable, dimension(:) :: radius !variable change from 0 to infty
        !Legendre+Lebedev combination
        double precision, allocatable, dimension(:) :: weight, weight_vee, srweight, srweight_vee !total weight for each point (before neglecting points)
        !3-becke vi
        double precision, allocatable, dimension(:) :: w_beck
        !-starting grid points 
        double precision, allocatable, dimension(:,:) :: gr1, gr2, gr3 
        integer :: i, ia, ir, sm, smp, sma, smnn, j, k, smr, smpr, ngrid, i1
        integer :: np !number of points
        double precision, parameter :: pi=4.d0*datan(1.d0)
        double precision, parameter :: trsh=1.d-105, trsh2=1.d-106, tol=1.d-100!tol=dsqrt(epsilon(1.d0))
        double precision :: actual_dist
        integer, allocatable, dimension(:) :: temp_map 
    
        !count maximum number of points
        np=0
        do i=1,nquad
            np=np+nradc(i)*nAngc(i)
        end do
        maxgrid=np
        allocate(weight(np))
        !allocate(weight_vee(np))
        allocate(gr1(3,np)) !first grid points
        allocate(gr2(3,np)) !sym reduced grid points
        allocate(smn(nquad)) !counts number of grid points per quadrature
        sm=0 !count total grid points (1-->np)
        smp=0 !count last grid point of previous centre
        smr=0  !count sym reduced grid points
        smpr=0  !count sym reduced gp of previous centre
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i1=1,nquad       !loop over centres
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            nrad=nradc(i1)
            nAng=nAngc(i1)
            !!Obtain radial nodes and weights (Gauss-Legendre)!!!!!!!!!!!!!!!!!!!!!
            allocate(xl_i(nrad))
            allocate(wl_i(nrad))
            call sub_GauLeg(-1.d0,1.d0,xl_i,wl_i,nrad)   
            !!!Obtain angular nodes and weights (Gauss-Lebedev)!!!!!!!!!!!!!!!!!!!!!!!  
            allocate(x_lb(nang)); allocate(y_lb(nang)); allocate(z_lb(nang))
            allocate(Wlb(nang))
            Wlb=0.d0
            if (nAng.eq.6) call LD0006(x_lb,y_lb,z_lb,Wlb,nAng)  !6 grid points
            if (nAng.eq.14) call LD0014(x_lb,y_lb,z_lb,Wlb,nAng)  !14 grid points
            if (nAng.eq.26) call LD0026(x_lb,y_lb,z_lb,Wlb,nAng)  !26 grid points 
            if (nAng.eq.38) call LD0038(x_lb,y_lb,z_lb,Wlb,nAng)  !38 grid points
            if (nAng.eq.50) call LD0050(x_lb,y_lb,z_lb,Wlb,nAng)  !50 grid points
            if (nAng.eq.74) call LD0074(x_lb,y_lb,z_lb,Wlb,nAng)  !74 grid points
            if (nAng.eq.86) call LD0086(x_lb,y_lb,z_lb,Wlb,nAng)  !86 grid points
            if (nAng.eq.110) call LD0110(x_lb,y_lb,z_lb,Wlb,nAng)  !110 grid points
            if (nAng.eq.146) call LD0146(x_lb,y_lb,z_lb,Wlb,nAng)  !146 grid points
            if (nAng.eq.170) call LD0170(x_lb,y_lb,z_lb,Wlb,nAng) !170 grid points
            if (nAng.eq.194) call LD0194(x_lb,y_lb,z_lb,Wlb,nAng) !194 grid points 
            if (nAng.eq.230) call LD0230(x_lb,y_lb,z_lb,Wlb,nAng) !230 grid points
            if (nAng.eq.266) call LD0266(x_lb,y_lb,z_lb,Wlb,nAng) !266 grid points 
            if (nAng.eq.302) call LD0302(x_lb,y_lb,z_lb,Wlb,nAng) !302 grid points
            if (nAng.eq.350) call LD0350(x_lb,y_lb,z_lb,Wlb,nAng) !350 grid points
            if (nAng.eq.434) call LD0434(x_lb,y_lb,z_lb,Wlb,nAng) !434
            if (nAng.eq.590) call LD0590(x_lb,y_lb,z_lb,Wlb,nAng) !590
            if (nAng.eq.770) call LD0770(x_lb,y_lb,z_lb,Wlb,nAng) !770
            if (nAng.eq.974) call LD0974(x_lb,y_lb,z_lb,Wlb,nAng) !974
            if (nAng.eq.1202) call LD1202(x_lb,y_lb,z_lb,Wlb,nAng) !1202
            if (nAng.eq.1454) call LD1454(x_lb,y_lb,z_lb,Wlb,nAng) !1454
            if (nAng.eq.1730) call LD1730(x_lb,y_lb,z_lb,Wlb,nAng) !1730
            if (nAng.eq.2030) call LD2030(x_lb,y_lb,z_lb,Wlb,nAng) !2030
            if (nAng.eq.2354) call LD2354(x_lb,y_lb,z_lb,Wlb,nAng) !2354
            if (nAng.eq.2702) call LD2702(x_lb,y_lb,z_lb,Wlb,nAng) !2702
            if (nAng.eq.3074) call LD3074(x_lb,y_lb,z_lb,Wlb,nAng) !3074
            if (nAng.eq.3470) call LD3470(x_lb,y_lb,z_lb,Wlb,nAng) !3470
            if (nAng.eq.3890) call LD3890(x_lb,y_lb,z_lb,Wlb,nAng) !3890
            if (nAng.eq.4334) call LD4334(x_lb,y_lb,z_lb,Wlb,nAng) !4334
            if (nAng.eq.4802) call LD4802(x_lb,y_lb,z_lb,Wlb,nAng) !4802
            if (nAng.eq.5294) call LD5294(x_lb,y_lb,z_lb,Wlb,nAng) !5294
            if (nAng.eq.5810) call LD5810(x_lb,y_lb,z_lb,Wlb,nAng) !5810
            !!!!!!!!!!!!!Compute grid points of quadrature!!!!!!!!!!!
            nGrid=nAng*nrad
            allocate(radius(nrad))
            !compute radial points for the quadrature from 0 to infty or a to b
            allocate(wl2_i(nrad))  
            if (definite) then
                write(*,*) "Integral from", a, "to", b   
                do ir=1,nrad    
                    radius(ir)=(b-a)*0.5d0*xl_i(ir)+(a+b)*0.5d0
                    wl_i(ir)=(b-a)*0.5d0*wl_i(ir)*radius(ir)**2.d0 !radial integration
                    wl2_i(ir)=(b-a)*0.5d0*wl_i(ir)*radius(ir)   ! for Vee
                end do   
            else
                do ir=1,nrad    
                    radius(ir)=(1.d0+xl_i(ir))/(1.d0-xl_i(ir))*sfalpha(i1)
                    if ((nquad.gt.1).and.(radius(ir).gt.15.d0)) then
                        !write(*,*) "Removing grid point=", radius(ir), "in quadrature center", i1
                        wl2_i(ir)=0.d0
                        wl_i(ir)=0.d0
                    else          
                        wl_i(ir)=2*pi*(2.d0*sfalpha(i1)/((1.d0-xl_i(ir))**2.d0))*wl_i(ir)*radius(ir)*radius(ir) !for I(r)
                    end if
                end do   
            end if  
            !compute grid points for all the becke centers
            smn(i1)=0
            do ir=1,nrad
                do ia=1,nAng
                    sm=sm+1   !count total grid points    
                    gr1(1,sm)=radius(ir)*x_lb(ia)+cent(1,i1)
                    gr1(2,sm)=radius(ir)*y_lb(ia)+cent(2,i1)
                    gr1(3,sm)=radius(ir)*z_lb(ia)+cent(3,i1)
                    !store total grid points (all quadratures) 
                    if (nosym) then
                        smn(i1)=smn(i1)+1 
                    else    
                        if (gr1(3,sm).ge.-tol) then
                            smn(i1)=smn(i1)+1 !sum the number of sym reduced points of each quadrature                  
                        end if
                    end if    
                end do   
            end do  
            sm=smp   !start again from 1st point of center
            if (nosym) then
                do j=1,nrad
                    do k=1,nAng  
                        sm=sm+1
                        weight(sm)=Wlb(k)*Wl_i(j)
                       !weight_vee(sm)=Wlb(k)*wl2_i(j) !for Vee
                    end do             
                end do
            else 
                smr=smpr !start from 1st reduced point of center
                smnn=0   !count number of neglected points by symmetry
                do j=1,nrad
                    do k=1,nAng  
                        sm=sm+1
                        !sum total grid points
                        if (gr1(3,sm).ge.-tol) then        !reduce points by symmetry
                            smr=smr+1                       !sum sym reduced grid points   
                            gr2(:,smr)=gr1(:,sm)           !store sym reduced points
                            if (abs(gr1(3,sm)).lt.tol) then !z is 0, do not multiply by 2
                                weight(smr)=Wlb(k)*Wl_i(j)
                                !weight_vee(smr)=Wlb(k)*wl2_i(j) !for Vee
                            else
                                weight(smr)=2.d0*Wlb(k)*Wl_i(j) !z is positive, use sym (I(z)=I(-z))
                                !weight_vee(smr)=2.d0*Wlb(k)*wl2_i(j) !for Vee
                            end if
                        else 
                            !sym neglected point    
                            smnn=smnn+1
                        end if   
                    end do             
                end do   
                smpr=smr !store last reduced point of the quadrature
            end if    
            smp=sm   !store last point of the quadrature
            deallocate(radius) 
            deallocate(Wlb)
            deallocate(x_lb); deallocate(y_lb); deallocate(z_lb)
            deallocate(xl_i)
            deallocate(wl_i)
            deallocate(wl2_i)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end do !end loop over quadratures
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rgrid=sum(smn) !set total number of grid points
        allocate(gr3(3,rgrid))
        allocate(srweight(rgrid))
        !allocate(srweight_vee(rgrid))
        allocate(temp_map(rgrid))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (nosym) then
            sm=0
            do i=1,nquad                   !transfer points to new matrix 
                do j=1,smn(i)
                    sm=sm+1
                    gr3(:,sm)=gr1(:,sm)
                    srweight(sm)=weight(sm)
                    temp_map(sm)=i !store the center index for each point
                end do  
            end do 
            deallocate(gr1)
        else
            sm=0
            do i=1,nquad                   !store the sym reduced points in a smaller matrix
                do j=1,smn(i) 
                    sm=sm+1
                    gr3(:,sm)=gr2(:,sm)
                    srweight(sm)=weight(sm)
                    temp_map(sm)=i !store the center index for each point
                end do
            end do
            deallocate(gr1)
            deallocate(gr2)    
        end if   
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (nquad.gt.1) then         !we have more than one quadrature centre
            allocate(w_beck(rgrid)) !compute becke weights for each point    
            call becke(gr3,smn,rgrid,nquad,cent,w_beck,Ps)
            do i=1,rgrid
                srweight(i)=srweight(i)*w_beck(i) !store becke weight into the total one
            end do
            deallocate(w_beck)
            !reduce points with 0 weight
            sm=0
            sma=0
            do i=1,rgrid
                if ((srweight(i).le.trsh2)) then !not 1st center and little weight
                    sm=sm+1                     !sum number of neglected points  
                    if (i.gt.smn(1)) then
                        sma=sma+1 !count neglected points in 2nd center 
                    end if        
                end if        
            end do
            rrgrid=rgrid-sm !number of sym and weight reduced points
        else !only 1 quadrature centre--> No Becke    
            rrgrid=rgrid
        end if
    
        !store reduced points in final variables (in mod_quadratures.f90)
        allocate(rweight(rrgrid))     !reduced weight
        allocate(rweight_vee(rrgrid)) !for Vee
        allocate(rrrg(3,rrgrid))      !doubly reduced points
        allocate(point_center_map(rrgrid)) !store the center index for each reduced point
        sm=0
        do i=1,rgrid
            if ((srweight(i).gt.trsh2)) then !remove points with low (zero) weight
                sm=sm+1 
                rweight(sm)=srweight(i) 
                point_center_map(sm)=temp_map(i) !store the center index for each reduced point
                actual_dist=dsqrt(sum(gr3(:,i)**2)) !compute actual distance of the point to the origin
                rweight_vee(sm)=srweight(i)*(actual_dist)**(-1.d0) !for Vee, divide by r to get correct weight for Vee
                rrrg(:,sm)=gr3(:,i)           
            end if
        end do
        !write(*,*) "Grid points:"
        !open(4,file='gridpoints.dat',status='replace')
        !do i=1,rrgrid
        !    write(4,'(I6,3F15.8,F30.8)') i, rrrg(1,i), rrrg(2,i), rrrg(3,i), rweight(i)
        !end do
        deallocate(temp_map)
        deallocate(weight); deallocate(srweight)
        deallocate(gr3)
        close(4)
    end subroutine 

    !********************************************************************************
    !! Computes grid points for a radial plot with manually given radiis.
    subroutine gridpoints2(nblock,tart,step,n_an_per_part)                                                    !
        implicit none
        !global variables
        integer, intent(in) :: nblock !number of parts
        double precision, dimension(nblock), intent(in) :: step
        integer, intent(in), dimension(nblock) :: n_an_per_part !number of steps per part
        double precision, intent(in), dimension(2,nblock) :: tart !range of each Block
        !local variables
        !1-levedeb quadrature
        double precision, allocatable, dimension(:) :: x_lb, y_lb, z_lb 
        double precision, allocatable, dimension(:) :: Wlb
        !double precision, allocatable, dimension(:) :: theta, phi
        double precision, allocatable, dimension(:) :: pweight
        !grid points
        integer, allocatable,dimension(:) :: npb
        integer :: nr_total
        double precision, allocatable, dimension(:,:) :: gpt
        double precision :: z
        integer :: i, ia, ir, sm, smnn, j, n_an, smrad
        double precision, parameter :: pi=4.d0*datan(1.d0)
        double precision, parameter :: trsh=1.d-15, trsh2=1.d-16, tol=dsqrt(epsilon(1.d0))
        !double precision :: a,b
        double precision :: r_start,r_end
        allocate(npb(nblock))
        nr_total=0
        do i=1,nblock
            r_start=tart(1,i)
            r_end=tart(2,i)
            !compute number of points
            npb(i)=int((r_end-r_start)/step(i))+1
            nr_total = nr_total + npb(i)
        end do    
        nradi = nr_total
        allocate(radi(nradi))
        allocate(smn(nradi))
        radi = 0.d0
        smn = 0.d0
        sm = 0.d0
        rpgrid = 0
        do i = 1, nblock
            do j = 0, npb(i)-1
                sm = sm + 1
                radi(sm) = tart(1,i) + step(i)*dble(j)
            end do
            rpgrid = rpgrid + npb(i) * n_an_per_part(i)
        end do
        allocate(gpt(3,rpgrid))
        allocate(pweight(rpgrid))
        smrad=0
        sm=0 !sum for grid points
        smnn=0 !sum symmetry reduced points
        do i=1,nblock !loop for each radius fragment
            n_an=n_an_per_part(i)
            allocate(x_lb(n_an)); allocate(y_lb(n_an)); allocate(z_lb(n_an))
            allocate(Wlb(n_an))
            !!!Obtain angular nodes and weights (Gauss-Lebedev)!!!!!!!!!!!!!!!!!!!!!!!   
            if (n_an.eq.6) call LD0006(x_lb,y_lb,z_lb,Wlb,n_an)  !6 grid points
            if (n_an.eq.14) call LD0014(x_lb,y_lb,z_lb,Wlb,n_an)  !14 grid points
            if (n_an.eq.26) call LD0026(x_lb,y_lb,z_lb,Wlb,n_an)  !26 grid points 
            if (n_an.eq.38) call LD0038(x_lb,y_lb,z_lb,Wlb,n_an)  !38 grid points
            if (n_an.eq.50) call LD0050(x_lb,y_lb,z_lb,Wlb,n_an)  !50 grid points
            if (n_an.eq.74) call LD0074(x_lb,y_lb,z_lb,Wlb,n_an)  !74 grid points
            if (n_an.eq.86) call LD0086(x_lb,y_lb,z_lb,Wlb,n_an)  !86 grid points
            if (n_an.eq.110) call LD0110(x_lb,y_lb,z_lb,Wlb,n_an)  !110 grid points
            if (n_an.eq.146) call LD0146(x_lb,y_lb,z_lb,Wlb,n_an)  !146 grid points
            if (n_an.eq.170) call LD0170(x_lb,y_lb,z_lb,Wlb,n_an) !170 grid points
            if (n_an.eq.194) call LD0194(x_lb,y_lb,z_lb,Wlb,n_an) !194 grid points 
            if (n_an.eq.230) call LD0230(x_lb,y_lb,z_lb,Wlb,n_an) !230 grid points
            if (n_an.eq.266) call LD0266(x_lb,y_lb,z_lb,Wlb,n_an) !266 grid points 
            if (n_an.eq.302) call LD0302(x_lb,y_lb,z_lb,Wlb,n_an) !302 grid points
            if (n_an.eq.350) call LD0350(x_lb,y_lb,z_lb,Wlb,n_an) !350 grid points
            if (n_an.eq.434) call LD0434(x_lb,y_lb,z_lb,Wlb,n_an) !434
            if (n_an.eq.590) call LD0590(x_lb,y_lb,z_lb,Wlb,n_an) !590
            if (n_an.eq.770) call LD0770(x_lb,y_lb,z_lb,Wlb,n_an) !770
            if (n_an.eq.974) call LD0974(x_lb,y_lb,z_lb,Wlb,n_an) !974
            if (n_an.eq.1202) call LD1202(x_lb,y_lb,z_lb,Wlb,n_an) !1202
            if (n_an.eq.1454) call LD1454(x_lb,y_lb,z_lb,Wlb,n_an) !1454
            if (n_an.eq.1730) call LD1730(x_lb,y_lb,z_lb,Wlb,n_an) !1730
            if (n_an.eq.2030) call LD2030(x_lb,y_lb,z_lb,Wlb,n_an) !2030
            if (n_an.eq.2354) call LD2354(x_lb,y_lb,z_lb,Wlb,n_an) !2354
            if (n_an.eq.2702) call LD2702(x_lb,y_lb,z_lb,Wlb,n_an) !2702
            if (n_an.eq.3074) call LD3074(x_lb,y_lb,z_lb,Wlb,n_an) !3074
            if (n_an.eq.3470) call LD3470(x_lb,y_lb,z_lb,Wlb,n_an) !3470
            if (n_an.eq.3890) call LD3890(x_lb,y_lb,z_lb,Wlb,n_an) !3890
            if (n_an.eq.4334) call LD4334(x_lb,y_lb,z_lb,Wlb,n_an) !4334
            if (n_an.eq.4802) call LD4802(x_lb,y_lb,z_lb,Wlb,n_an) !4802
            if (n_an.eq.5294) call LD5294(x_lb,y_lb,z_lb,Wlb,n_an) !5294
            if (n_an.eq.5810) call LD5810(x_lb,y_lb,z_lb,Wlb,n_an) !5810 
            !!!!!!!!Compute grid points!!!!!!!!!!!!!!!!!!!!!!!!!!    
            do ir=1,npb(i) !loop for each radius inside the block
                smrad=smrad+1 !sum over radius
                do ia=1,n_an  !sum over angular points
                    sm=sm+1 !sum grid points
                    if (nosym) then
                        gpt(1,sm)=radi(smrad)*x_lb(ia)!dcos(phi(ia))*dsin(theta(ia))
                        gpt(2,sm)=radi(smrad)*y_lb(ia)!dsin(phi(ia))*dsin(theta(ia))
                        gpt(3,sm)=radi(smrad)*z_lb(ia)!dcos(theta(ia))
                        pweight(sm)=Wlb(ia)
                    else    
                        z=radi(smrad)*z_lb(ia) !compute z
                        if (z.ge.-tol) then
                            smnn=smnn+1 !sum the number of reduced points of each quadrature
                            smn(smrad)=smn(smrad)+1 !sum number of points per radius
                            gpt(1,smnn)=radi(smrad)*x_lb(ia)!dcos(phi(ia))*dsin(theta(ia))
                            gpt(2,smnn)=radi(smrad)*y_lb(ia)!dsin(phi(ia))*dsin(theta(ia))
                            gpt(3,smnn)=z
                            if (abs(z) < tol) then
                                pweight(smnn)=Wlb(ia)
                            else
                                pweight(smnn)=2.d0*Wlb(ia)
                            end if   
                        end if
                    end if    
                end do   
            end do  
            deallocate(Wlb)
            deallocate(x_lb); deallocate(y_lb); deallocate(z_lb)
        end do   !end loop for radius fragment
        allocate(rpg(3,smnn))
        allocate(w_ang(smnn))
        w_ang = 0.0d0
        rpg = 0.0d0 
        do i=1,smnn  !store the data in reduced size matrix
            rpg(:,i)=gpt(:,i) 
            w_ang(i)=pweight(i)          
        end do
        rgrid=smnn
        end subroutine 
        
    !********************************************************************************
    !! Computes uniform 3D grid points for vectorial cubefile generation.    
        subroutine gridpoints3(center_i,step_i,np_i)
            implicit none  
            double precision, intent(in), dimension(3) :: center_i, step_i
            integer, intent(in), dimension(3) :: np_i    
            !local variables
            integer :: i,j,k,sm
            double precision, dimension(3) :: xm  
            rgrid=np_i(1)*np_i(2)*np_i(3)
            allocate(rg(3,rgrid))
            rg=0.d0
            do i=1,3
               if (dmod(dble(np_i(i)),2.d0).eq.0.d0) then !even number
                  xm(i)=center_i(i)-(step_i(i)/2.d0)-((dble(np_i(i))-2.d0)/2.d0)*step_i(i)
               else !odd number  
                  xm(i)=center_i(i)-((dble(np_i(i))-1.d0)/2.d0)*step_i(i)
               end if
            end do
            sm=0
            do i=1,np_i(1)         !Depending on a calculates the point with a different function
                do j=1,np_i(2)
                    do k=1,np_i(3)
                        sm=sm+1
                        rg(1,sm)=xm(1)+step_i(1)*(i-1)
                        rg(2,sm)=xm(2)+step_i(2)*(j-1)
                        rg(3,sm)=xm(3)+step_i(3)*(k-1)
                    end do   
                end do
            end do
        end subroutine gridpoints3    

    !********************************************************************************
    !! Computes uniform 3D grid points for vectorial cubefile generation with corner reference.                  
        subroutine gridpoints3_corner(center_i, step_i, np_i)
            implicit none
            double precision, intent(in), dimension(3) :: center_i, step_i
            integer, intent(in), dimension(3) :: np_i
            integer :: i, j, k, sm
            double precision, dimension(3) :: xm           
            ! Total number of grid points
            rgrid = np_i(1) * np_i(2) * np_i(3)
            allocate(rg(3, rgrid))
            rg = 0.d0           
            ! Starting coordinate (corner reference)
            xm(:) = center_i(:)
           
            ! Fill grid
            sm = 0
            do i = 1, np_i(1)
                do j = 1, np_i(2)
                   do k = 1, np_i(3)
                      sm = sm + 1
                      rg(1, sm) = xm(1) + step_i(1) * (i - 1)
                      rg(2, sm) = xm(2) + step_i(2) * (j - 1)
                      rg(3, sm) = xm(3) + step_i(3) * (k - 1)
                   end do
                end do
            end do
        end subroutine gridpoints3_corner

    !********************************************************************************
    !! Computes Becke weights for each grid point and quadrature center, and stores the final weight for each point in w_beck.
    !! Implements Salvador's TFVC scheme (JCP 139, 071103 2013).        
        subroutine becke(rrg,sumq,rgrid,nquad,cent,w_beck,Ps)                                      
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
                            ! 2. TFVC CHANGE 2: Salvador's Transformation (Eq 9 from Salvador's paper)
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
end module build_grid        