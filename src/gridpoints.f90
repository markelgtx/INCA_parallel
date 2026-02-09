!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gridpoints(nradc,nAngc,sfalpha,nquad,cent,Ps) 
! Compute grid points to perform a multicenter radial integral from 0 to infinity    
!or from a to b (specify in the input file).
! Radial quadrature:Gauss-Legendre ; Angular quadrature: Gauss-Levedeb.
! Using Becke's weights to compute multicenter integrals.
! Neglects grid points with negative z using intracule symmetry, I(r)=I(-r)
! Neglects grid points with 0 weight.
! The total (reduced) grid ponints are stored in a 3*rrgrid matrix
use quadratures !contains output variables rg, rrg, and rgrid that will be used in intracule.f90       !
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
double precision, parameter :: trsh=1.d-15, trsh2=1.d-16, tol=dsqrt(epsilon(1.d0))
double precision :: actual_dist
!count maximum number of points
np=0
do i=1,nquad
  np=np+nradc(i)*nAngc(i)
end do
maxgrid=np
write(*,*) "We are in subroutine gridpoints"
write(*,*) "Error?"
allocate(weight(np))
allocate(weight_vee(np))
allocate(gr1(3,np)) !first grid points
allocate(gr2(3,np)) !sym reduced grid points
allocate(smn(nquad)) !counts number of grid points per quadrature
write(*,*) "Error2="
sm=0 !count total grid points (1-->np)
smp=0 !count last grid point of previous centre
smr=0  !count sym reduced grid points
smpr=0  !count sym reduced gp of previous centre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i1=1,nquad       !loop over centres
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nrad=nradc(i1)
    nAng=nAngc(i1)
    write(*,*) "nrad=,nang=", nrad, nang
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
                write(*,*) "Removing grid point=", radius(ir), "in quadrature center", i1
                wl2_i(ir)=0.d0
                wl_i(ir)=0.d0
            else          
                !wl2_i(ir)=2*pi*(2.d0*sfalpha(i1)/((1.d0-xl_i(ir))**2.d0))*wl_i(ir)*radius(ir)   ! for Vee
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
    !write(*,*) "Total number of gp after center",i1,"=", sm, ngrid
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
    !    write(*,*) "*********Neglecting points by symmetry*******************"
        do j=1,nrad
            do k=1,nAng  
                sm=sm+1
                !sum total grid points
                if (gr1(3,sm).ge.-tol) then        !reduce points by symmetry
                    smr=smr+1                       !sum sym reduced grid points   
                    gr2(:,smr)=gr1(:,sm)           !store sym reduced points
                    if (abs(gr1(3,sm)).lt.tol) then !z is 0, do not multiply by 2
                        weight(smr)=Wlb(k)*Wl_i(j)
                        weight_vee(smr)=Wlb(k)*wl2_i(j) !for Vee
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
       ! write(*,*) smnn, "points have been neglected in center number", i1
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
!write(*,*) "There are", rgrid, "grid points"
allocate(gr3(3,rgrid))
allocate(srweight(rgrid))
allocate(srweight_vee(rgrid))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (nosym) then
    sm=0
    do i=1,nquad                   !transfer points to new matrix 
        do j=1,smn(i)
            sm=sm+1
            gr3(:,sm)=gr1(:,sm)
            srweight(sm)=weight(sm)
            !srweight_vee(sm)=weight_vee(sm) !for Vee
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
            !srweight_vee(sm)=weight_vee(sm) !for Vee
        end do
    end do
    deallocate(gr1)
    deallocate(gr2)    
end if   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (nquad.gt.1) then         !we have more than one quadrature centre
    !write(*,*) "number of points per centre", smn(:)
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
sm=0
!write(*,*) "Final number of grid points=", rrgrid
do i=1,rgrid
    if ((srweight(i).gt.trsh2)) then !remove points with low (zero) weight
        sm=sm+1 
        rweight(sm)=srweight(i) 
        !rweight_vee(sm)=srweight_vee(i) !for Vee
        actual_dist=dsqrt(sum(gr3(:,i)**2)) !compute actual distance of the point to the origin
        rweight_vee(sm)=srweight(i)*(actual_dist)**(-1.d0) !for Vee, divide by r to get correct weight for Vee
        rrrg(:,sm)=gr3(:,i)           
    end if
end do
write(*,*) "Grid points:"
open(4,file='gridpoints.dat',status='replace')
do i=1,rrgrid
    write(4,'(I6,3F15.8,F30.8)') i, rrrg(1,i), rrrg(2,i), rrrg(3,i), rweight(i)
end do
deallocate(weight); deallocate(srweight)
deallocate(weight_vee); deallocate(srweight_vee)
deallocate(gr3)
close(4)
end subroutine 

