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
implicit none
logical, intent(in) :: normalize_dm2p
!!!!!!!!!!!!!!!local variables from the paper!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: kk1, kk2 !integers we do not want to read from .dm2
integer :: dfact    !double factorial
integer :: ii,iii,sum, ig, ir, summ, np, sm, smm !some integers
double precision :: Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,Xl,Yl,Zl !primitive centers
double precision :: alfi,alfj,alfk,alfl !primitive exponents
integer :: ti,mi,ni,tj,mj,nj,tk,mk,nk,tl,ml,nl !angular momenta of primitives
double precision :: aik, ajl, eik, ejl            !variables from eqn 10
double precision ::  Xik, Yik, Zik, Xjl, Yjl, Zjl !variables from eqn 10
!double precision :: zeta, eijkl, alfijkl,sqe      !variables from eqn 12 in module intrastuff
double precision :: Xijkl, Yijkl, Zijkl          !variables from eqn 12
double precision :: Rik2, Rjl2  !R_ik=(R_i-R_k)^2 !for equation A8
double precision :: Jik, Jjl !J of equation A8 for screening
double precision :: Aijkl !grid independent part of the intracule (eqn 18)
double precision, allocatable, dimension(:) :: rh, wh !nodes and weights for gauss hermite (eqn16)
integer :: nn !number of nodes for gauss hermite (eqn 16) output from gauherm, input to polycoef
integer, parameter :: Lmax=21 !maximum angular momentum for primitives (eq 16)
integer, parameter :: nmax = (Lmax+1)/2 + 1
integer, parameter :: npmax = Lmax + 1
double precision, allocatable, dimension(:) :: Cx, Cy, Cz !coefficients of V (eqns 16,42) 
double precision, allocatable, dimension(:) :: wmarray !eq A10, for 2nd integral screening (eq 43)
double precision :: Ux, Uy, Uz !coef. for 2nd integral screening (eqn. 43)
double precision :: Vx, Vy, Vz !eq 16 
!!!!!!!!!!primitive quartet variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: i, j, k, l !indices for primitive functions
double precision :: nprimt !total normalization factor
double precision :: DMval !value of DM2
double precision :: screen1, screen2, lim !integral screenings
integer, dimension(3) :: Lrtot  !total angular momentum and number of nodes
!integer, allocatable, dimension(:) :: ipiv !for dgesv in intrastuff
!GRID POINTS
integer :: ngrid !number of grid points
double precision, allocatable, dimension(:,:) :: r !grid points
double precision, dimension(3) :: rp  !prima points(X', Y', Z')
!intracules
double precision :: rintegral                         !radial integral
double precision, allocatable, dimension(:) :: rintra !radial intracule
double precision, allocatable, dimension(:) :: Ivec   !intracule at a point
!check accuracy of calculations
double precision :: traceDM2prim, trDM2 !normalized and not normalized DM2prim
integer :: npairs      !number of electron pairs
double precision :: T1, T2, T3, T4, TT1, TT2, TT3, TT4, TT5 !time check
double precision :: Tread, T1screen, T2screen, Tgrid
double precision :: vee
double precision :: intraculezero
!primitive normalization
double precision, allocatable, dimension(:) :: normprim
!set counter for primitive quartets
integer :: quartetcount, refval
!block summation stuff
!integer, parameter :: BUFSIZE = 128   ! or 128 for even better stability
!double precision, allocatable :: buf(:,:)
!integer, allocatable :: bufcnt(:)


call cpu_time(T1)
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
!allocate(buf(BUFSIZE, ngrid))
!allocate(bufcnt(nGrid))
!buf = ZERO
!bufcnt = 0
allocate(Ivec(ngrid))
Ivec=ZERO
call cpu_time(T2)    
traceDM2prim=ZERO !sum of all the normalized DM2 terms.
!Time check
Tread=ZERO
T1screen=ZERO
T2screen=ZERO
Tgrid=ZERO
allocate(rh(nmax)) 
rh=ZERO
allocate(wh(nmax))
allocate(Cx(npmax)); allocate(Cy(npmax)); allocate(Cz(npmax))
allocate(ipiv(npmax))
!!!!!!!Precompute w_m values for 2nd integral screening!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(wmarray(npmax))
do ii = 1, npmax
    if (ii == 1) then
        wmarray(ii) = ONE
    else
        wmarray(ii) = (dble(ii) * HALF)**(dble(ii) * HALF) * dexp(-dble(ii) * HALF)
    end if
end do
!perform primitive normalization for DM2 if requested before primitive loop
allocate(normprim(nprim))
if (normalize_dm2p) then
    do ii=1,nprim
        normprim(ii)=(TWO*Alpha(ii)/pi)**(THREEQUARTER)&
        *dsqrt(((FOUR*Alpha(ii))**(dble(TMN(ii,1)+TMN(ii,2)+TMN(ii,3))))*&
        dble(dfact(2*TMN(ii,1)-1)*dfact(2*TMN(ii,2)-1)*dfact(2*TMN(ii,3)-1))**(-ONE))  
    end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!Start loop over primitive quartets!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!Following Cioslowski and Liu algorithm!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
write(*,*) "Starting loop over primitive quartets..."
quartetcount=0
lim=thresh*(dble(nprim)*(dble(nprim)+ONE)*HALF)**(-ONE) !limit for the 1st integral screening
open(unit=5,file=dm2name, form='unformatted',access='stream') !open binary file
do while (.true.)  !loop for primitive quartets.
    call cpu_time(TT1)    
    !read(5,end=100, err=200) kk1,i,j,k,l,DMval,kk2!read a line from binary file .dm2
    read(5,end=100, err=200) i,j,k,l,DMval
    if (i.lt.1 .or. i.gt.nprim) goto 200
    if (j.lt.1 .or. k.lt.1 .or. l.lt.1) goto 200
    quartetcount=quartetcount+1
    !if (i.eq.0) goto 100 !file is finished            
    call cpu_time(TT2) 
    Tread=Tread+(TT2-TT1)
    trDM2=trDM2+DMval
    smm=smm+1
    if (normalize_dm2p) then
        nprimt=normprim(i)*normprim(j)*normprim(k)*normprim(l)
        DMval=nprimt*DMval 
    end if     
    traceDM2prim=traceDM2prim+DMval !sum all the DM2 quartets to check accuracy
    !load data for i,j,k,l primitives (only one for speed-up)
    alfi=Alpha(i); alfj=Alpha(j); alfk=Alpha(k); alfl=Alpha(l)   
    ti=TMN(i,1); mi=TMN(i,2); ni=TMN(i,3)
    tj=TMN(j,1); mj=TMN(j,2); nj=TMN(j,3)
    tk=TMN(k,1); mk=TMN(k,2); nk=TMN(k,3)
    tl=TMN(l,1); ml=TMN(l,2); nl=TMN(l,3)
    Xi=Xn(i); Yi=Yn(i); Zi=Zn(i)
    Xj=Xn(j); Yj=Yn(j); Zj=Zn(j)
    Xk=Xn(k); Yk=Yn(k); Zk=Zn(k)    
    Xl=Xn(l); Yl=Yn(l); Zl=Zn(l)
    !compute the first variables
    aik=alfi+alfk                  !eqn. 10
    ajl=alfj+alfl                  !eqn. 10                         
    eik=alfi*alfk*aik**(-ONE)
    ejl=alfj*alfl*ajl**(-ONE)   
    if ((aik.lt.1.d-4).or.(ajl.lt.1.d-4)) then !skip this quartet if exponents are too small
        write(*,*) "one exponent too small", aik, ajl, i, j, k, l
    end if
    Rik2=(Xi-Xk)**TWO+(Yi-Yk)**TWO+(Zi-Zk)**TWO    
    Rjl2=(Xj-Xl)**TWO+(Yj-Yl)**TWO+(Zj-Zl)**TWO                
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1st integral screening!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Jik=JA8(ti,tk,mi,mk,ni,nk,alfi,alfk,xi,xk,yi,yk,zi,zk,Rik2,aik,eik)
    Jjl=JA8(tj,tl,mj,ml,nj,nl,alfj,alfl,xj,xl,yj,yl,zj,zl,Rjl2,ajl,ejl)                                   
    screen1=dabs(DMval)*dsqrt(Jik*Jjl) !eqn.39           
    call cpu_time(TT3)
    T1screen=T1screen-(TT3-TT2)         
    if (screen1.ge.lim) then    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
        !compute the other variables (eqn. 12)
        zeta=aik+ajl
        !and also careful with this!
        if (zeta.lt.1.d-4) then
            write(*,*) "zeta_ijkl too small", zeta, i, j, k, l
        end if
        eijkl=(aik*ajl)*zeta**(-ONE)
        sqe=dsqrt(eijkl)
        !careful with these expressions!
        Xik=(alfi*Xi+alfk*Xk)*(aik**(-ONE)) !eq. 11
        Yik=(alfi*Yi+alfk*Yk)*(aik**(-ONE))
        Zik=(alfi*Zi+alfk*Zk)*(aik**(-ONE))
        Xjl=(alfj*Xj+alfl*Xl)*(ajl**(-ONE))
        Yjl=(alfj*Yj+alfl*Yl)*(ajl**(-ONE))
        Zjl=(alfj*Zj+alfl*Zl)*(ajl**(-ONE))
        Xijkl=(aik*Xik+ajl*Xjl)*zeta**(-ONE) !eq. 12
        Yijkl=(aik*Yik+ajl*Yjl)*zeta**(-ONE)
        Zijkl=(aik*Zik+ajl*Zjl)*zeta**(-ONE)
        !do ii=1,3                   !X_ik,Y_ik,Z_ik,...
        !    R_ik(ii)=(Alpha(i)*Cartes(Ra(i),ii)+alfk*Cartes(Ra(k),ii))*aik**(-ONE)
        !    R_jl(ii)=(Alpha(j)*Cartes(Ra(j),ii)+Alpha(l)*Cartes(Ra(l),ii))*ajl**(-ONE)
        !    R_ijkl(ii)=(aik*R_ik(ii)+ajl*R_jl(ii))*zeta**(-ONE)
        !end do
        invz=dsqrt(zeta)**(-ONE) !zeta^(-1)                     
        alfijkl=HALF*(aik-ajl)*zeta**(-ONE)
        !compute Aijkl (the grid independent part)  !eq.18
        Aijkl=DMval*(zeta)**(-ONEANDHALF)*dexp(-eik*Rik2-ejl*Rjl2) 
        !!!!!!!!!!!!Calculate coeficients of V!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        sm=0
        ! --- x direction ---
        !Lrtot(1)=TMN(i,1)+TMN(j,1)+TMN(k,1)+TMN(l,1) !Lrtot=degree of eqn 15 (I don't use Lmax)
        Lrtot(1)=ti+tj+tk+tl 
        call gauherm(Lrtot(1),nn,rh,wh) !obtain Gauss-Hermite nodes(rh) and weights(wh) (2L+1)
        np=Lrtot(1)+1              !np is the number of coefficients to represent the polynomial V
        call polycoef(Xi,Xj,Xk,Xl,Xik,Xjl,Xijkl,ti,tj,tk,tl,np,nn,rh,wh,Cx) !obtain coefficients of V in x direction
        Ux=dot_product(Cx(1:np),wmarray(1:np)) !compute U coefficients for 2nd integral screening (eq.43 of the paper)
        ! --- y direction ---
        !Lrtot(2)=TMN(i,2)+TMN(j,2)+TMN(k,2)+TMN(l,2) !Lrtot=degree of eqn 15 (I don't use Lmax)
        Lrtot(2)=mi+mj+mk+ml
        call gauherm(Lrtot(2), nn, rh, wh) !obtain Gauss-Hermite nodes(rh) and weights(wh) (2L+1)
        np=Lrtot(2)+1              !np is the number of coefficients to represent the polynomial V
        call polycoef(Yi,Yj,Yk,Yl,Yik,Yjl,Yijkl,mi,mj,mk,ml,np,nn,rh,wh,Cy) !obtain coefficients of V in y direction
        Uy=dot_product(Cy(1:np),wmarray(1:np)) !compute U coefficients for 2nd integral screening (eq.43 of the paper)  
        ! --- z direction ---
        Lrtot(3)=ni+nj+nk+nl

        !Lrtot(3)=TMN(i,3)+TMN(j,3)+TMN(k,3)+TMN(l,3) !Lrtot=degree of eqn 15 (I don't use Lmax)
        call gauherm(Lrtot(3), nn, rh, wh) !obtain Gauss-Hermite nodes(rh) and weights(wh) (2L+1)
        np=Lrtot(3)+1              !np is the number of coefficients to represent the polynomial V
        call polycoef(Zi,Zj,Zk,Zl,Zik,Zjl,Zijkl,ni,nj,nk,nl,np,nn,rh,wh,Cz) !obtain coefficients of V in z direction
        Uz=dot_product(Cz(1:np),wmarray(1:np)) !compute U coefficients for 2nd integral screening (eq.43 of the paper)             
        !!!!!!!!!!2nd Integral Screening!!!!!!!!!!!!!!!!!!!!!                      
        screen2=dabs(Aijkl*Ux*Uy*Uz) !eqn.40
        call CPU_time(TT4)
        T2screen=T2screen+(TT4-TT3)                
        if (screen2.ge.lim) then 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                             
            !loop over grid points
            do ig=1,nGrid                                             
                Vx=ZERO; Vy=ZERO; Vz=ZERO  
                !rp(1)=sqe*(r(1,ig)+r_ik(1)-r_jl(1)) !compute R'(eq. 19)
                !rp(2)=sqe*(r(2,ig)+r_ik(2)-r_jl(2))
                !rp(3)=sqe*(r(3,ig)+r_ik(3)-r_jl(3)) 
                rp(1)=sqe*(r(1,ig)+Xik-Xjl) !compute R'(eq. 19)
                rp(2)=sqe*(r(2,ig)+Yik-Yjl)
                rp(3)=sqe*(r(3,ig)+Zik-Zjl)
                ! ----Horner's method for x ----
                do iii=1,Lrtot(1)+1   !Ltot+1 is the number of coefficients we have
                    Vx=Vx*rp(1)+Cx(iii)
                end do  !end loop over the polynomial coefficients
                ! ----Horner's method for y ----
                do iii=1,Lrtot(2)+1   !Ltot+1 is the number of coefficients we have
                    Vy=Vy*rp(2)+Cy(iii)
                end do  !end loop over the polynomial coefficients
                ! ----Horner's method for z ----
                do iii=1,Lrtot(3)+1   !Ltot+1 is the number of coefficients we have
                    Vz=Vz*rp(3)+Cz(iii)
                end do  !end loop over the polynomial coefficients
                !!!!!!!!!!!!!Calculate intracule at the grid point!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                Ivec(ig)=Ivec(ig)+Aijkl*dexp(-(rp(1)**TWO+rp(2)**TWO+rp(3)**TWO))*Vx*Vy*Vz
                !term = Aijkl * dexp(-(rp(1)**TWO + rp(2)**TWO + rp(3)**TWO)) * Vx * Vy * Vz
                ! ---- Block summation ----
                !cnt = bufcnt(ig) + 1
                !bufcnt(ig) = cnt
                !buf(cnt,ig) = term
                !if (cnt == BUFSIZE) then
                !    Ivec(ig) = Ivec(ig) + pairwise_sum(buf(1:BUFSIZE,ig), BUFSIZE)
                !    bufcnt(ig) = 0
                !end if
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
write(*,*) "Reached end of .dm2p file."
write(*,*) "Total number of quartets read:", quartetcount
goto 300
200 continue
write(*,*) "REad error or corrupted record near quartet", quartetcount+1
write(*,*) "Ended loop for primitives"
write(*,*) "Stopping safely"
goto 300

300 continue

write(*,*) "Loop over primtives completed"
deallocate(Cx)
deallocate(Cy)   !deallocate polynomial coefficients
deallocate(Cz) 
!check if we read a symmetry reduced dmp2prim or the full dm2prim
refval=(((nprim*(nprim+1)/2)+1)*(nprim*(nprim+1)/2))/2
write(*,*) "last quartet", i,j,k,l
write(*,*) "nprim", nprim, "refval", refval, "quartet_count", quartetcount
if (quartetcount.gt.refval) then
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
    allocate(rintra(nradi)) 
    ig=0
    sm=0
    rintra=ZERO
    ir=0
    do i=1,nradi   
        write(*,*) i, radi(i), smn(i)
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
write(*,*) "Sum of all DM2prim terms=", trDM2, traceDM2prim
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
!contains
!recursive function pairwise_sum(arr, n) result(res)
!    implicit none
!    double precision, intent(in) :: arr(n)
!    integer, intent(in) :: n
!    double precision :: res
!    integer :: m, i
!    double precision, allocatable :: tmp(:)
!    if (n == 1) then
!        res = arr(1)
!        return
!    end if
!    m = (n + 1) / 2
!    allocate(tmp(m))
!    do i = 1, n-1, 2
!        tmp((i+1)/2) = arr(i) + arr(i+1)
!    end do
!    if (mod(n,2) == 1) tmp(m) = arr(n)
!    res = pairwise_sum(tmp, m)
!    deallocate(tmp)
!end function pairwise_sum
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!