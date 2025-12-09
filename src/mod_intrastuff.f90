module intrastuff 
use geninfo
use quadratures
use numbers
implicit none
!global variables that depens on all ijkl primitives
double precision :: alfijkl, zeta, eijkl
double precision :: invz, sqe
integer, allocatable, dimension(:) :: ipiv
!primtive parameters
!double precision :: Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,Xl,Yl,Zl !center of each of the primitives
!double precision :: ti,mi,ni,tj,mj,nj,tk,mk,nk,tl,ml,nl !angular momenta of each primitive
!double precision :: alfi,alfj,alfk,alfl !exponents of each primitive
!Cioslowski's parameters
!double precision :: sqe !dsqrt(e_ijkl) save time and compute it once
!Grid independent variables for the intracule:
!double precision ::  zeta 
!double precision ::  eijkl 
!double precision, dimension(3) :: R_ik, R_jl, R_ijkl !X_ik, ...
!double precision :: alfijkl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!Functions for the first integral screening!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        function JA8(ti,tk,mi,mk,ni,nk,alfi,alfk,Xi,Xk,Yi,Yk,Zi,Zk,Rik2,aik,eik) !J upper bound for the first integral screening
        implicit none                    !eqn A8 
        double precision :: JA8          
        integer, intent(in) :: ti,tk,mi,mk,ni,nk
        double precision, intent(in) :: alfi,alfk
        double precision, intent(in) :: Xi,Xk,Yi,Yk,Zi,Zk
        double precision, intent(in) :: Rik2
        double precision, intent(in) :: aik, eik
        !local variables
        double precision, dimension(3) :: x_max
        integer :: n_nodes
        !x_max values depending on the number of nodes
        n_nodes=ti+tk+1
        if (n_nodes.eq.1) then
                x_max(1)=ZERO
        else if (n_nodes.eq.2) then
                x_max(1)=0.707106781186548d0 
        else if (n_nodes.eq.3) then 
                x_max(1)=1.224744871391589d0 
        else if (n_nodes.eq.4) then
                x_max(1)=1.650680123885785d0
        else if (n_nodes.eq.5) then 
                x_max(1)=2.020182870456086d0 
        else if (n_nodes.eq.6) then 
                x_max(1)=2.350604973674492d0
        else if (n_nodes.eq.7) then 
                x_max(1)=2.651961356835233d0 
        else if (n_nodes.eq.8) then 
                x_max(1)=2.930637420257244d0
        else if (n_nodes.eq.9) then 
                x_max(1)=3.190993201781528d0  
        else if (n_nodes.eq.10) then 
                x_max(1)=3.436159118837738d0
        end if
        !y
        n_nodes=mi+mk+1
        if (n_nodes.eq.1) then
                x_max(2)=ZERO
        else if (n_nodes.eq.2) then
                x_max(2)=0.707106781186548d0 
        else if (n_nodes.eq.3) then 
                x_max(2)=1.224744871391589d0 
        else if (n_nodes.eq.4) then
                x_max(2)=1.650680123885785d0
        else if (n_nodes.eq.5) then 
                x_max(2)=2.020182870456086d0 
        else if (n_nodes.eq.6) then 
                x_max(2)=2.350604973674492d0
        else if (n_nodes.eq.7) then 
                x_max(2)=2.651961356835233d0 
        else if (n_nodes.eq.8) then 
                x_max(2)=2.930637420257244d0
        else if (n_nodes.eq.9) then 
                x_max(2)=3.190993201781528d0  
        else if (n_nodes.eq.10) then 
                x_max(2)=3.436159118837738d0
        end if
        !z
        n_nodes=ni+nk+1
        if (n_nodes.eq.1) then
                x_max(3)=ZERO
        else if (n_nodes.eq.2) then
                x_max(3)=0.707106781186548d0 
        else if (n_nodes.eq.3) then 
                x_max(3)=1.224744871391589d0 
        else if (n_nodes.eq.4) then
                x_max(3)=1.650680123885785d0
        else if (n_nodes.eq.5) then 
                x_max(3)=2.020182870456086d0 
        else if (n_nodes.eq.6) then 
                x_max(3)=2.350604973674492d0
        else if (n_nodes.eq.7) then 
                x_max(3)=2.651961356835233d0 
        else if (n_nodes.eq.8) then 
                x_max(3)=2.930637420257244d0
        else if (n_nodes.eq.9) then 
                x_max(3)=3.190993201781528d0  
        else if (n_nodes.eq.10) then 
                x_max(3)=3.436159118837738d0
        end if
        JA8=pi**(ONEANDHALF)*(TWO*aik)**(-(dble(ti+tk+mi+mk+ni+nk)+ONEANDHALF))&
        *dexp(-TWO*eik*(Rik2))&
        *((x_max(1)+alfk*dsqrt(TWO/aik)*dabs(Xk-Xi))**(TWO*dble(ti)))&
        *((x_max(1)+alfi*dsqrt(TWO/aik)*dabs(Xk-Xi))**(TWO*dble(tk)))&
        *((x_max(2)+alfk*dsqrt(TWO/aik)*dabs(Yk-Yi))**(TWO*dble(mi)))&
        *((x_max(2)+alfi*dsqrt(TWO/aik)*dabs(Yk-Yi))**(TWO*dble(mk)))&
        *((x_max(3)+alfk*dsqrt(TWO/aik)*dabs(Zk-Zi))**(TWO*dble(ni)))&
        *((x_max(3)+alfi*dsqrt(TWO/aik)*dabs(Zk-Zi))**(TWO*dble(nk)))
        end function


        !JA8=pi**(ONEANDHALF)*(TWO*a_ik)**(-(dble(TMN(i,1)+TMN(k,1)+TMN(i,2)+TMN(k,2)+TMN(i,3)+TMN(k,3))+ONEANDHALF))&
        !*dexp(-TWO*e_ik*(Rik2))&
        !*((x_max(1)+Alpha(k)*dsqrt(TWO/a_ik)*dabs(Cartes(Ra(k),1)-Cartes(Ra(i),1)))**(TWO*dble(TMN(i,1))))&
        !*((x_max(1)+Alpha(i)*dsqrt(TWO/a_ik)*dabs(Cartes(Ra(k),1)-Cartes(Ra(i),1)))**(TWO*dble(TMN(k,1))))&
        !*((x_max(2)+Alpha(k)*dsqrt(TWO/a_ik)*dabs(Cartes(Ra(k),2)-Cartes(Ra(i),2)))**(TWO*dble(TMN(i,2))))&
        !*((x_max(2)+Alpha(i)*dsqrt(TWO/a_ik)*dabs(Cartes(Ra(k),2)-Cartes(Ra(i),2)))**(TWO*dble(TMN(k,2))))&
        !*((x_max(3)+Alpha(k)*dsqrt(TWO/a_ik)*dabs(Cartes(Ra(k),3)-Cartes(Ra(i),3)))**(TWO*dble(TMN(i,3))))&
        !*((x_max(3)+Alpha(i)*dsqrt(TWO/a_ik)*dabs(Cartes(Ra(k),3)-Cartes(Ra(i),3)))**(TWO*dble(TMN(k,3))))
        !end function
 
       ! function J_jl(i,l,R_j_l_2) !J upper bound for the first integral screening
       ! double precision :: J_jl
       ! integer, intent(in) :: j,l
       ! double precision, intent(in) :: R_j_l_2
       ! double precision, dimension(3) :: x_max 
       ! integer :: n_nodes, ii
       ! n_nodes=0  
       ! do ii=1,3 !loop for J_ik(x_i(*), (y) and (z)
       !         n_nodes=TMN(j,ii)+TMN(l,ii) + 1
       !         if (n_nodes.eq.1) then 
       !                 x_max(ii)=ZERO
       !         else if (n_nodes.eq.2) then
       !                 x_max(ii)=0.707106781186548d0 
       !         else if (n_nodes.eq.3) then 
       !                 x_max(ii)=1.224744871391589d0 
       !         else if (n_nodes.eq.4) then 
       !                 x_max(ii)=1.650680123885785d0
       !         else if (n_nodes.eq.5) then 
       !                 x_max(ii)=2.020182870456086d0 
       !         else if (n_nodes.eq.6) then 
       !                 x_max(ii)=2.350604973674492d0
       !         else if (n_nodes.eq.7) then 
       !                 x_max(ii)=2.651961356835233d0 
       !         else if (n_nodes.eq.8) then 
       !                 x_max(ii)=2.930637420257244d0
       !         else if (n_nodes.eq.9) then 
       !                 x_max(ii)=3.190993201781528d0  
       !         else if (n_nodes.eq.10) then 
       !                 x_max(ii)=3.436159118837738d0
       !         end if
       ! end do
  
        !J_jl=pi**(ONEANDHALF)*(TWO*a_jl)**(-(dble(TMN(j,1)+TMN(l,1)+TMN(j,2)+TMN(l,2)+TMN(j,3)+TMN(l,3))+ONEANDHALF))&
        !*dexp(-TWO*e_jl*(R_j_l_2))&
        !*((x_max(1)+Alpha(l)*dsqrt(TWO/a_jl)*dabs(Cartes(Ra(l),1)-Cartes(Ra(j),1)))**(TWO*dble(TMN(j,1))))&
        !*((x_max(1)+Alpha(j)*dsqrt(TWO/a_jl)*dabs(Cartes(Ra(l),1)-Cartes(Ra(j),1)))**(TWO*dble(TMN(l,1))))&
        !*((x_max(2)+Alpha(l)*dsqrt(TWO/a_jl)*dabs(Cartes(Ra(l),2)-Cartes(Ra(j),2)))**(TWO*dble(TMN(j,2))))&
        !*((x_max(2)+Alpha(j)*dsqrt(TWO/a_jl)*dabs(Cartes(Ra(l),2)-Cartes(Ra(j),2)))**(TWO*dble(TMN(l,2))))&
        !*((x_max(3)+Alpha(l)*dsqrt(TWO/a_jl)*dabs(Cartes(Ra(l),3)-Cartes(Ra(j),3)))**(TWO*dble(TMN(j,3))))&
        !*((x_max(3)+Alpha(j)*dsqrt(TWO/a_jl)*dabs(Cartes(Ra(l),3)-Cartes(Ra(j),3)))**(TWO*dble(TMN(l,3))))
        !end function
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!eq 15!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!function Wijkl(i,j,k,l,rhat,r,ax)  !for 1st integral screening. Equation 15
        subroutine polycoef(Ri,Rj,Rk,Rl,Rik,Rjl,Rijkl,tmni,tmnj,tmnk,tmnl,np,n,rh,wh,Cijkl)
                implicit none
                !==== Arguments ====
                double precision, intent(in) :: Ri,Rj,Rk,Rl !primitive centers
                integer, intent(in) :: tmni, tmnj, tmnk, tmnl !angular momenta
                double precision, intent(in) :: Rik,Rjl,Rijkl
                integer, intent(in) :: np, n !np is the number of points, n is 
                double precision, intent(out) :: Cijkl(:) !coefficients of the polynomial
                double precision, intent(in) :: wh(:), rh(:)
                !integer, intent(inout) :: ipiv(:)   
                !==== Local variables ====
                double precision, dimension(np,np) :: M  
                double precision, dimension(np)    :: xx
                integer :: Ltot
                integer :: i1, j1, l1, INFO            
                if (np == 1) then
                    Cijkl = wh(1)
                else
                    Ltot = np - 1
                    ! Generate points
                    do i1 = 1, np
                        xx(i1) = -HALF + dble(i1)
                    end do
            
                    ! Build M
                    do i1 = 1, np
                        l1 = Ltot
                        do j1 = 1, np
                            M(i1,j1) = (sqe*(xx(i1)+Rik-Rjl))**dble(l1)
                            l1 = l1 - 1
                        end do
                    end do
            
                    ! Build RHS vector Cijkl
                    Cijkl = 0.0d0
                    do i1 = 1, np
                        do j1 = 1, n
                            Cijkl(i1) = Cijkl(i1) + wh(j1)*Wijkl(Ri,Rj,Rk,Rl,Rijkl,tmni,tmnj,tmnk,tmnl,rh(j1),xx(i1))
                        end do
                    end do
            
                    ! Solve
                    call dgesv(np,1,M,np,ipiv,Cijkl,np,INFO)
                    if (INFO /= 0) write(*,*) "ERROR, cannot solve linear equations"
                end if
            
            contains
            
                !===========================================================
                ! Internal Wijkl function
                !===========================================================
                function Wijkl(Ri,Rj,Rk,Rl,Rijkl,tmni,tmnj,tmnk,tmnl,rhat,r) result(Wv)
                    implicit none
                    double precision, intent(in) :: Ri,Rj,Rk,Rl !primitive centers
                    integer, intent(in) :: tmni, tmnj, tmnk, tmnl !angular momenta
                    double precision, intent(in) :: rijkl
                    double precision :: Wv
                    double precision, intent(in) :: rhat, r
                    double precision :: invzrhat, alfijklplushalfr, alfijklminushalfr
                    invzrhat=invz*rhat
                    alfijklplushalfr=(alfijkl+HALF)*r 
                    alfijklminushalfr=(alfijkl-HALF)*r
                    Wv=(invzrhat+alfijklminushalfr+(rijkl-Ri))**dble(tmni)* &
                         (invzrhat+alfijklplushalfr+(rijkl-Rj))**dble(tmnj)* &
                         (invzrhat+alfijklminushalfr+(rijkl-Rk))**dble(tmnk)* &
                         (invzrhat+alfijklplushalfr+(rijkl-Rl))**dble(tmnl)
                    !Wv = ((dsqrt(zeta)**(-ONE)*rhat)+(alfijkl-HALF)*r+(r_ijkl(ax)-Cartes(Ra(i),ax)))**dble(tmni) * &
                    !     ((dsqrt(zeta)**(-ONE)*rhat)+(alfijkl+HALF)*r+(r_ijkl(ax)-Cartes(Ra(j),ax)))**dble(tmnj) * &
                    !     ((dsqrt(zeta)**(-ONE)*rhat)+(alfijkl-HALF)*r+(r_ijkl(ax)-Cartes(Ra(k),ax)))**dble(tmnk) * &
                    !     ((dsqrt(zeta)**(-ONE)*rhat)+(alfijkl+HALF)*r+(r_ijkl(ax)-Cartes(Ra(l),ax)))**dble(tmnl)
                end function Wijkl
            
            end subroutine polycoef

subroutine gauherm(Lrtot,n,rh,wh) 
        !performs gauss-hermite quadrature                          
        !we will obtain the nodes and weights depending on the degree of the polinomial (2n-1) --> (n) 
        !(n is nx ny or nz in intracule.f90)
        implicit none
        integer, intent(in) :: Lrtot !degree of polynomial
        double precision,intent(out) :: rh(:), wh(:) !nodes
        integer, intent(out) :: n    !number of gauss-hermite nodes
        !local variables
        integer :: factn, i2
        double precision :: Cons
        !double precision, parameter :: pi=FOUR*datan(ONE)
        if (Lrtot.gt.0) then
        !compute the number of nodes for exact Hermite quadrature
             if(MOD(Lrtot,2).eq.0) then !even 
                   n=int((dble(Lrtot)*HALF)+ONE) !nx->number of nodes
             else 
                   n=int((dble(Lrtot)+ONE)*HALF) 
             end if
         else
          !only one node
            n=1
         end if  
         !store rh values depending on pol. degree
          if (n.eq.1) then !1 node-->degree of the pol is 0 or 1 
           rh(1)=ZERO 
          else if (n.eq.2) then !2 nodes-->degree of the pol is 3 or 2                
           rh(2)=0.7071067811865475d0 
           rh(1)=-0.7071067811865475d0 
          else if (n.eq.3) then !3 nodes-->degree of pol is 5 or 4
           rh(2)=ZERO
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
           rh(4)=ZERO
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
           rh(5)=ZERO
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
            rh(6)=ZERO
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
          Cons=((TWO)**(dble(n-1)) * dble(factn)* dsqrt(pi))/dble(n*n)
          do i2=1,n
           wh(i2)=cons*(ONE/(Hermite(rh(i2),n))**TWO)  
          end do  
          contains
            function Hermite(x,n) !
            double precision :: Hermite
            integer, intent(in) :: n !numer of nodes of pol
            double precision, intent(in) :: x !roots of Hermite pol
           
           if (n.eq.1) then !0 deg
             Hermite=ONE
           else if (n.eq.2) then
             Hermite=TWO*x 
           else if (n.eq.3) then
             Hermite=FOUR*(x**TWO) -TWO
           else if (n.eq.4) then
             Hermite=EIGHT*(x**THREE) -12.d0*x
           else if (n.eq.5) then
             Hermite=16.d0*(x**FOUR)-48.d0*(x**TWO)+12.d0  
           else if (n.eq.6) then
             Hermite=32.d0*(x**FIVE)-160.d0*(x**THREE)+120.d0*x  
           else if (n.eq.7) then
             Hermite=64.d0*(x**SIX)-480.d0*(x**FOUR)+720.d0*(x**TWO)-120.d0
           else if (n.eq.8) then
             Hermite=128.d0*(x**SEVEN)-1344.d0*(x**FIVE)+3360.d0*(x**THREE)-1680.d0*x    
           else if (n.eq.9) then
             Hermite=256.d0*(x**EIGHT)-3584.d0*(x**SIX)+13440.d0*(x**FOUR)-13440.d0*(x**TWO)+1680.d0  
           else if (n.eq.10) then
             Hermite=512.d0*(x**NINE)-9216.d0*(x**SEVEN)+48384.d0*(x**FIVE)-80640.d0*(x**THREE)+30240.d0*x
           else if (n.eq.11) then
           Hermite=1024.d0*(x**TEN)-23040.d0*(x**EIGHT)+161280.d0*(x**SIX)-403200.d0*(x**FOUR)+302400.d0*(x**TWO)-30240.d0  
           end if       
           end function
        
        end subroutine gauherm 
      
end module intrastuff
