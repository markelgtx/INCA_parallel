subroutine c1hole(nrad,sfalpha)
    use wfxinfo
    use geninfo
    implicit none
    !radial quadrature!!!!!!!!
    integer, intent(in) :: nrad
    double precision, intent(in) :: sfalpha
    !Local variables
    double precision, dimension(nrad) :: r, rad, wr
    !angular quadrature!!!!!!!
    double precision, dimension(1202) :: wa, lbx, lby, lbz, phi,theta
    double precision, dimension(1202) :: cosin, sinsin, cost
    integer :: nang
    double precision :: xs
    !grid points!!!!!!!!!!!!!!
    double precision :: px,py,pz
    double precision :: x,y,z
    !functions!!!!!!!!!!!!!!!!
    double precision :: rDM1, beckex, Density, rdm1_alf
    double precision, dimension(45) :: hc1
    double precision :: dens, xchole
    !!!!!!!!!!!!!!!!!!!!!!!!!/
    double precision :: rdmv1, rho2, rho2s, rdmv2 !vector 1rdm and pair dens
    integer :: i,j,k,s
    double precision :: kk, xcint, xcint2, radint, radint2!densint, rdm1int, radint2
    double precision :: b,alf
    character*40 :: outputname
    integer :: i1, npoints
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!obtain radial nodes and weights for integrations!!!!!!!!!!!

    call sub_GauLeg(-1.d0,1.d0,rad,wr,nrad)
    do i=1,nrad
         r(i)=(1.d0+rad(i))/(1.d0-rad(i))*sfalpha
         wr(i)=(2.d0*sfalpha/((1.d0-rad(i))**2.d0))*wr(i)
    end do

    !obtain angular nodes (angles) and weights
    lbx=0.d0
    lby=0.d0 
    lbz=0.d0
    nang=1202
    !call LD0006(lbx,lby,lbz,Wa,nAng)  !6 grid points
    !call LD0110(lbx,lby,lbz,wa,nang)
    !call LD0590(lbx,lby,lbz,wa,nAng) !590
    call LD1202(lbx,lby,lbz,wa,nAng)

    do i=1,nang
         theta(i)=dacos(lbz(i))
         if (dsin(theta(i)).ne.0.d0) then
              xs=lbx(i)/dsin(theta(i))
              if (xs.gt.1.d0) xs=1.d0
              if (xs.lt.-1.d0) xs=-1.d0
              phi(i)=dacos(xs)
              if (lby(i).lt.0.d0) phi(i)=-phi(i)
         else
              phi(i)=0.d0
         end if   
         wa(i)=wa(i)*4.d0*pi 
    end do

    open(unit=3, file="angular_grid")
    do i=1,nang
           x=dcos(phi(i))*dsin(theta(i))
           y=dsin(phi(i))*dsin(theta(i))
           z=dcos(theta(i))
           write(3,*) x,y,z
           cosin(i)=x
           sinsin(i)=y
           cost(i)=z
    end do
    close(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!Sanity checks!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Test 1: integrate the density and 1rdm diagonal (must be = nelec)
    radint=0.d0
    radint2=0.d0
    !do i=1,nrad
    !   densint=0.d0
    !   rdm1int=0.d0
    !   do k=1,nang
    !       x=r(i)*dcos(phi(k))*dsin(theta(k))
    !       y=r(i)*dsin(phi(k))*dsin(theta(k))
    !       z=r(i)*dcos(theta(k))
    !       densint=densint+wa(k)*Density(x,y,z)
    !       rdm1int=rdm1int+wa(k)*rDM1(x,y,z,x,y,z)
    !   end do    
    !   radint=radint+wr(i)*densint*r(i)**2.d0   
    !   radint2=radint2+wr(i)*rdm1int*r(i)**2.d0
   ! end do  
   ! write(*,*) "----------Test 1-------------" 
   ! write(*,*) "integral of density=", radint,"=?", nelec
   ! write(*,*) "integral of diagonal 1rdm=", radint2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!end of test 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     open(unit=3, file='refpoints.inp')
     read(3,*) npoints
     do i1=1,npoints      
     read(3,*) x,y,z, outputname
    !!!!!!Test 2: Compute exchange hole at reference point:!!!!!!!!!!!!!
    !!!!!!!!!!-With Becke model
    !!!!!!!!! -With 1rdm^2/density
      open(unit=4, file=outputname)
      write(*,*) "Reference point"
      write(*,*) x,y,z
      write(*,*) "------------------"
      call becke1(x,y,z,b,alf)
      Dens=Density(x,y,z)/2.d0
      write(*,*) "b=", b
      write(*,*) "a=", alf
      kk=0.0001d0 !interelectronic distance
      write(*,*) "Loop over distances"
      do i=1,1000
          kk=kk+0.005d0
          xchole=beckex(b,alf,kk)   !exchange hole model
          rdmv1=0.d0
          do k=1,nang !angular integral of the squared 1rdm
                px=x+kk*cosin(k)  !dcos(phi(k))*dsin(theta(k))
                py=y+kk*sinsin(k) !dsin(phi(k))*dsin(theta(k))
                pz=z+kk*cost(k)   !dcos(theta(k))
                rdmv1=rdmv1+wa(k)*(rDM1_alf(x,y,z,px,py,pz)*rDM1_alf(px,py,pz,x,y,z))
          end do
          write(4,*) kk, xchole, rdmv1/((Dens)*4.d0*pi) !write both models
      end do
      STOP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!end of test 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Test 3: Integrate xchole, becke model and 1rdm^2/dens!!!!!!!!!!!!!!!!
      write(*,*) "End of test 2"
      xcint=0.d0
      xcint2=0.d0
      do i=1,nrad
          rdmv1=0.d0
          rdmv2=0.d0
          do k=1,nang
             px=x+r(i)*dcos(phi(k))*dsin(theta(k))
             py=y+r(i)*dsin(phi(k))*dsin(theta(k))
             pz=z+r(i)*dcos(theta(k))
             rdmv1=rdmv1+wa(k)*(rDM1(x,y,z,px,py,pz)*rDM1(px,py,pz,x,y,z))
          end do
          xcint2=xcint2+wr(i)*r(i)**2.d0*(rdmv1)
          xchole=beckex(b,alf,r(i))
          xcint=xcint+wr(i)*xchole*(r(i)**2.d0)*4.d0*pi
      end do
      write(*,*) "----------Test 3-------------"
      write(*,*) "Integral of exchange hole="
      write(*,*) "Becke Model=", xcint
      if (uhf) then 
            write(*,*) "1rdm^2/dens=", xcint2/(dens)
      else        
            write(*,*) "1rdm^2/dens", xcint2/(dens/2.d0)
      end if        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    close(4)  
    !see the plot of the exchange hole vs the distance at (x,y,z) point.
    end do
    close(3)
    write(*,*) "-----Test xchole done-----------"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Compute C1 part of Coulomb Hole!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(unit=4, file="c1hole.out")
    kk=0.00001d0   !interelectronic distance
    do s=1,45 !for 20 interelectronic distances (EXAMPLE)
         hc1(s)=0.d0     
         do i=1,nrad       !loop for radial nodes (integrate over r reference dist)
              rho2s=0.d0   !spherically averaged pair density at r dist and s dist.
              do j=1,nang  !loop for angular nodes (to average r)
                   x=r(i)*dcos(phi(j))*dsin(theta(j))  !compute x,y,z points (\Vec{r})
                   y=r(i)*dsin(phi(j))*dsin(theta(j))
                   z=r(i)*dcos(theta(j))
                   rdmv1=0.d0 !spherically averaged density at s distance at r point. 
                   do k=1,nang
                        !compute interelectronic separation vector (px,py,pz) (to average s)
                        px=x+kk*dcos(phi(k))*dsin(theta(k)) !u are averaging coord 2, not distance s!!!
                        py=y+kk*dsin(phi(k))*dsin(theta(k))
                        pz=z+kk*dcos(theta(k))                                           
                        !spherical average of s
                        rdmv1=rdmv1+wa(k)*rDM1_alf(x,y,z,px,py,pz)*rDM1_alf(px,py,pz,x,y,z)
                   end do
                  ! write(*,*) "Point:", x,y,z
                   call becke1(x,y,z,b,alf)
                  ! write(*,*) "b=,alf=", b, alf
                   xchole=beckex(b,alf,kk)
                   write(*,*) "Xchole,edmv1=", xchole, rdmv1
                   rho2=(rdmv1/(4.d0*pi)-xchole*((Density(x,y,z)/2.d0)))*2.d0 !c1 part of Coulomb hole: hc(\Vec{r},s)
                   rho2s=rho2s+wa(j)*rho2   !spherical average of \Vec{r}
                   !write(*,*) "rho2s=", rho2s
              end do
              !write(*,*) "-----rho2s=", rho2s
              !write(*,*) "---hc1=", hc1(s)
              hc1(s)=hc1(s)+wr(i)*rho2s*r(i)**2.d0 !radial average of r
         end do
         write(4,*) kk, hc1(s)*kk**2.d0
         kk=kk+0.1d0
    end do 
end subroutine c1hole

subroutine becke1(x,y,z,b,alf) !give a point and returns a and b from eq 17 of Becke-Rousell
  use fractions
  use geninfo
  implicit none
  double precision, intent(in) :: x,y,z
  double precision, intent(out) :: b,alf
  !!!!!local variables!!!!!!!!!!!!!!!!!!
  double precision :: xx ! x in Becke-Rousell paper. 
  double precision :: dens, E_kin, Grad, Lapl, newtr
  double precision, parameter :: trsh=1d-16
 ! write(*,*) "In Becke 1"
  call dens_ops(x,y,z,Dens,E_kin,Grad,Lapl) !subroutine to compute densities (one spin)
 ! write(*,*) "Dens, E_kin, Grad, Lapl" !this quantities are for one spin
 ! write(*,*) Dens, E_kin, Grad, Lapl
     write(*,*) "x,y,z point"
     write(*,*) x,y,z
     xx=newtr(Dens,E_kin,Grad,Lapl)
     if (xx.gt.0.d0) then
        b=(xx**(3.d0)*dexp(-xx)*0.125d0*((pi*Dens)**(-1.d0)))**(hr)
        alf=xx*(b**(-1.d0))
     else
        alf=0.d0
        b=0.d0
     end if        
        write(*,*) "alf=", alf
end subroutine becke1        

function beckex(b,alf,s) !Becke-Roussel Exchange Hole, eqn 16
    use fractions
    use geninfo
    implicit none    
    double precision :: beckex
    double precision, intent(in) :: b,alf,s  !r and s are distances
    !!!local variables!!!!!!!!!!!!!!!!!!
    double precision :: p1, p2, p3
        if (b.eq.0.d0) then
                beckex=0.d0
        else        
        p1=0.0625d0*alf*(pi*b*s)**(-1.d0)
        p2=(alf*dabs(b-s)+1.d0)*dexp(-alf*dabs(b-s))
        p3=(alf*dabs(b+s)+1.d0)*dexp(-alf*dabs(b+s))
        beckex=p1*(p2-p3)
end if
end function beckex 

function newtr(Dens,E_kin,Grad,Lapl) !performs Newton-Raphson to find x
   use fractions  
   use geninfo   
   implicit none
   double precision :: newtr
   double precision, intent(in) :: Dens, E_kin, Grad, Lapl
   !!!local variables!!!!!
   double precision :: dif
   double precision :: x0,x
   double precision :: qval
   double precision :: f, fp, fpr, fr,k
   double precision :: q
   integer, parameter :: maxit=100
   double precision, parameter :: trsh=1d-11
   integer :: i,j
   double precision :: sf !scaling factor (step size)
   qval=q(Dens,E_kin,Grad,Lapl) !compute eq 20b
   if (qval.eq.0.d0) then
           qval=1d-15
   end if        
   write(*,*) "Q=", qval
   !if (dens.eq.0.d0) then
   !        k=0.d0
   !else        
           k=bih*pi**(bih)*Dens**(boh)*qval**(-1.d0) !right hand part of eq 21
   !end if        
   write(*,*) "k=", k
   x0=2.d0
   sf=1.d0
   if (dabs(k).lt.1d-16) then
         newtr=0.d0
   else      
     if (k.lt.0.d0) then
          do i=1,16
              sf=sf*0.1d0
              x=x0-sf
              fpr=fp(x)
              fr=f(x,k)
              if (fr.lt.0.d0) then
                      write(*,*) "x=", x
                      do j=1,maxit
                          fpr=fp(x)
                          fr=f(x,k)
                          newtr=x-fr*(fpr**(-1.d0))
                          dif=dabs(newtr-x)
                          if (dif.lt.trsh) then
                             write(*,*) "Converged after", j, "iterations"
                             write(*,*) "xx=", newtr
                             goto 100
                          end if
                          if (j.eq.maxit) then
                               write(*,*) "Newton-Raphson1 not converged"
                               write(*,*) "i=", x
                               write(*,*) "k=", k
                               STOP
                          end if
                          x=newtr
                          write(*,*) x
                      end do    
              end if  
          end do      
      else if (k.gt.0.d0) then
          do i=1,16
              sf=sf*0.1d0
              x=x0+sf
              fpr=fp(x)
              fr=f(x,k)
              if (fr.gt.0.d0) then
                      write(*,*) "x=", x
                     do j=1,maxit
                          fpr=fp(x)
                          fr=f(x,k)
                          newtr=x-fr*(fpr**(-1.d0))
                          dif=dabs(newtr-x)
                          if (dif.lt.trsh) then
                             write(*,*) "Converged after", j, "iterations"
                             write(*,*) "xx=", newtr
                             goto 100
                          end if
                          if (j.eq.maxit) then
                               write(*,*) "Newton-Raphson1 not converged"
                               write(*,*) "x=", x
                               write(*,*) "k=", k
                               STOP
                          end if
                          x=newtr
                          write(*,*) x
                     end do    
              end if
          end do  
        else
         write(*,*) "k is zero"
        newtr=0.d0   
       end if
     end if   
 100 CONTINUE   
end function

function f(x,k) !equation 21 of Becke-Roussel
    use fractions 
    use geninfo   
    implicit none  
    double precision, parameter :: trsh=1d-8
    double precision, intent(in) :: x,k 
    double precision :: f1,f 
    f1=x*dexp(-bih*x)*(x-2.d0)**(-1.d0)                             !nomes cal calcular q una vegada
    f=f1-k
end function

function fp(x) !derivative of equation 22
    use fractions
    use geninfo !for pi variable
    implicit none    
    double precision :: fp
    double precision, intent(in) :: x
    fp=-bih*dexp(-bih*x)*(x**2.d0-2.d0*x+3.d0)*(x-2.d0)**(-2.d0)

end function

function q(Dens,E_kin,Grad,Lapl) !equation 20b
    use fractions
    implicit none    
    double precision, intent(in) :: Dens,E_kin,Grad,Lapl
    double precision :: D
    double precision :: q
    double precision, parameter :: gmm=0.8d0
    q=bs*(Lapl-2.d0*gmm*D(Dens,E_kin,Grad))
end function

function D(Dens,E_kin,Grad) !equation 13b
    implicit none
    double precision,intent(in) :: Dens,E_kin,Grad
    double precision :: D
    double precision, parameter :: trsh=1d-15   
    if (dabs(Dens).lt.trsh) then
           D=E_kin
    else       
           D=E_kin-0.25d0*(Grad)*Dens**(-1.d0) !Grad-->scalar product of gradient
    end if       
  !  write(*,*) "D=", D
end function

subroutine dens_ops(x,y,z,Dens,E_kin,Grad,Lapl)
    !Compute Density, kinetic energy density, gradient, and laplacian
    !at x,y,z point (for what?????)
    use geninfo
    use wfxinfo
    implicit none
    !!!global variables!!!!!!!!!!!!!!!!!!
    double precision, intent(in) :: x,y,z
    double precision, intent(out) :: Dens,E_kin,Grad,Lapl
    !!!local variables!!!!!!!!!!!!!!!!!!!!!!
    !primitive derivatives (matrices)
    double precision, dimension(nprim) :: pval !value of primitive
    double precision, dimension(nprim) ::  dpxval, dpyval, dpzval !value of derivative of primitive
    double precision, dimension(nprim) :: fxval,fyval,fzval !value of fragment of derivative of prim.
    double precision, dimension(nprim) :: dpxx,dpyy,dpzz !value of double derivative of primitive
    !primitive derivatives (on the fly values)
    double precision :: dpx,dpy,dpz !fragment of derivative of primitives
    double precision :: dx,dy,dz !derivative of primitives
    double precision :: dxx,dyy,dzz !double derivative of primitive
    double precision :: xa,ya,za    !center of primitives
    double precision :: pr !primitive
    !functions
    double precision :: density, prim
    !double precision :: E_kn
    
    
    double precision :: gradx,grady,gradz
    double precision :: laplx,laply,laplz, laplx1, laplx2!, laply1, laply2, laplz1, laplz2
    !double precision :: Ek_x,Ek_y,Ek_z   !x,y,z components of kinetic energy density
    double precision :: Ekin_x,Ekin_y,Ekin_z
    double precision :: Tmult!, Tmult2
    integer :: i,j,k!,l,m,i2
    double precision :: grad2!, gx1,gx2,gy1,gy2,gz1,gz2
    
    
    Dens=Density(x,y,z) !compute 1spin density
    
    do i=1,nprim !compute and store value of primitives and its derivatives
        xa=cartes(Ra(i),1)
        ya=cartes(Ra(i),2)
        za=cartes(Ra(i),3)
        pr=Prim(x,y,z,i)
        pval(i)=pr !value of primitive i   
        call dprim(i,x,y,z,xa,ya,za,pr,dpx,dpy,dpz,dx,dy,dz) 
        fxval(i)=dpx
        fyval(i)=dpy
        fzval(i)=dpz
        dpxval(i)=dx
        dpyval(i)=dy
        dpzval(i)=dz
        call ddprim(i,x,y,z,xa,ya,za,pval(i),dx,dy,dz,dpx,dpy,dpz,dxx,dyy,dzz)
        dpxx(i)=dxx
        dpyy(i)=dyy
        dpzz(i)=dzz
        if (dyy.gt.100) then
           write(*,*) "dyy was too big"
           write(*,*) TMN(i,1), TMN(i,2), TMN(i,3)     
           write(*,*) "Primitive",i, pval(i), tmn(i,1), tmn(i,2), tmn(i,3)
           write(*,*) dyy
           write(*,*) "--------------------------------------------------"
        end if   
        if (dxx.gt.100) then
           write(*,*) "dyy was too big"
           write(*,*) TMN(i,1), TMN(i,2), TMN(i,3)
           write(*,*) "Primitive",i, pval(i), tmn(i,1), tmn(i,2), tmn(i,3)
           write(*,*) dxx
           write(*,*) "--------------------------------------------------"
        end if
    
    end do
    
    E_kin=0.d0
    Gradx=0.d0
    Grady=0.d0
    Gradz=0.d0
    Laplx=0.d0
    Laplx1=0.d0
    Laplx2=0.d0
    Laply=0.d0
    Laplz=0.d0
    Lapl=0.d0
    Ekin_x=0.d0
    Ekin_y=0.d0
    Ekin_z=0.d0
    Grad2=0.d0
    do i=1,noccmo
       if (Occ(i).gt.0.d0) then
        do j=1,nprim
          do k=1,nprim
             Tmult=T(i,j)*T(i,k)*Occ(i)
             !do i2=1,noccmo
             ! if (Occ(i).gt.0.d0) then
             !  do l=1,nprim
             !     do m=1,nprim
             !         Tmult2=T(i2,k)*T(i2,l)*Occ(i)
             !         gx1=dpxval(j)*pval(k)+pval(j)*dpxval(k)
             !        gx2=dpxval(l)*pval(m)+pval(l)*dpxval(m)
             !         gy1=dpyval(j)*pval(k)+pval(j)*dpyval(k)
             !         gy2=dpyval(l)*pval(m)+pval(l)*dpyval(m)
             !        gz1=dpzval(j)*pval(k)+pval(j)*dpzval(k)
             !         gz2=dpzval(l)*pval(m)+pval(l)*dpzval(m)
             !         Grad2=Grad2+Tmult*Tmult2*(gx1*gx2+gy1*gy2+gz1*gz2)
             !      end do
             !  end do
             ! end if
             !end do   
             Gradx=Gradx+Tmult*(dpxval(j)*pval(k)+pval(j)*dpxval(k))
             Grady=Grady+Tmult*(dpyval(j)*pval(k)+pval(j)*dpyval(k))
             Gradz=Gradz+Tmult*(dpzval(j)*pval(k)+pval(j)*dpzval(k))
             Laplx=Laplx+Tmult*((2.d0*dpxval(j)*dpxval(k))+pval(j)*dpxx(k) &
                        +dpxx(j)*pval(k))
             !Laplx1=Laplx1+Tmult*(2.d0*dpxval(j)*dpxval(k)) 
             !Laplx2=Laplx2+Tmult*(pval(j)*dpxx(k)+dpxx(j)*pval(k))
    
            ! Laply1=Laply1+Tmult*(2.d0*dpyval(j)*dpyval(k))
            ! Laply2=Laply2+Tmult*(pval(j)*dpyy(k)+dpyy(j)*pval(k))
    
            ! Laplz1=Laplz1+Tmult*(2.d0*dpzval(j)*dpzval(k))
            ! Laplz2=Laplz2+Tmult*(pval(j)*dpzz(k)+dpzz(j)*pval(k))
    
             Laply=Laply+Tmult*((2.d0*dpyval(j)*dpyval(k))+pval(j)*dpyy(k) &
                        +dpyy(j)*pval(k))
             Laplz=Laplz+Tmult*((2.d0*dpzval(j)*dpzval(k))+pval(j)*dpzz(k) &
                        +dpzz(j)*pval(k))
            ! write(*,*) "i=",i, "j=",j, "k=",k
           !  write(*,*) "pval(j)*dpyy(k)", pval(j), dpyy(k)
            ! write(*,*) "pval(k)*dpyy(j)", pval(k), dpyy(j)
             if (dpyy(k).gt.100) then
                     write(*,*) "dpyy very big"
                     write(*,*) dpyy(k), k, Tmult
                     write(*,*) "INFO ABOUT PRIMITIVE"
                     write(*,*) "PVAL=", pval(k)
                     write(*,*) "Coeff", TMN(k,1), TMN(k,2), TMN(k,3)
             end if
             if (dpzz(k).gt.100) then
                     write(*,*) "dpzz very big"
                     write(*,*) dpzz(k), k, Tmult
                     write(*,*) "Pval=", pval(k)
                     write(*,*) "Coeff", TMN(k,1), TMN(k,2), TMN(k,3)
             end if
           ! Lapl=Lapl+Tmult*(2.d0*(dpxval(j)*dpxval(k)+dpyval(j)*dpyval(k)+dpzval(j)*dpzval(k))+ &
           !                 pval(j)*(dpxx(k)+dpyy(k)+dpzz(k))+pval(k)*(dpxx(j)+dpyy(j)+dpzz(j)))
             Ekin_x=Ekin_x+Tmult*dpxval(j)*dpxval(k)
             Ekin_y=Ekin_y+Tmult*dpyval(j)*dpyval(k)
             Ekin_z=Ekin_z+Tmult*dpzval(j)*dpzval(k)
         end do
       end do
      end if 
    end do 
    write(*,*) "DENSITY=", Dens
    write(*,*) "GRADIENT="
    write(*,*) Gradx, grady, gradz
    Grad=(gradx/2.d0)**(2.d0)+(grady/2.d0)**(2.d0)+(gradz/2.d0)**(2.d0) !scalar product of gradient (one spin)
    write(*,*) Grad
    write(*,*) "Gradient product"
    write(*,*) Grad2
    write(*,*) "LAPLACIAN="
    write(*,*) Laplx, Laply, Laplz
    !write(*,*) Laplx1, Laplx2
    !write(*,*) Laply1, Laply2
    !write(*,*) Laplz1, Laplz2
    Lapl=Laplx+Laply+Laplz
    write(*,*) Lapl
    write(*,*) Laplx+Laply+Laplz
    write(*,*) "Kinetic energy"
    write(*,*) Ekin_x, Ekin_y, Ekin_z
    write(*,*) Ekin_x/2.d0,Ekin_y/2.d0, Ekin_z/2.d0
    E_kin=(Ekin_x+Ekin_y+Ekin_z)
    write(*,*) E_kin
    
    !compute all one spin
    if (clsh) then
            write(*,*) "Closed shell, Alpha=Beta"
            Dens=Dens/2.d0
           !Grad=Grad/2.d0
            Lapl=Lapl/2.d0
            E_kin=E_kin/2.d0
    else
            write(*,*) "Caution: Open shell"
    end if        
    end subroutine
    
    subroutine dprim(mu,x,y,z,xa,ya,za,pr,dpx,dpy,dpz,dx,dy,dz) 
    use geninfo
    implicit none
    double precision, intent(in) :: x,y,z,xa,ya,za,pr
    double precision, intent(out) :: dpx, dpy, dpz,dx,dy,dz
    integer, intent(in) :: mu
    double precision :: fx1,fy1,fz1
    double precision :: expr
    dpx=0.d0
    dpy=0.d0
    dpz=0.d0
    if (dabs(x-xa).gt.0d-5) then
            fx1=dble(tmn(mu,1))*(x-xa)**(-1.d0) !avoid div 0
    else
            fx1=0.d0
            if (tmn(mu,1).eq.1) fx1=1.d0         !see the development in the notebook
    end if        
    if (dabs(y-ya).gt.0d-5) then
            fy1=dble(tmn(mu,2))*(y-ya)**(-1.d0)
    else
            fy1=0.d0
            if (tmn(mu,2).eq.1) fy1=1.d0
    end if
    if (dabs(z-za).gt.0d-5) then
            fz1=dble(tmn(mu,3))*(z-za)**(-1.d0)
    else
            fz1=0.d0
            if (tmn(mu,3).eq.1) fz1=1.d0
    end if  
     dpx=fx1-2.d0*Alpha(mu)*(x-xa)
     dpy=fy1-2.d0*Alpha(mu)*(y-ya)
     dpz=fz1-2.d0*Alpha(mu)*(z-za)
    
    !compute derivative of the primitive 
    if (((dabs(x-xa).eq.0.d0)).and.(tmn(mu,1).eq.1)) then !x=xA and n=1
            dx=dpx*(y-ya)**(dble(TMN(mu,2)))*(z-za)**(dble(TMN(mu,3)))*expr(x,y,z,xa,ya,za,mu)
    else
            dx=dpx*pr
    end if
    if (((dabs(y-ya).eq.0.d0)).and.(tmn(mu,2).eq.1)) then !y=yA and m=1
            dy=dpy*(x-xa)**(dble(TMN(mu,1)))*(z-za)**(dble(TMN(mu,3)))*expr(x,y,z,xa,ya,za,mu)
    else
            dy=dpy*pr  
    end if  
    if (((dabs(z-za).eq.0.d0)).and.(tmn(mu,3).eq.1)) then !z=zA and l=1
            dz=dpz*(x-xa)**(dble(TMN(mu,1)))*(y-ya)**(dble(TMN(mu,2)))*expr(x,y,z,xa,ya,za,mu)
    else
            dz=dpz*pr
    end if
    
    end subroutine dprim
    
    subroutine ddprim(mu,x,y,z,xa,ya,za,pval,dx,dy,dz,fxval,fyval,fzval,dxx,dyy,dzz)
    use geninfo  
    implicit none
    double precision, intent(in) :: x,y,z,xa,ya,za
    double precision, intent(in) :: pval !value of the primitive
    double precision, intent(in) :: dx,dy,dz !value of x y and z derivatives of prim
    double precision, intent(in) :: fxval,fyval,fzval !part of the derivative of prim
    integer, intent(in) :: mu
    double precision, intent(out) :: dxx,dyy,dzz
    double precision :: px, py, pz
    double precision :: expr
    if (dabs(x-xa).gt.0.d0) then
            px=-pval*TMN(mu,1)*(x-xa)**(-2.d0) !avoid div 0
            dxx=dx*fxval+pval*(-2.d0*Alpha(mu))+px
    else
            if (TMN(mu,1).eq.1) then
               !  dxx=-2.d0*Alpha(mu)*dx
                  dxx=0.d0
            else if (TMN(mu,1).eq.2) then
                    dxx=2.d0*(y-ya)**(dble(TMN(mu,2)))*(z-za)**(dble(TMN(mu,3)))*expr(x,y,z,xa,ya,za,mu)
            else if (TMN(mu,1).eq.0) then 
                    dxx=pval*(-2.d0*Alpha(mu))  
               
            else
                    dxx=0.d0        
            end if
    end if
    if (dabs(y-ya).gt.0.d0) then
            py=-pval*TMN(mu,2)*(y-ya)**(-2.d0)
            dyy=dy*fyval+pval*(-2.d0*Alpha(mu))+py
    else
            if (TMN(mu,2).eq.1) then
                 !   dyy=-2.d0*Alpha(mu)*dy
                    dyy=0.d0
            else if (TMN(mu,2).eq.2) then
                    dyy=2.d0*(x-xa)**(dble(TMN(mu,1)))*(z-za)**(dble(TMN(mu,3)))*expr(x,y,z,xa,ya,za,mu)
                   
     
            else if (TMN(mu,2).eq.0) then
                    dyy=pval*(-2.d0*Alpha(mu))     
                
            else
                   ! dyy=pval*(-2.d0*Alpha(mu))     
                    dyy=0.d0
            end if
    end if
    if (dabs(z-za).gt.0.d0) then
            pz=-pval*TMN(mu,3)*(z-za)**(-2.d0)
            dzz=dz*fzval+pval*(-2.d0*Alpha(mu))+pz
    else
            if (TMN(mu,3).eq.1) then
                   write(*,*) "dzz not zero", TMN(mu,1),TMN(mu,2),TMN(mu,3)
                  ! dzz=-2.d0*Alpha(mu)*dz
                  ! write(*,*) dzz,pval
                   dzz=0.d0
            else if (TMN(mu,3).eq.2) then
                   dzz=2.d0*(x-xa)**(dble(TMN(mu,1)))*(y-ya)**(dble(TMN(mu,2)))*expr(x,y,z,xa,ya,za,mu)      
            else if (TMN(mu,3).eq.0) then
                   dzz=pval*(-2.d0*Alpha(mu))
                  
            else
                  ! dzz=pval*(-2.d0*Alpha(mu))
                   dzz=0.d0
            end if
    end if        
    end subroutine
    
