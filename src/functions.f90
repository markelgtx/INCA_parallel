!*******************************************************************************
 function Prim(x,y,z,npr) 
   use geninfo
   implicit none
   double precision :: Prim
   double precision, intent(in) :: x, y, z
   integer, intent(in) :: npr
   
   Prim=(x-(cartes(Ra(npr),1)))**(dble(TMN(npr,1)))* &
        (y-(cartes(Ra(npr),2)))**(dble(TMN(npr,2)))* &
        (z-(cartes(Ra(npr),3)))**(dble(TMN(npr,3)))* &
        dexp(-Alpha(npr)*((x-cartes(Ra(npr),1))**2.d0+(y-(cartes(Ra(npr),2)))**2.d0+ & 
                         (z-(cartes(Ra(npr),3)))**2.d0))
              
 end function

  function expr(x,y,z,xa,ya,za,npr) 
   use geninfo
   implicit none
   double precision :: expr
   double precision, intent(in) :: x, y, z,xa,ya,za
   integer, intent(in) :: npr
   expr=dexp(-Alpha(npr)*((x-xa)**2.d0+(y-ya)**2.d0+(z-za)**2.d0))
   write(*,*) "Expr=-------------", expr
 end function

   
!**************************************************************************
function sd(x,y,z,npr) !second derivative of primitives analitically
use geninfo
implicit none
double precision :: sd, Prim
double precision, intent(in) :: x, y, z
integer, intent(in) :: npr
double precision :: laplx, laply, laplz
double precision :: t,m,n
double precision :: xa,ya,za

t=0.d0 
m=0.d0 
n=0.d0
xa=0.d0
ya=0.d0 
za=0.d0
laplx=0.d0 
laply=0.d0
laplz=0.d0
sd=0.d0

t=dble(TMN(npr,1))
m=dble(TMN(npr,2))
n=dble(TMN(npr,3))
xa=cartes(Ra(npr),1)
ya=cartes(Ra(npr),2)
za=cartes(Ra(npr),3)

laplx=(Prim(x,y,z,npr)/(x**t))*(t*(t-1.0d0)*x**(t-2.0d0) &
+4.0d0*t*(x**(t-1.0d0))*Alpha(npr)*(xa-x)+((x**t)*Alpha(npr) &
*(4.0d0*Alpha(npr)*(xa-x)-2.0d0)))

laply=(Prim(x,y,z,npr)/(y**m))*(m*(m-1.0d0)*y**(m-2.0d0) &
+4.0d0*m*(y**(m-1.0d0))*Alpha(npr)*(ya-y)+((y**m)*Alpha(npr) &
*(4.0d0*Alpha(npr)*(ya-y)-2.0d0)))

laplz=(Prim(x,y,z,npr)/(z**n))*(n*(n-1.0d0)*z**(n-2.0d0) &
+4.0d0*n*(z**(n-1.0d0))*Alpha(npr)*(za-z)+((z**n)*Alpha(npr) &
*(4.0d0*Alpha(npr)*(za-z)-2.0d0)))

sd=laplx+laply+laplz

end function
    
      
!**************************************************************************
   function AO(x,y,z,cao)
   use geninfo
   use loginfo      !cao from 1 to nao, 1->1S, 2->2S, 3->2Px, 4->2Py, ...
   implicit none
   double precision :: AO
   double precision :: Prim
   double precision, intent(in) :: x, y, z
   integer, intent(in) :: cao !choosen ao in the input file
   integer :: smm !sum for all the primitives
   integer :: i
   integer :: sm, kk, k, kp
   smm=0
   AO=0.d0
   sm=0
   k=1
   if (cao.eq.1) then 
     k=1
   else if (cao.eq.2) then
     kp=1
     k=k+npao(1)
   else  
     kp=cao-1    !kp=previous AO 
     do i=1,kp  !set k starting primitive, Flg(k)
      k=k+npao(i)  !sum the primitives until the current AO
     end do
   end if
     
   kk=k+npao(cao)-1 !kk=last primitive of the current AO (cao)
    do i=k,kk        
      AO= AO + Flg(i) * Prim(x,y,z,i)    
    end do
 
  end function AO
   
!*************************************************************************************   
   function MoOr(x,y,z,mo)
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MoOr
   double precision :: Prim
   integer, intent(in) :: mo
   integer :: i    
   MoOr=0.d0
   do i=1,nprim
      MoOr= MoOr+ T(mo,i)*Prim(x,y,z,i)
   end do   
   end function
 !**************************************************************************
  function MO_a(x,y,z,mo) 
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x,y,z
   integer, intent(in) :: mo
   double precision :: MO_a
   double precision :: Prim
   integer :: i
   MO_a=0.d0
   do i=1,nprim
     MO_a= MO_a + T_a(mo,i) * Prim(x,y,z,i)
   end do
   end function
   
   function MO_b(x,y,z,mo)
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x,y,z
   integer, intent(in) :: mo
   double precision :: MO_b
   double precision :: Prim
   integer :: i
   MO_b=0.d0
   do i=1,nprim
     MO_b=MO_b+ T_b(mo,i)*Prim(x,y,z,i)
   end do
   end function   
    
 !**********************************************************************
   function Density(x,y,z)   !from wfx file
   use wfxinfo
   use geninfo 
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MoOr, Mo_a, Mo_b
   double precision :: Density, dens_a, dens_b
   integer :: i
    Density=0.d0
    if (corr) then
      do i=1,noccmo
        Density= Density + Occ(i) * MoOr(x,y,z,i)**2.d0  
      end do        
    !else if (uhf) then
    !     if (dble(nalfae).gt.0.1d0) then
    !         do i=1,nalfae
    !            dens_a=dens_a+MO_a(x,y,z,i)**2.d0
    !         end do
    !     end if
    !       if (dble(nbetae).gt.0.1d0) then
    !          do i=1,nbetae
    !            dens_b=dens_b+MO_b(x,y,z,i)**2.d0
    !          end do
    !       end if            
    !       density=dens_a+dens_b
    else        
       do i=1,noccmo 
           Density=Density+ Occ(i)*MoOr(x,y,z,i)**2.d0
       end do   
    end if
    end function
!*******************************************************************     
   function Dens_a(x,y,z) 
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MO_a, MoOr
   double precision :: Dens_a
   integer :: i
   if ((uhf).or.(opsh)) then
     Dens_a=0.d0
     do i=1,nalfae
       Dens_a=Dens_a+MO_a(x,y,z,i)**2.0d0
     end do
   else if ((clsh).and.(rhf)) then
      do i=1,noccmo
         Dens_a=Dens_a+Occ(i)*MoOr(x,y,z,i)**2.d0
      end do  
      Dens_a=Dens_a/2.d0
   end if   
    end function  
 !****************************************************************
   function Dens_b(x,y,z) 
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MO_b, MoOr
   double precision :: Dens_b
   integer :: i
   if ((uhf).or.(opsh)) then
     Dens_b=0.d0
     do i=1,nalfae
       Dens_b=Dens_b+MO_b(x,y,z,i)**2.0d0
     end do
    else if ((clsh).and.(rhf)) then
          do i=1,noccmo
         Dens_b=Dens_b+Occ(i)*MoOr(x,y,z,i)**2.d0
      end do
      Dens_b=Dens_b/2.d0
     end if  
    end function  
!*****************************************************************
    
function lapl(x,y,z) !computes Laplacian numerically
use wfxinfo
implicit none 
double precision :: lapl                           !this could be a function!!
double precision, intent(in) :: x,y,z
double precision :: spx, spy, spz !laplacian, square partial x, s.p.y, spz
double precision :: density,h

h=0.0000001d0

spx=(Density(x+h,y,z)+Density(x-h,y,z)-2.d0*Density(x,y,z))/(h**2.0d0)
spy=(Density(x,y+h,z)+Density(x,y-h,z)-2.d0*Density(x,y,z))/(h**2.0d0)
spz=(Density(x,y,z+h)+Density(x,y,z-h)-2.d0*Density(x,y,z))/(h**2.0d0)

lapl=spx+spy+spz

end function

!************************************************************************* 

!function Laplacian(x,y,z)  !laplacian analitically using fd and sd functions (first derivative and second derivative of primitives  !no funtziona
!use wfxinfo
!use geninfo
!implicit none
!double precision :: Laplacian, Laplx, Laply, Laplz, dp, d2p
!double precision :: prim
!double precision, intent(in) :: x, y, z  
!integer :: i,j,k

!Laplacian=0.d0                      !consider the UHF case!
!Laplx=0.d0
!Laply=0.d0
!Laplz=0.d0

!do i=1,noccmo    
!          !current molecular orbital value (also valid for UHF)
!   do j=1,nprim
!    do k=1,nprim    
!     Laplx=Laplx+2*(T(i,j)*T(i,k)*d2p(x,y,z,j,1)*Prim(x,y,z,k)+(T(i,j)*dp(x,y,z,j,1))**2)   
!     Laply=Laply+2*(T(i,j)*T(i,k)*d2p(x,y,z,j,2)*Prim(x,y,z,k)+(T(i,j)*dp(x,y,z,j,2))**2)   
!     Laplz=Laplz+2*(T(i,j)*T(i,k)*d2p(x,y,z,j,3)*Prim(x,y,z,k)+(T(i,j)*dp(x,y,z,j,3))**2)            
!     end do                
!    end do   
!end do
!Laplacian=laplx+laply+laplz
!end function

function dp(x,y,z,j,part)   !no funtziona
use geninfo
implicit none 
double precision, intent(in) :: x,y,z
double precision :: dp, prim
integer, intent(in) :: part !1->x, 2->2, 3->z
integer, intent(in) :: j !number of primitive
double precision :: ti, m, n, xa, ya, za

  ti=dble(TMN(j,1))
  m=dble(TMN(j,2))
  n=dble(TMN(j,3))
  xa=cartes(Ra(j),1)
  ya=cartes(Ra(j),2)
  za=cartes(Ra(j),3)

 if (part.eq.1) then !first partial derivative of prim with respect to x
   dp=(Prim(x,y,z,j)*(ti*x**(ti-1.0d0))/(x**ti))+Prim(x,y,z,j)*2.0d0*Alpha(j)*(xa-x)
 else if (part.eq.2) then !y
   dp=(Prim(x,y,z,j)*(m*y**(m-1.0d0))/(y**m))+Prim(x,y,z,j)*2.0d0*Alpha(j)*(ya-y)
 else if (part.eq.3) then !z
   dp=(Prim(x,y,z,j)*(n*z**(n-1.0d0))/(z**n))+Prim(x,y,z,j)*2.0d0*Alpha(j)*(za-z)
 end if   
 
end function

function d2p(x,y,z,j,part)   !no funtziona
use geninfo
implicit none 
double precision, intent(in) :: x,y,z
double precision :: d2p, prim
integer, intent(in) :: part !1->x, 2->2, 3->z
integer, intent(in) :: j !number of primitive
double precision :: ti, m, n, xa, ya, za
  ti=dble(TMN(j,1))
  m=dble(TMN(j,2))
  n=dble(TMN(j,3))
  xa=cartes(Ra(j),1)
  ya=cartes(Ra(j),2)
  za=cartes(Ra(j),3)

 if (part.eq.1) then !first partial derivative of prim with respect to x
   d2p=(Prim(x,y,z,j)/(x**ti))*(ti*(ti-1.0d0)*x**(ti-2.0d0) &
    +4.0d0*ti*(x**(ti-1.0d0))*Alpha(j)*(xa-x)+((x**ti)*Alpha(j) &
    *(4.0d0*Alpha(j)*(xa-x)-2.0d0)))
 else if (part.eq.2) then !y
   d2p=(Prim(x,y,z,j)/(y**m))*(m*(m-1.0d0)*y**(m-2.0d0) &
    +4.0d0*m*(y**(m-1.0d0))*Alpha(j)*(ya-y)+((y**m)*Alpha(j) &
    *(4.0d0*Alpha(j)*(ya-y)-2.0d0)))
 else if (part.eq.3) then !z
   d2p=(Prim(x,y,z,j)/(z**n))*(n*(n-1.0d0)*z**(n-2.0d0) &
    +4.0d0*n*(z**(n-1.0d0))*Alpha(j)*(za-z)+((z**n)*Alpha(j) &
    *(4.0d0*Alpha(j)*(za-z)-2.0d0)))
 end if   

end function
              
!*******************************************************************

subroutine gradient(x,y,z,gradx,grady,gradz) !computes the gradient numerically
                                             !now we won't use that but could be useful
 implicit none
 double precision, intent(in) :: x,y,z
 double precision, intent(out) :: gradx, grady, gradz
 double precision :: density, h

 h=0.0000000000001d0 !change in the derivative

 gradx=(Density(x+h,y,z)-Density(x,y,z))/h
 grady=(Density(x,y+h,z)-Density(x,y,z))/h
 gradz=(Density(x,y,z+h)-Density(x,y,z))/h
  
end subroutine
    
!*****************************************************

  function rDM1(x1,y1,z1,x2,y2,z2) !1-electron reduced density matrix (HF case)
                                   !for alpha and beta electrons
   use wfxinfo                     
   use geninfo
   implicit none
   double precision :: rDM1
   double precision :: rDM1_a, rDM1_b
   double precision :: MoOr, MO_a, MO_b
   double precision, intent(in) :: x1, y1, z1, x2, y2, z2
   integer :: i  
       if (uhf) then
         if (nelec.eq.1) then   !only one electron 
           rDM1=0.d0
           do i=1,noccmo
                rDM1=rDM1+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
           end do
         else        
           rDM1_a=0.d0
           rDM1_b=0.d0  
           if (dble(nalfae).gt.0.1d0) then  
             do i=1,nalfae
                rDM1_a=rDM1_a+MO_a(x1,y1,z1,i)*MO_a(x2,y2,z2,i)
             end do
           end if
           if (dble(nbetae).gt.0.1d0) then  
              do i=1,nbetae
                rDM1_b=rDM1_b+MO_b(x1,y1,z1,i)*MO_b(x2,y2,z2,i)
              end do
           end if  
           rDM1=rDM1_a+rDM1_b
         end if  
      else if (clsh) then !restricted Hartree Fock
           do i=1,noccmo
             if (occ(i).gt.0.d0) then
                rDM1=rDM1+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
             end if   
           end do
       end if    
   end function

   function rDM1_alf(x1,y1,z1,x2,y2,z2)
   use wfxinfo
   use geninfo
   implicit none
   double precision :: rdm1_alf
   double precision :: MoOr, MO_a, MO_b
   double precision, intent(in) :: x1, y1, z1, x2, y2, z2
   integer :: i
    rDM1_alf=0.d0
   if (uhf) then 
        if (nelec.eq.1) then   !only one electron
           rDM1_alf=0.d0
           do i=1,noccmo
                rDM1_alf=rDM1_alf+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
           end do
        else
           rDM1_alf=0.d0
           if (dble(nalfae).gt.0.1d0) then
             do i=1,nalfae
                rDM1_alf=rDM1_alf+MO_a(x1,y1,z1,i)*MO_a(x2,y2,z2,i)
             end do
           else 
                write(*,*) "No alpha electrons, rdm1alf=0"   
                rDM1_alf=0.d0
           end if
        end if  
    else if (clsh) then
        rdm1_alf=0.d0   
        do i=1,noccmo
              if (occ(i).gt.0.d0) then
                 rDM1_alf=rDM1_alf+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
              end if
        end do
        rDM1_alf=rDM1_alf/2.d0
    end if 
    end function   
        
   function rDM1_bet(x1,y1,z1,x2,y2,z2)
   use wfxinfo
   use geninfo
   implicit none
   double precision :: rdm1_bet
   double precision :: MoOr, MO_a, MO_b
   double precision, intent(in) :: x1, y1, z1, x2, y2, z2
   integer :: i
   if (uhf) then
        if (nelec.eq.1) then   !only one electron
           rDM1_bet=0.d0
           do i=1,noccmo
                rDM1_bet=rDM1_bet+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
           end do
        else
           rDM1_bet=0.d0
           if (dble(nalfae).gt.0.1d0) then
             do i=1,nalfae
                rDM1_bet=rDM1_bet+MO_a(x1,y1,z1,i)*MO_a(x2,y2,z2,i)
             end do
           else
                write(*,*) "No alpha electrons, rdm1alf=0"
                rDM1_bet=0.d0
           end if
        end if
    else if (clsh) then
        rdm1_bet=0.d0
        do i=1,noccmo
              if (occ(i).gt.0.d0) then
                 rDM1_bet=rDM1_bet+MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
              end if
        end do
    end if
    end function



