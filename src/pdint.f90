subroutine pdint()
!right now we are not using this function for anything...    
!computes SD pair density and its integral
use geninfo !contains cartes and natoms
use intrainfo !contains nquad and cent
implicit none
integer :: i,j,k,l,nqmax,sm, sm2, ir,ngrid,nrd,nan
double precision, allocatable, dimension(:) :: x0, Alp, N
double precision, parameter :: trsh=1.d-1, zero=1.d-15, trs=5.d-1
double precision :: Density, rDM1, Gauss, Dens_a, Dens_b
double precision :: Val_r,Dens_1,Dens_2, r
double precision :: Pk,bn
double precision :: x1,x2,y1,y2,z1,z2
double precision :: px1,py1,pz1,px2,py2,pz2
!variabls for integrals
double precision, allocatable, dimension(:) :: xl_i, wl_i, wlb, radius, wtot
double precision, allocatable, dimension(:,:) :: r_lb, rr
double precision, allocatable, dimension(:) :: phi, theta !angles
double precision :: PD, PD2, xs
double precision :: a, b !limits of the integral
 nqmax=(natoms*(natoms-1))/2
 allocate(x0(nqmax))
 allocate(Alp(nqmax))
 allocate(N(nqmax))
 write(*,*) "In subroutine pdint"  
 write(*,*) "natoms=", natoms, nqmax  
 sm=0
 write(*,*) "start loop"
 nrd=5
 nan=110
 ngrid=nrd*nan
 allocate(rr(3,ngrid))
 allocate(xl_i(nrd))
 allocate(wl_i(nrd))
 allocate(radius(nrd))
 allocate(r_lb(3,nan))
 allocate(wlb(nan))
 allocate(phi(nan))
 allocate(theta(nan))
 a=0.d0
 b=0.00006
 !compute Gauss-Legendre nodes and weights
 call sub_GauLeg(-1.d0,1.d0,xl_i,wl_i,nrd)
 !compute radial nodes from a to b
 do ir=1,nrd
     radius(ir)=(b-a)*0.5d0*xl_i(ir)+(a+b)*0.5d0
     wl_i(ir)=(b-a)*0.5d0*wl_i(ir)   !radial integration (see pdf) 
     write(*,*) radius(ir), wl_i(ir)
 end do
 !compute Gauss-Levedev nodes and weights
 call LD0110(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAn)
  do i=1,nAn   !compute angles     
       theta(i)= dacos(r_lb(3,i))
       if (dsin(theta(i)).ne.0.d0) then
            xs=r_lb(1,i)/dsin(theta(i))
            if (xs.gt.1.d0) xs=1.d0
            if (xs.lt.-1.d0) xs=-1.d0
            phi(i)=dacos(xs)
            if (r_lb(2,i).lt.0.d0) phi(i)=-phi(i)
       else
            phi(i)=0.d0
       end if
       Wlb(i)=Wlb(i)*4.d0*pi !store 4pi factor on the weight (see solid angle integral)
  end do
 open(unit=3, file="check")
 allocate(wtot(ngrid)) 
 !compute total grid points
 sm=0
 do i=1,nrd
   do j=1,nan
       sm=sm+1
       rr(1,sm)=radius(i)*dcos(phi(j))*dsin(theta(j))
       rr(2,sm)=radius(i)*dsin(phi(j))*dsin(theta(j))
       rr(3,sm)=radius(i)*dcos(theta(j))
       wtot(sm)=wl_i(i)*wlb(j)
   end do
 end do  
 !loop for atomic pairs
 write(*,*) "Loop for atomic pairs"
 sm=0
 do i=1,natoms-1
    do j=i+1,natoms
         !center of atoms
         x1=cartes(i,1)
         y1=cartes(i,2)
         z1=cartes(i,3)
         x2=cartes(j,1)
         y2=cartes(j,2)
         z2=cartes(j,3)
         sm=sm+1
         !reference distance
         x0(sm)=dsqrt((x1-x2)**2.d0+(y1-y2)**2.d0+(z1-z2)**2.d0)
         !compute pair density integral
         sm2=0  
         Pk=0      
         do k=1,nrd
            PD2=0.d0
            do l=1,nan
            !compute pair density
                sm2=sm2+1
                px1=rr(1,sm2)+x1
                px2=rr(1,sm2)+x2
                py1=rr(2,sm2)+y1
                py2=rr(2,sm2)+y2
                pz1=rr(3,sm2)+z1
                pz2=rr(3,sm2)+z2
                if (uhf) then
                    Dens_1=Dens_a(px1,py1,pz1)+Dens_b(px1,py1,pz1)
                    Dens_2=Dens_a(px2,py2,pz2)+Dens_b(px2,py2,pz2)
                else        
                    Dens_1=Density(px1,py1,pz1)
                    Dens_2=Density(px2,py2,pz2)
                end if
                PD=Dens_1*Dens_2-(rDM1(px1,py1,pz1,px2,py2,pz2)**2.d0)
                PD2=PD2+wtot(l)*PD
            end do 
            Pk=Pk+PD2
         end do   
         !Pk=Dens_A*Dens_B-(rDM1(x1,y1,z1,x2,y2,z2)**2.d0)
         
         write(3,*) "Atoms", i, "and", j
         write(3,*) "x0=", x0(sm)
         write(3,*) "Peak=", Pk
         Pk=Pk*(x0(sm)**2.d0)
         bn=an(i)*an(j)
         Alp(sm)=pi*0.0625d0*((Pk*dexp(1.d0))**2.d0)*(bn**(-2.d0))
         N(sm)=dsqrt(pi)*0.25d0*Alp(sm)**(-1.5d0) !normalization factor if integral is 1
         N(sm)=bn/N(sm)                           !full factor for exponential
         write(3,*) "Alpha=", Alp(sm)
         write(3,*) "Integral value=", bn
    end do
 end do
 close(3)
 !
 open(unit=3,file='radial')
 r=0.d0
 Val_r=0.d0
 do i=1,500
     r=r+0.03d0
     do j=1,nqmax
        !write(*,*) "x0,Alp,Nj"
        !write(*,*) x0(j), Alp(j), N(j)
        Val_r=Val_r+Gauss(r,x0(j),Alp(j),N(j))
     end do
     write(3,*) r, Val_r
     Val_r=0.d0
 end do      
 
 close(3)
 deallocate(rr)
 deallocate(xl_i)
 deallocate(wl_i)
 deallocate(radius)
 deallocate(r_lb)
 deallocate(wlb)
 deallocate(phi)
 deallocate(theta)
 deallocate(x0)
 deallocate(Alp)
 deallocate(N)                 
end subroutine pdint               

function Gauss(r,x0,Alp,N)
  double precision :: Gauss
  double precision :: r,x0,Alp,N
  double precision :: minim, a1  
    a1=(dsqrt(Alp))**(-1.d0)
    minim=x0-a1
    if (r.ge.minim) then
        Gauss=N*((r-x0+a1)**2.d0)*dexp(-Alp*(r-x0+a1)**2.d0)
    else
        Gauss=0.d0 !neglect the peak from the left
    end if
end function


