!*********************************************************************
!!Creates a cubefile (.cube extension) which can be opened on vmd    !
!!and see the 3D representation of the magnitude we want i.e.        !
!!density, gradient of the density, laplacian, MOs, AOs, etc         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cubefile(a,cubename) 
!creates a cubefile with a magnitude we choose f(a)
use geninfo
use cubeinfo
implicit none
character*40 :: cubename
character*5 :: cubepart
integer, intent(in) :: a !a=1 primitive, a=2 AO, a=3 MO, a=5 Density, a=6 Laplacian
integer :: i, j, aa 
double precision, parameter :: zero=0.d0

aa=0

if ((uhf).and.((a.eq.3).or.(a.eq.5))) then    !if we have an unrestricted wavefunction
       do i=1,3                              !do a loop for alpha, beta and total wf.
                aa=aa+1                       !3 cubefiles will be generated
                if (i.eq.1) cubepart="alpha"  !insert description in the name of the .cube
                if (i.eq.2) cubepart="beta_"  !so that we know what we have represented
                if (i.eq.3) cubepart="all__"
                open(unit=2,file=cubepart//cubename) 
                write(2,*) "CUBE FILE"
                write(2,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
                write(2,*) natoms, center(1), center(2), center(3)
                write(2,*) np(1), step(1), zero, zero
                write(2,*) np(2), zero, step(2), zero
                write(2,*) np(3), zero, zero, step(3)
                do j=1,natoms
                     write(2,*) an(j), chrg(j), cartes(j,1), cartes(j,2), cartes(j,3)
                end do
                call fwrite(a,aa)  !writes the values of the function in the cubepoints
                close(2)
       end do
else  
       open(unit=2,file=cubename)
       write(2,*) "CUBE FILE"
       write(2,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
       write(2,*) natoms, center(1), center(2), center(3)
       write(2,*) np(1), step(1), zero, zero
       write(2,*) np(2), zero, step(2), zero
       write(2,*) np(3), zero, zero, step(3)
       do i=1,natoms
            write(2,*) an(i), chrg(i), cartes(i,1), cartes(i,2), cartes(i,3)
        end do
        call fwrite(a,aa) !writes the value of a function in each point
        close(2)
end if
end subroutine cubefile


subroutine fwrite(a,aa)    !evaluates a given function in the specified points.
use geninfo
use cubeinfo
implicit none
integer, intent(in) :: a, aa
integer :: i,j,k
double precision :: Prim, AO, MoOr, MO_a, MO_b, Density, Dens_a, Dens_b, Lapl !Laplacian
double precision :: Dens, MO_t
double precision, dimension(3) :: xm !minimum point when np odd or even
double precision :: x,y,z
logical :: num

num=.true. !interruptor manual per numerical or analitical !change this to choose in the input
x=0.0d0
y=0.0d0
z=0.0d0
xm=0.d0
do i=1,3
    if (mod(np(i),2).eq.0)then !odd number
         xm(i)=center(i)-(step(i)/2.0d0)-((dble(np(i))-2.0d0)/2.0d0)*step(i)
    else !even number  !Calculates the starting point acording to cubeinfo data
         xm(i)=center(i)-((dble(np(i))-1.0d0)/2.0d0)*step(i)
    end if  
end do

do i=1,np(1)         !Depending on a calculates the point with a different function
     x=xm(1)+step(1)*dble(i-1)
     do j=1,np(2) 
          y=xm(2)+step(2)*dble(j-1)
          do k=1,np(3)
              z=xm(3)+step(3)*dble(k-1)
              if (a.eq.1) then                                           
                    write(2,40) Prim(x,y,z,npr)
              else if (a.eq.2) then
                    write(2,40) AO(x,y,z,cao)   
              else if (a.eq.3) then
                    if (uhf) then
                         if (aa.eq.1) write(2,40) MO_a(x,y,z,mo)
                         if (aa.eq.2) write(2,40) MO_b(x,y,z,mo)
                         if (aa.eq.3) then
                              MO_t=0.d0         
                              MO_t=MO_a(x,y,z,mo)+MO_b(x,y,z,mo)
                              write(2,40) MO_t
                         end if
                    else
                         write(2,40) MoOr(x,y,z,mo) 
                    end if 
              else if (a.eq.5) then
                    if (uhf) then
                         if (aa.eq.1) write(2,40) Dens_a(x,y,z)
                         if (aa.eq.2) write(2,40) Dens_b(x,y,z)
                         if (aa.eq.3) then 
                              Dens=0.d0
                              Dens=Dens_a(x,y,z)+Dens_b(x,y,z)
                              write(2,40) Dens
                         end if            
                    else
                         write(2,40) Density(x,y,z)    
                    end if 
              else if (a.eq.6) then !laplacian
                    if (num) then
                        write(2,40) Lapl(x,y,z) !numerically
                    else
                    !write(2,40) Laplacian(x,y,z) !analitically (does not work)
                        write(*,*) "Not available"       
                    end if                
              else   
                    write(*,*) "Error, no cubefile will be generated"
              end if       
          end do
     end do     
end do  

40 format(6(E16.6E3))

end subroutine fwrite


