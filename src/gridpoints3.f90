 subroutine gridpoints3(center_i,step_i,np_i)
  use quadratures
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
   write(*,*) "xm", xm(:)
   write(*,*) "np_i", np_i(:)
   do i=1,rgrid
     write(7,*) rg(:,i)
   end do  
 end subroutine gridpoints3    

 subroutine gridpoints3_corner(center_i, step_i, np_i)
   use quadratures
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
 
   ! Diagnostic output
   write(*,*) "Grid corner (origin):", xm(:)
   write(*,*) "Number of points per direction:", np_i(:)
   do i = 1, rgrid
      write(7,*) rg(:, i)
   end do
 
 end subroutine gridpoints3_corner