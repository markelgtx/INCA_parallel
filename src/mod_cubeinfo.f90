module cubeinfo
   implicit none
   integer :: mo !molecular orbital to represent in the cubefile
   integer :: npr  !primitive to represent in the cubefile
   integer :: cao  !atomic orbital to represent in the cubefile
   double precision, dimension(3) :: center !cube centered in (x,y,z)
   double precision, dimension(3) :: step !distance between points for each axis
   integer, dimension(3) :: np !number of points in the cube for each axis
end module cubeinfo

