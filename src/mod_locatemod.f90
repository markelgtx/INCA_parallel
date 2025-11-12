module locatemod
  implicit none
contains
  !--------------------------------------------
  subroutine locate(iunit,string)
  !--------------------------------------------
  ! Search for string in unit. If not found, abort.
  integer, intent(in) :: iunit
  character(*), intent(in) :: string
  character(len=200) :: linia
  integer :: ios

  rewind(iunit)
  do
     read(iunit,"(A)",iostat=ios) linia
     if (ios /= 0) exit
     if (index(linia,string) /= 0) return
  end do

  write(*,*) "ERROR: Section not found -> ", trim(string)
  stop 1
  end subroutine locate

  !--------------------------------------------
  logical function located(iunit,string)
  !--------------------------------------------
  ! Search for string in unit. Returns .true. if found.
  integer, intent(in) :: iunit
  character(*), intent(in) :: string
  character(len=200) :: linia
  integer :: ios

  rewind(iunit)
  located = .false.
  do
     read(iunit,"(A)",iostat=ios) linia
     if (ios /= 0) exit
     if (index(linia,string) /= 0) then
        located = .true.
        write(*,*) "Located section -> ", trim(string)
        return
     end if
  end do
  end function located

  !--------------------------------------------
  subroutine getline_with(iunit,string,line)
   !--------------------------------------------
   ! Search for string in unit.
   ! If found, returns the line containing it.
   ! Leaves file pointer positioned *after* that line.
   integer, intent(in) :: iunit
   character(*), intent(in) :: string
   character(len=*), intent(out) :: line
   character(len=500) :: linia
   integer :: ios
 
   rewind(iunit)
   do
      read(iunit,"(A)",iostat=ios) linia
      if (ios /= 0) exit
      if (index(linia,string) /= 0) then
         line = linia
         return
      end if
   end do
 
   write(*,*) "ERROR: Section not found -> ", trim(string)
   stop 1
   end subroutine getline_with
end module locatemod
