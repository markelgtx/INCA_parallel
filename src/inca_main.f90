!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program inca !main program
!!!!!reads .wfx and .log files and performs calculations about this info
use inputdat  !information about the calculations we want to do
implicit none
integer :: a !defines to subroutine cubefile what function we want to represent: prim, ao, mo, dens
logical :: normalize_dm2p

call readinput()  

if (readwfx) call filewfx(wfxfilename)  !reads info from a wfx file 
if (readfchk) call filefchk(fchkfilename) !reads info from a fchk file
if (readbas) call filebas(basfilename)  !reads info from a bas file
if (readlog) call filelog(logfilename)  !reads info from a log file

if (cube) then
  if (primcube) then  
      a=1
      call cubefile(a,nameprim) !Generate cubefile with a Primitive
  end if
  if (aocube) then 
      a=2
      call cubefile(a,nameao)  !Generate cubefile with AO
  end if
  if (mocube) then 
      a=3
      call cubefile(a,namemo) !Generate cubefile with a MO 
  end if  
  if (denscube) then
      a=5
      call cubefile(a,namedens) !Generate cubefile with density from MO
  end if
  !if (gradient) then 
  !  call gradient(0,0,0,gradx,grady,gradz) !not impletented
  !end if
  if (laplacian) then 
     a=6     
     call cubefile(a,namelap) !generate a cubefile with laplacian
  end if
end if

if (intracalc) then !compute the intracule
  if (.not.readbas) then
    write(*,*) "CAUTION: Primitive info comes from .wfx or .fchk file"
    write(*,*) " i primtive may not coincide with that of the .dm2p file"
  end if
  call readintra() !read input information about intracule
  normalize_dm2p=.true. 
  call intracule(normalize_dm2p)  
end if

if (c1calc) then !Becke-Roussel (on halt)
  call c1hole(70,1.d0)
end if

end program inca







