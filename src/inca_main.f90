!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program inca 
  !! Main driver program for INCA.
  !! Reads .wfx and .log files and performs intracule calculations.
  use inputdat    ! global flags
  use intrainfo   ! intracule specific parameters
  use input_reader! namelist parsing
  use read_files
  use cubefile
  use intracule
  use c1hole
  implicit none

  character(len=80) :: input_file_name
  integer :: cube_mode ! defines to subroutine cubefile what function we want to represent
  logical :: normalize_dm2p

  ! ==========================================================
  ! SAFETY CHECK: Did the user provide an input file?
  ! ==========================================================
  if (command_argument_count() == 0) then
     write(*,*) "====================================================="
     write(*,*) " INCA: Intracule Calculator"
     write(*,*) "====================================================="
     write(*,*) " ERROR: No input file provided."
     write(*,*) " USAGE: ./inca.exe <input_file.txt>"
     write(*,*) "====================================================="
     stop
  end if

  ! get the name of the input file from the command line 
  call getarg(1,input_file_name) 
  input_file_name = trim(input_file_name)

  ! call the routine to read the input file and populate variables 
  call read_all_input(input_file_name)

  if (readwfx) call filewfx(wfxfilename)  
  if (readfchk) call filefchk(fchkfilename) 
  if (readbas) call filebas(basfilename)  
  if (readlog) call filelog(logfilename)  

  if (cube) then
    if (primcube) then  
        cube_mode=1
        call cubegen(cube_mode,nameprim) 
    end if
    if (aocube) then 
        cube_mode=2
        call cubegen(cube_mode,nameao)  
    end if
    if (mocube) then 
        cube_mode=3
        call cubegen(cube_mode,namemo) 
    end if  
    if (denscube) then
        cube_mode=5
        call cubegen(cube_mode,namedens) 
    end if
    if (laplacian) then 
       cube_mode=6     
       call cubegen(cube_mode,namelap) 
    end if
  end if

  if (intracalc) then 
    if (.not.readbas) then
      write(*,*) "CAUTION: Primitive info comes from .wfx or .fchk file"
      write(*,*) " i primtive may not coincide with that of the .dm2p file"
    end if
    
    ! Set up the grids now that atomic numbers are loaded 
    call setup_grids()
    
    normalize_dm2p=.true.
    call calculate_intracule(normalize_dm2p)  
  end if

  if (c1calc) then  !Becke Roussel approximation
    !call c1hole(70,1.d0)
    call c1holeapp(70,1.d0)  !70 is the number of points in the grid, 1.d0 is the scaling factor for the grid
  end if
end program inca