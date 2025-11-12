!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program wavefunction !main program
!!!!!reads .wfx and .log files and performs calculations about this info
use inputdat  !information about the calculations we want to do
implicit none
integer :: a !defines to subroutine cubefile what function we want to represent: prim, ao, mo, dens
!double precision :: gradx, grady, gradz !gradient
logical :: normalize_dm2p

readwfx=.false. 
readlog=.false. 

call readinput()  

if (readwfx) call filewfx(wfxfilename)  !reads info from a wfx file 
if (readfchk) call filefchk(fchkfilename)
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
  call readintra()
  normalize_dm2p=.true.
  call intracule(normalize_dm2p)  
end if

if (c1calc) then !Becke-Roussel
   call c1hole(70,1.d0)
end if

end program wavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   In this subroutine we read an input file where we specify what type of
!   calculations we want to do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readinput() 
use inputdat
use cubeinfo
use locatemod
implicit none
character*80 :: name
integer :: i
call getarg(1,name)  !gets the name of the input file (writen in the prompt)
name=trim(name)      !remove the blank spaces of the string 'name'
open(unit=3,file=name,status='OLD') 
!set up logical variables as false
primcube=.false.
aocube=.false.
MOcube=.false.
denscube=.false.
gradient=.false.
laplacian=.false.
intracalc=.false.
!read input file
if (located(3,"$wfxfile")) then
    readwfx=.true.
    read(3,*) wfxfilename
else if (located(3,"$fchkfile")) then
    readfchk=.true.
    read(3,*) fchkfilename
else
    write(*,*) "Warning! You must provide at least a .wfx or .fchk file for intracule calculations"
end if
rewind 3
if (located(3,"$basfile")) then
    readbas=.true.
    read(3,*) basfilename
end if    
if (located(3,"$logfile")) then
    readlog=.true.
    read(3,*) logfilename
end if
rewind 3    

if (located(3,"$cubefile")) then
    read(3,*) (center(i), i=1,3) !cube centered in (x,y,z)
    read(3,*) (step(i), i=1,3) !distance between points in the axis
    read(3,*) (np(i), i=1,3) !number of points for each axis
    rewind 3
      ! --- Cubefile printing options ---
    if (located(3,"$density")) then
        denscube = .true.
    end if
    if (located(3,"$primitive")) then
         primcube = .true.
         read(3,*) npr
    end if     
    rewind(3)
    if (readlog .and. located(3,"$AO")) then
         aocube = .true.
         read(3,*) cao
    end if 
    rewind(3)
    if (located(3,"$MO")) then
         MOcube = .true.
         read(3,*) mo
    end if     
    rewind(3)     
    if (located(3,"$gradient")) then
         gradient = .true.
    end if 
    rewind(3)     
    if (located(3,"$laplacian")) then
        laplacian = .true.
    end if
end if
 
if (located(3,"$intracule")) then
    intracalc=.true.
end if
 
if (located(3,"$HFhole")) then
    c1calc=.true.
end if
end subroutine readinput 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readintra()
!read information from the intracule function     
use intrainfo
use locatemod
use quadratures
implicit none
character*80 :: name
integer :: i, j
character*20 :: ccent
call getarg(1,name)  !gets the name of the input file (writen in the prompt)
name=trim(name)      !remove the blank spaces of the string 'name'
open(unit=3,file=name,status='OLD')  
        !read input for intracule
        call locate(3,'$Integral screening threshold') !threshold for integral screenings
        read(3,*) thresh  
        rewind(3)
        dm2name = ''
        call locate(3,'$DM2P')
        read(3,*) dm2name           !name of the dm2p file
        read(3,*) trsh1, trsh2      !thresholds used in DM2prim
        read(3,*) outname           !name of the output file
        if (located(3,'$nosym')) then
            nosym=.true.
            write(*,*) 'Deactivating symmetry in intracule calculations'
        else
            nosym=.false.
        end if
        rewind(3)
        if (located(3,'$radial_integral')) then
            radial_integral=.true.
            if (located(3,'$definite')) then
                definite=.true.
                read(3,*) a,b
            end if
            if (located(3,'$Multicenter')) then
                if (located(3,'$manual_grid')) then !number of quadrature centers(not automatically calculated)
                    read(3,*) nquad
                    allocate(cent(3,nquad))
                    do i=1,nquad
                        read(3,*) (cent(j,i), j=1,3) !read position of centers
                    end do
                    allocate(nradc(nquad))
                    allocate(nangc(nquad))
                    allocate(sfalpha(nquad)) 
                    call locate(3,'$Gauss-Legendre')
                    read(3,*) (nradc(i), i=1,nquad)
                    read(3,*) sfalpha(:)
                    call locate(3,'$Gauss-Lebedev')
                    read(3,*) (nangc(i), i=1,nquad)
                    call locate(3,'$Center weights')
                    allocate(Ps(nquad)) !allocate weight of each centre
                    read(3,*) Ps(:) !the weigth of a positive center must be equal to the neg.
                    rewind(3)
                else !automatically calculate the parameters (default)
                    if (located(3,'$Gauss-Legendre')) then
                        read(3,*) nrad
                    else
                        nrad=50
                    end if
                    if (located(3,'$Gauss-Lebedev')) then
                        read(3,*) nang
                    else
                        nang=590
                    end if
                    if (located(3,'$hydro')) then
                        nohydro=.false.
                    else
                        nohydro=.true. !do not include H atoms as centres per default
                    end if
                    call centercalc() !calculate the number of centers and their positions
                end if
            else !single center quadrature
                write(*,*) 'Using single center quadrature'
                nquad=1
                allocate(cent(3,1))
                allocate(nradc(1)); allocate(nangc(1)); allocate(sfalpha(1)); allocate(Ps(nquad))
                cent(1:3,1)=0.0d0
                if (located(3,'$Gauss-Legendre')) then
                    write(*,*) 'Reading number of radial points and scaling factor for radial integration'
                    read(3,*) nradc(1)
                    read(3,*) sfalpha(1)
                else !default values
                    nradc(1)=50
                    sfalpha(1)=1.0d0
                end if
                if (located(3,'$Gauss-Lebedev')) then
                    read(3,*) nangc(1)
                else !default values
                    write(*,*) 'Using default value of 590 points for angular integration'
                    nangc(1)=590
                end if
            end if    
        else if (located(3,'$radial_plot')) then
            radial_plot=.true.
            read(3,*) r_plot_name
            read(3,*) nblock
            allocate(n_an_per_part(nblock))
            allocate(tart(2,nblock))
            allocate(stp(nblock))
            do i=1,nblock
                read(3,*) tart(:,i), stp(i), n_an_per_part(i)
            end do
        else if (located(3,'$Vectorial_plot')) then     
            cubeintra=.true.  
            read(3,*) cubeintraname
            read(3,*) (center_i(i), i=1,3) !cube centered in (x,y,z)
            read(3,*) (step_i(i), i=1,3) !distance between points in the axis
            read(3,*) (np_i(i), i=1,3) !number of points for each axis
        else if (located(3,'$Intracule_at_zero')) then
            intracule_at_zero=.true.
        else if (located(3,'$Intracule_two_points')) then
            intracule_two_points=.true.   
            read(3,*) x_point1, y_point1, z_point1
            read(3,*) x_point2, y_point2, z_point2
        else
            write(*,*) 'Warning! You must provide at least a radial or vectorial plot option for intracule calculations'    
        end if
        rewind(3)
 close(3) 
end subroutine readintra





