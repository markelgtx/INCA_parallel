module geninfo           !general information, can be obtained from .wfx files.
  implicit none
  logical :: uhf, rhf    !unrestricted-restricted Hartree Fock???
  logical :: rdens,udens !relaxed or unrelaxed density
  logical :: opsh, clsh  !open or closed shell
  integer :: natoms      !number of atoms
  integer :: nelec
  integer :: nprim       !number of primitives
  double precision :: netch !net charge
  integer :: mspin       !electronic spin multiplicity
  integer :: nalfae      !number of alpha electrons
  integer :: nbetae      !number of beta electrons
  integer, allocatable, dimension(:) :: Ra    !atomic centers of primitives
  integer, allocatable, dimension(:) :: Ptyp  !primitive type
  double precision, allocatable, dimension(:) :: Alpha !Primitive exponents
  integer, allocatable, dimension(:,:) :: TMN !matrix with t,m,n coeficients
  double precision, allocatable,  dimension(:,:) :: cartes !nuclear cartesian coordinates
  integer, allocatable, dimension(:) :: an !atomic number
  double precision, allocatable, dimension(:) :: chrg !nuclear charge  
  double precision, parameter :: pi=4.d0*datan(1.d0)  
end module

