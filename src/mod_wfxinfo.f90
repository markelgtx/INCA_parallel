module wfxinfo !especific data of wfx files 
  implicit none
  logical :: corr   !correlated method
  integer :: noccmo !number of occupied MOs
  double precision :: toteng
  double precision, allocatable, dimension(:,:) :: T   !MOs in primitives coeficients (matriu T)
  double precision, allocatable, dimension(:,:) :: T_a !Alpha MOs coeficients
  double precision, allocatable, dimension(:,:) :: T_b !Beta MOs coeficients
  double precision, allocatable, dimension(:) :: Occ   !ocupancies of orbitals
end module wfxinfo

