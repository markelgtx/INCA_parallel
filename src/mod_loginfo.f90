 module loginfo
   implicit none 
   double precision, allocatable, dimension(:) :: Flg, N_prim !fixed expansion coeficient for the contraction of AOs into primitives, normalization constant for the primitives
   double precision, allocatable, dimension(:,:) :: Ckalk !expansion coeficient of MOs in AOs
   integer, allocatable, dimension(:) :: npao !number of primitives for each atomic orbital
   integer :: nao !number of atomic orbitals in the molecule
   integer, allocatable, dimension(:) :: aotyp !type of atomic orbital: s->0, px, py, pz ->1,  ...   
end module loginfo   

