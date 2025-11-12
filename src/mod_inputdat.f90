module inputdat
implicit none
 logical :: readwfx, readfchk, readlog, readbas, cube, primcube, aocube, MOcube, denscube, gradient, laplacian, intracalc
 logical :: c1calc
 character*40 :: wfxfilename
 character*40 :: fchkfilename
 character*40 :: logfilename
 character*40 :: basfilename
 character*40 :: nameprim
 character*40 :: nameao
 character*40 :: namemo   !name of the cubefile we will generate
 character*40 :: namedens
 character*40 :: namegrad
 character*40 :: namelap
 character*40 :: nameint
 character*40 :: namec1
end module inputdat

