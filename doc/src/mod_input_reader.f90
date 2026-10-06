module input_reader
   use inputdat   
   use intrainfo  
   use quadratures 
   use cubefile   
   use autogrid
   implicit none
   
   ! --- Temporary static arraays for NAMELIST reading ---
   ! Max 20 blocks for radial scan
   real(8), dimension(20) :: scan_start = 0.0d0
   real(8), dimension(20) :: scan_end   = 0.0d0
   real(8), dimension(20) :: scan_step  = 0.0d0
   integer, dimension(20) :: scan_nang  = 590
   ! Max 50 centers for multicenter manual grid
   integer :: nquad_in=1
   real(8), dimension(3,50) :: cent_in = 0.0d0
   integer, dimension(50) :: nradc_in = 50, nangc_in = 590
   real(8), dimension(50) :: sfalpha_in = 1.0d0, Ps_in = 1.0d0

   ! --- Module Variables ---
   logical :: multicenter = .false.  ! Moved here so setup_grids can see it

contains

   subroutine read_all_input(filename)
      character(len=*), intent(in) :: filename
      integer :: i, ios, dot_idx      
      ! --- Define the comprehenive NAMELIST variables ---
      namelist /input_files/ &
      ! Strings / Files
         wfxfilename, fchkfilename, basfilename, logfilename, dm2name
      namelist /output_files/ &
         outname, cubeintraname, r_plot_name, &
      ! Main Job Flags
         intracalc, c1calc, cube, &
      !Standard Cube Options   
         denscube, primcube, aocube, mocube, gradient, laplacian, &
         npr, cao, mo, center, step, np, &
      ! Intracule calculation options
      namelist /intracule_job/ &   
         radial_integral, radial_plot, cubeintra, &
         intracule_at_zero, intracule_two_points, &
      !  Primitive quartet options
         thresh, nosym, &
      !  Two points for the intracule
         x_point1, y_point1, z_point1, x_point2, y_point2, z_point2, &
      !  Radial integral grid configuration flags
         multicenter, manual_grid, nrad, nang, &
      !  Manual gtid arrays
         nquad_in, cent_in, nradc_in, nangc_in, sfalpha_in, Ps_in, &
      !  Radial Scan Block Arrays
         nblock, scan_start, scan_end, scan_step, scan_nang, & 
      !  Vectorial plot options
         center_i, step_i, np_i,            
      
      ! ==========================================================
      ! Set Defaults
      ! ==========================================================
      thresh = 1.0d-8
      intracalc = .false.; c1calc = .false.; cube = .false.
      denscube = .false.; primcube = .false.; aocube = .false.; mocube = .false.; gradient = .false.; laplacian = .false.      
      radial_integral = .false.; radial_plot = .false.; cubeintra = .false.
      intracule_at_zero = .false.; intracule_two_points = .false.
      definite = .false.; multicenter = .true.; manual_grid = .false.
      betaone = .false.; nohydro= .true.; nosym = .false.
      nrad=50; nang=590; nblock=1   

      !Read the file
      open(unit=10, file=filename, status='old', iostat=ios)
      if (ios /= 0) then
         write(*,*) "ERROR: Could not open input file: ", trim(filename)
         stop
      end if
      read(10, nml=input_files, iostat=ios)
      if (ios /= 0) then
         write(*,*) "ERROR: Malformed NAMELIST in ", trim(filename)
         stop
      end if
      close(10)

      ! --- Post-processing string and logicals ---
      readwfx = (len_trim(wfxfilename) > 0)
      readfchk = (len_trim(fchkfilename) > 0)
      readbas  = (len_trim(basfilename) > 0)
      readlog  = (len_trim(logfilename) > 0)

      if (.not.readwfx .and. .not.readfchk) then
         write(*,*) "ERROR: No wavefunction file specified. Please provide a .wfx or .fchk file."
         stop
      end if  

      ! --- Smart output Filename Generation ---
      if (len_trim(outname) == 0) then
         dot_idx = index(trim(filename), '.', back=.true.)
         if (dot_idx > 0) then
            outname = filename(1:dot_idx-1) // ".out"
         else
            outname = trim(filename) // ".out"
         end if
      end if

      if (len_trim(r_plot_name) == 0) then
         dot_idx = index(trim(filename), '.', back=.true.)
         if (dot_idx > 0) then
             r_plot_name = filename(1:dot_idx-1) // ".rad"
         else
             r_plot_name = trim(filename) // ".rad"
         end if
      end if
      
      ! --- Allocate Radial Scan Dynamic Arrays ---
      if (radial_plot) then
         allocate(n_an_per_part(nblock))
         allocate(tart(2, nblock))
         allocate(stp(nblock))
         do i=1,nblock
            tart(1, i) = scan_start(i)
            tart(2, i) = scan_end(i)
            stp(i) = scan_step(i)
            n_an_per_part(i) = scan_nang(i)
         end do
      end if
      !character(len=20) :: calc_type = 'none'                                
      !wfxfilename = ''
      !fchkfilename = ''
      !basfilename = ''
      !logfilename = ''
      !dm2name = ''
      !outname = ''
      !cubeintraname = ''
      !r_plot_name = ''
      ! ==========================================================
      ! --- Open and read the file ---
      !open(unit=10, file=filename, status='old', action='read')
      !read(10, nml=files, iostat=ios)
      !if (ios /= 0 .and. ios > 0) write(*,*) "WARNING: Formatting issue in &files namelist!"
      !read(10, nml=intracule_job, iostat=ios)
      !if (ios /= 0 .and. ios > 0) write(*,*) "WARNING: Formatting issue in &intracule_job namelist!"
      !close(10)      
      ! Restoring your default outname logic
      !if (len_trim(outname) == 0) then
      !   outname = trim(filename) // '.out'
         ! --- Clean Master Output Filename ---
      !   block
      !      integer :: dot_idx
      !      dot_idx = index(trim(filename), '.', back=.true.)
      !      if (dot_idx > 0) then
      !         outname = filename(1:dot_idx-1) // ".out"
      !      else
      !         outname = trim(filename) // ".out"
      !      end if
      !   end block
      !end if
      ! --- Translation Layer: Map defaults and logicals ---
      !readwfx  = (len_trim(wfxfilename) > 0)
      !readfchk = (len_trim(fchkfilename) > 0)
      !readbas  = (len_trim(basfilename) > 0)
      !readlog  = (len_trim(logfilename) > 0)      
      !if (len_trim(calc_type) > 0 .and. trim(calc_type) /= 'none') intracalc = .true.      
      !radial_plot = (trim(calc_type) == 'radial_plot')
      !radial_integral = (trim(calc_type) == 'radial_integral')
      !cubeintra = (trim(calc_type) == 'vectorial_plot')  
      ! --- Allocate and populate dynamic arrays for plots ---
      !if (radial_plot) then
      !   allocate(n_an_per_part(nblock))
      !   allocate(tart(2, nblock))
      !   allocate(stp(nblock))  
      !   do i = 1, nblock
      !      tart(1, i) = tmp_start(i)
      !      tart(2, i) = tmp_end(i)
      !      stp(i) = tmp_step(i)
      !      n_an_per_part(i) = tmp_nang(i)
      !   end do
      !end if
      
   end subroutine read_all_input
   ! ==========================================================
   ! Called AFTER the molecule is loaded
   ! ==========================================================
   subroutine setup_grids()
      if (radial_integral) then
         if (multicenter) then
            if (definite) then
               write(*,*) "ERROR: Multicenter radial integrals not implemented for finite limits. Please set definite=.false. or multicenter=.false."
               stop
            end if
            if (manual_grid) then  
               !Load manual parameters from namelist static arrays
               nquad = nquad_in
               allocate(cent(3,nquad), nradc(nquad), nangc(nquad), sfalpha(nquad), Ps(nquad))
               do i=1,nquad
                  cent(1:3, i) = cent_in(1:3, i)
                  nradc(i)     = nradc_in(i)
                  nangc(i)     = nangc_in(i)
                  sfalpha(i)   = sfalpha_in(i)
                  Ps(i)        = Ps_in(i)
               end do
            end if   
         else   
            call centercalc()
         end if   
      else
         ! Single center radial integral. Just set up one center at the origin with the specified nrad and nang
         nquad = 1
         allocate(cent(3,1), nradc(1), nangc(1), sfalpha(1), Ps(1))
         cent(1:3,1) = 0.0d0
         nradc(1) = nradc_in(1)   ! Falls back to 50 if not set
         nangc(1) = nangc_in(1)   ! Falls back to 590 if not set
         sfalpha(1) = sfalpha_in(1)
         Ps(1) = 1.0d0
      end if
   end subroutine setup_grids

end module input_reader