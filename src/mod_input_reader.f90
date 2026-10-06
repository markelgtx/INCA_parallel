module input_reader
   use inputdat   
   use intrainfo  
   use quadratures 
   use cubefile   
   use autogrid
   implicit none
   
   ! --- Temporary static arrays for NAMELIST reading ---
   ! Max 20 blocks for radial scan
   real(8), dimension(20) :: scan_start = 0.0d0
   real(8), dimension(20) :: scan_end   = 0.0d0
   real(8), dimension(20) :: scan_step  = 0.0d0
   integer, dimension(20) :: scan_nang  = 590
   
   ! Max 50 centers for multicenter manual grid
   integer :: nquad_in = 1
   real(8), dimension(3,50) :: cent_in = 0.0d0
   integer, dimension(50) :: nradc_in = 50, nangc_in = 590
   real(8), dimension(50) :: sfalpha_in = 1.0d0, Ps_in = 1.0d0

   ! --- Module Variables ---
   ! Moved here so setup_grids can see them
   logical :: multicenter = .false.  
   logical :: manual_grid = .false.

contains

   subroutine read_all_input(filename)
      character(len=*), intent(in) :: filename
      integer :: i, ios, dot_idx      
      
      ! ==========================================================
      ! Define the NAMELIST groups
      ! ==========================================================
      namelist /files/ wfxfilename, fchkfilename, basfilename, &
                       logfilename, dm2name, outname, cubeintraname, r_plot_name

      namelist /intracule_job/ intracalc, c1calc, radial_integral, radial_plot, &
                               cubeintra, intracule_at_zero, intracule_two_points, &
                               thresh, nosym, x_point1, y_point1, z_point1, &
                               x_point2, y_point2, z_point2, multicenter, &
                               manual_grid, nrad, nang, nquad_in, cent_in, &
                               nradc_in, nangc_in, sfalpha_in, Ps_in, nblock, &
                               scan_start, scan_end, scan_step, scan_nang, &
                               center_i, step_i, np_i

      namelist /cube_job/ cube, denscube, primcube, aocube, mocube, gradient, &
                          laplacian, npr, cao, mo, center, step, np
      
      ! ==========================================================
      ! Set Defaults
      ! ==========================================================
      thresh = 1.0d-8
      intracalc = .false.; c1calc = .false.; cube = .false.
      denscube = .false.; primcube = .false.; aocube = .false.; mocube = .false.
      gradient = .false.; laplacian = .false.      
      radial_integral = .false.; radial_plot = .false.; cubeintra = .false.
      intracule_at_zero = .false.; intracule_two_points = .false.
      definite = .false.; multicenter = .true.; manual_grid = .false.
      betaone = .false.; nohydro = .true.; nosym = .false.
      nrad = 50; nang = 590; nblock = 1   
      wfxfilename = ""; fchkfilename = ""; basfilename = ""; logfilename = ""
      dm2name = ""; outname = ""; cubeintraname = ""; r_plot_name = ""
      ! ==========================================================
      ! Read the file
      ! ==========================================================
      open(unit=10, file=filename, status='old', iostat=ios)
      if (ios /= 0) then
         write(*,*) "ERROR: Could not open input file: ", trim(filename)
         stop
      end if

      open(unit=10, file=filename, status='old')
      
      read(10, nml=files)
      rewind(10)
      
      read(10, nml=intracule_job)
      rewind(10)
      
      read(10, nml=cube_job)
      close(10)

      ! Read blocks and check for syntax errors (ios > 0 means syntax error)
      read(10, nml=files, iostat=ios)
      if (ios > 0) write(*,*) "WARNING: Syntax error in &files block!"
      rewind(10)
      
      read(10, nml=intracule_job, iostat=ios)
      if (ios > 0) write(*,*) "WARNING: Syntax error in &intracule_job block!"
      rewind(10)
      
      read(10, nml=cube_job, iostat=ios)
      if (ios > 0) write(*,*) "WARNING: Syntax error in &cube_job block!"
      close(10)

      ! ==========================================================
      ! Post-processing string and logicals
      ! ==========================================================
      readwfx  = (len_trim(wfxfilename) > 0)
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
      
      if (len_trim(cubeintraname) == 0) then
         dot_idx = index(trim(filename), '.', back=.true.)
         if (dot_idx > 0) then
             cubeintraname = filename(1:dot_idx-1) // "_vector.cube"
         else
             cubeintraname = trim(filename) // "_vector.cube"
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
      
   end subroutine read_all_input

   ! ==========================================================
   ! Called AFTER the molecule is loaded
   ! ==========================================================
   subroutine setup_grids()
      integer :: i ! <-- Fixed: 'i' needed to be declared here for the loop
      
      if (radial_integral) then
         if (multicenter) then
            if (definite) then
               ! Fixed line truncation error
               write(*,*) "ERROR: Multicenter radial integrals not implemented " // &
                          "for finite limits. Please set definite=.false. or multicenter=.false."
               stop
            end if
            if (manual_grid) then  
               ! Load manual parameters from namelist static arrays
               nquad = nquad_in
               allocate(cent(3,nquad), nradc(nquad), nangc(nquad), sfalpha(nquad), Ps(nquad))
               do i=1,nquad
                  cent(1:3, i) = cent_in(1:3, i)
                  nradc(i)     = nradc_in(i)
                  nangc(i)     = nangc_in(i)
                  sfalpha(i)   = sfalpha_in(i)
                  Ps(i)        = Ps_in(i)
               end do
            else
               call centercalc()
            end if   
         else   
            ! Single center radial integral
            nquad = 1
            allocate(cent(3,1), nradc(1), nangc(1), sfalpha(1), Ps(1))
            cent(1:3,1) = 0.0d0
            nradc(1) = nradc_in(1)   ! Falls back to 50 if not set
            nangc(1) = nangc_in(1)   ! Falls back to 590 if not set
            sfalpha(1) = sfalpha_in(1)
            Ps(1) = 1.0d0
         end if
      end if
   end subroutine setup_grids

end module input_reader