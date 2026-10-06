title: User Manual

# INCA User Manual

INCA uses a standard Fortran `NAMELIST` input structure. The input file must begin with `&input` and end with `/`. Variables are order-independent, case-insensitive, and separated by commas or newlines. 

## 1. General Configuration

### File I/O
* **`wfxfilename`** *(string)*: Path to the input `.wfx` file.
* **`fchkfilename`** *(string)*: Path to the input `.fchk` file.
* **`basfilename`** *(string)*: Path to custom basis set `.bas`.
* **`dm2name`** *(string)*: Path to `.dm2p` stream file.
* **`outname`** *(string)*: Custom master output log filename.

### Task Flags
* **`intracalc`** *(logical)*: Triggers the intracule evaluation engine.
* **`cube`** *(logical)*: Triggers standard 1-electron 3D cubefile generation.
* **`c1calc`** *(logical)*: Triggers Exchange-Correlation hole evaluation.

---

## 2. Intracule Options (`intracalc = .true.`)

### Calculation Types
* **`radial_integral`**: Computes the total integrated radial intracule and $V_{ee}$.
* **`radial_plot`**: Evaluates the radial intracule $I(s)$ over a 1D scan grid.
* **`cubeintra`**: Evaluates the vectorial intracule on a 3D grid.
* **`intracule_at_zero`**: Evaluates the probability solely at $s = 0$.
* **`intracule_two_points`**: Evaluates the vector intracule between exactly two coordinates. (Requires `x_point1`, `y_point1`, `z_point1` and `x_point2`, `y_point2`, `z_point2`).

### General Grid & Math Parameters
* **`thresh`** *(real)*: Threshold for Cioslowski-Liu integral screening (Default: `1.0d-12`).
* **`nosym`** *(logical)*: Disables spatial symmetry optimization.


*Note: Angular integrations utilize high-precision Lebedev spherical quadratures.*

### Radial Plot Block Settings
Used when `radial_plot = .true.`. Allows for varying grid resolutions (e.g., fine grid near $s=0$, coarse grid at long range).
* **`nblock`** *(integer)*: Number of scan regions.
* **`scan_start(i)`**, **`scan_end(i)`**, **`scan_step(i)`** *(real)*: Boundaries and step size for block `i`.
* **`scan_nang(i)`** *(integer)*: Number of Lebedev angular points for block `i` (e.g., 6, 110, 590, 1202).

---

## 3. Advanced Grid Topologies (For Developers)

Used when `radial_integral = .true.`. 
By default, INCA uses a **Single Center** expansion at the origin. 

### Automatic Multicenter Grid
* **`multicenter = .true.`**: Enables automated multicenter partitioning (Salvador TFVC scheme).
* **`nrad`**, **`nang`** *(integers)*: Global radial and angular points per center.
* **`nohydro`** *(logical)*: Excludes Hydrogen atoms as integration centers (Default: `.true.`).
* **`betaone`** *(logical)*: Forces alternate Becke partition stiffness.

### Manual Multicenter Grid
* **`multicenter = .true.`** AND **`manual_grid = .true.`**: Allows explicit definition of integration centers.
* **`nquad_in`** *(integer)*: Total number of centers.
* **`cent_in(1:3, i)`** *(real)*: Cartesian X, Y, Z for center `i`.
* **`nradc_in(i)`**, **`nangc_in(i)`** *(integers)*: Grid points specific to center `i`.
* **`sfalpha_in(i)`**, **`Ps_in(i)`** *(reals)*: Scaling factors and center weights.

### Definite Integrals
* **`definite = .true.`**: Evaluates the integral only within bounds `a` and `b` (Must be single-center).

---

## 4. Cubefile Options (`cube = .true.`)

* **`denscube`**, **`primcube`**, **`aocube`**, **`mocube`**, **`gradient`**, **`laplacian`** *(logicals)*: Selects the property to export.
* **`npr`**, **`cao`**, **`mo`** *(integers)*: Selects the specific primitive, AO, or MO index to plot.
* **`center(1:3)`**, **`step(1:3)`**, **`np(1:3)`**: Defines the 3D bounding box origin, step size, and number of voxels per axis.
