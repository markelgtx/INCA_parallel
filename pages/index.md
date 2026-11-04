title: User Manual

# INCA User Manual

INCA uses a standard Fortran `NAMELIST` input structure. The input file contains three namelist blocks, `&files`, `&intracule_job` and `&cube_job` (each ended by `/`; a block that is left out keeps its defaults, and an unknown keyword stops the run with an error). Variables are order-independent, case-insensitive, and separated by commas or newlines. 

## 1. General Configuration

### File I/O
* **`wfxfilename`** *(string)*: Path to the input `.wfx` file.
* **`fchkfilename`** *(string)*: Path to the input `.fchk` file.
* **`basfilename`** *(string)*: Path to the primitive basis file `.bas`.
* **`dm2name`** *(string)*: Path to the `.dm2p` file (DM2 in the primitive basis).

**Intracule calculations need three files:** the `.fchk` (or `.wfx`), the `.bas` and the `.dm2p`. INCA does not generate the `.dm2p`; it has to be produced by an external program (for example DMN).

* The **`.bas`** file holds the primitive Gaussian basis, one line per primitive after a header line starting with `#`: atomic number of the centre, centre coordinates X, Y, Z (bohr), exponent, and the Cartesian angular powers nx, ny, nz.
* The **`.dm2p`** file holds the second-order reduced density matrix in the primitive basis: one binary record per primitive quartet, with the four primitive indices i, j, k, l (4-byte integers) and the DM2 value for that quartet (8-byte real, with respect to normalized primitives).
* Both files must use the **same primitive ordering**: primitive *n* in the `.dm2p` is the primitive on line *n* of the `.bas`. Without a `.bas`, INCA takes the primitives from the `.fchk` / `.wfx`, whose order does not in general match the `.dm2p`.

The full format description, with an example, is in the README.
* **`outname`** *(string)*: Custom master output log filename.

### Task Flags
* **`calc_type`** *(string)*: Triggers the intracule evaluation engine; selects `radial_integral`, `radial_plot` or `cubeintra` (see below).
* **`cube`** *(logical)*: Triggers standard 1-electron 3D cubefile generation.
* **`c1calc`** *(logical)*: Triggers Exchange-Correlation hole evaluation.

---

## 2. Intracule Options (`calc_type`, `intracule_at_zero`, `intracule_two_points`)

### Calculation Types
* **`calc_type = 'radial_integral'`**: Computes the total integrated radial intracule and $V_{ee}$.
* **`calc_type = 'radial_plot'`**: Evaluates the radial intracule $I(s)$ over a 1D scan grid.
* **`calc_type = 'cubeintra'`**: Evaluates the vectorial intracule on a 3D grid.
* **`intracule_at_zero`** *(logical)*: Evaluates the intracule solely at $s = 0$. It is an alternative to `calc_type` and does not need it. The result is written to the `.out` file as `INTRACULE AT ZERO` (the raw value divided by two).
* **`intracule_two_points`** *(logical)*: Evaluates the vector intracule at exactly two coordinates. It is an alternative to `calc_type` and does not need it. Requires `x_point1`, `y_point1`, `z_point1` and `x_point2`, `y_point2`, `z_point2` (bohr). The values are written to the `.out` file as `Point 1` and `Point 2` (raw values, not divided by two, so the value at the origin is twice `INTRACULE AT ZERO`).

Example, intracule at two points:

```fortran
&intracule_job
  intracule_two_points = .true.,
  x_point1 = 0.0d0, y_point1 = 0.0d0, z_point1 = 0.0d0,
  x_point2 = 1.0d0, y_point2 = 0.0d0, z_point2 = 0.0d0,
  thresh = 1.0d-12, multicenter = .false., nrad = 50, nang = 110
/
```

### General Grid & Math Parameters
* **`thresh`** *(real)*: Threshold for Cioslowski-Liu integral screening (Default: `1.0d-8`).
* **`nosym`** *(logical)*: Intracule symmetry $I(s)=I(-s)$. By default it is used and only half of the grid points are computed; `nosym = .true.` disables it and uses all points.


*Note: Angular integrations utilize high-precision Lebedev spherical quadratures.*

### Radial Plot Block Settings
Used when `calc_type = 'radial_plot'`. Allows for varying grid resolutions (e.g., fine grid near $s=0$, coarse grid at long range).
* **`nblock`** *(integer)*: Number of scan regions.
* **`scan_start(i)`**, **`scan_end(i)`**, **`scan_step(i)`** *(real)*: Boundaries and step size for block `i`.
* **`scan_nang(i)`** *(integer)*: Number of Lebedev angular points for block `i` (e.g., 6, 110, 590, 1202).

---

## 3. Advanced Grid Topologies (For Developers)

Used when `calc_type = 'radial_integral'`. 
By default, INCA uses the automatic **multicenter** grid (`multicenter = .true.`). We recommend always stating the grid explicitly in the input: `multicenter = .true.` for the multicenter grid, or `multicenter = .false.` for a single-center expansion at the origin.

### Automatic Multicenter Grid
* **`multicenter = .true.`**: Enables automated multicenter partitioning (Salvador TFVC scheme). This is the default; set `multicenter = .false.` for a single-center expansion at the origin.
* **`nrad`**, **`nang`** *(integers)*: Radial and angular points per center (Defaults: `50` and `590`). They are also used by the single-center grid (`multicenter = .false.`), whose radial scaling factor is `sfalpha_in(1)` (Default: `1.0` bohr).
* The radial scaling factor of each center is computed automatically from Bragg-Slater atomic radii (Slater, J. Chem. Phys. 41, 3199 (1964), with H = 0.35 Å as in Becke 1988; H-Kr are tabulated, other elements use 1.0 Å): the average of the radii of the two atoms of a pair for interatomic centers, and the average weighted by the number of intra-atomic pairs for the origin.
* **`nohydro`** *(logical)*: Drops the centers formed with hydrogen atoms (Default: `.true.`; set `nohydro = .false.` to keep them).
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
