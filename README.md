# INCA: Intracule Calculator (parallel version)

INCA evaluates radial and vector electron-electron intracules from quantum
chemical data (`.fchk` / `.wfx` files plus a two-particle density matrix
`.dm2p`). It uses the Cioslowski-Liu algorithm for the primitive integrals,
Lebedev angular quadratures, and optimized multicenter integration grids built
with Salvador's TFVC partitioning scheme. Loops are parallelised with OpenMP.

Reference implementation accompanying the paper:

> M. Ylla, J. M. Ugalde, E. Matito and E. Ramos-Cordoba, "Numerical integration
> of intracule pair densities with optimized multicenter grids",
> *J. Chem. Phys.* (submitted).

Test data for the larger examples: https://doi.org/10.5281/zenodo.23184299

If you use INCA, please cite the paper above.

## Requirements

- `gfortran` with OpenMP (developed with GNU Fortran 11.4)
- LAPACK and BLAS (the code calls `dgesv`)
- `make`
- `python3` (test runner and `make deps` only; standard library only)
- Linux or macOS. On Windows use WSL2.

```bash
# Ubuntu / Debian
sudo apt install gfortran make liblapack-dev libblas-dev python3
```

## Build

```bash
git clone https://github.com/markelgtx/INCA_parallel.git
cd INCA_parallel
make              # builds ./inca.exe
```

The Makefile compiles with `-O3 -march=native -fopenmp`. `-march=native` means
the executable is tuned to the machine that built it; rebuild it on each
machine rather than copying it. A debug flag set (bounds checking, warnings) is
provided, commented out, in the `Makefile`.

Other targets: `make clean`, `make test`, `make deps` (regenerate `deps.mk`
after adding or renaming source files).

## Run

```bash
./run_inca.sh input.inp          # or: ./inca.exe input.inp
```

`run_inca.sh` sets `OMP_NUM_THREADS` to all cores (unless you set it yourself)
and forces BLAS/LAPACK to a single thread to avoid oversubscription.

### Files needed for an intracule calculation

An intracule calculation needs three files, all given in the `&files` block:

| File | Keyword | Content |
|------|---------|---------|
| `.fchk` (or `.wfx`) | `fchkfilename` (or `wfxfilename`) | Wavefunction file: molecular geometry, atomic numbers and number of primitives |
| `.bas` | `basfilename` | Primitive Gaussian basis: one line per primitive with its centre, exponent and angular powers |
| `.dm2p` | `dm2name` | Second-order reduced density matrix (DM2) in the primitive basis: one record per primitive quartet with its DM2 value |

The `.bas` and `.dm2p` files must use the same primitive ordering: primitive
*n* in the `.dm2p` is the primitive on line *n* of the `.bas`. INCA does not
generate the `.dm2p`; it has to be produced by an external program (for
example DMN). Without a `.bas`, INCA takes the primitives from the `.fchk` /
`.wfx` and prints a warning, because that order does not in general match the
`.dm2p`.

**`.bas` format.** A text file. Everything before the first line starting with
`#` is skipped; the `#` line is a header. After it come `nprim` lines (the
number of primitives in the `.fchk` / `.wfx`), one per primitive, in the order
used in the `.dm2p`:

```
# Z              X_A             Y_A             Z_A             Exp             Quant(nx,ny,nz)
  6     0.00000000E+00  0.00000000E+00  0.00000000E+00  0.71616837E+02          0       0       0
  ...
  1     0.11967636E+01  0.11967636E+01  0.11967636E+01  0.34252509E+01          0       0       0
```

The columns are the atomic number of the centre, the Cartesian coordinates of
the centre (bohr, matching the geometry of the `.fchk` / `.wfx`), the Gaussian
exponent, and the Cartesian angular powers *n<sub>x</sub>, n<sub>y</sub>,
n<sub>z</sub>* (e.g. `1 0 0` for a p<sub>x</sub> primitive). See
`tests/data/ch4_ops.bas`.

**`.dm2p` format.** A binary file with one record per primitive quartet:
four integers *i, j, k, l* (4-byte) with the primitive indices, followed by
the DM2 value for that quartet (8-byte real). The indices refer to the lines
of the `.bas`. The DM2 values are taken with respect to normalized primitives;
INCA multiplies them by the primitive normalization factors. The file can be
a flat stream or a sequential unformatted Fortran file with 4-byte record
markers (detected automatically). It can contain all quartets or only the
symmetry-unique ones; INCA reports which one it read. See `tests/data/ch4.dm2p`.

### Input

The input is a Fortran namelist file with up to three blocks: `&files`,
`&intracule_job` and `&cube_job`. A block that is left out keeps its defaults;
an unknown keyword or other syntax error stops the run with a message. Each
block ends with `/`. Examples are in `tests/inputs/`.

Radial integral with an automatic multicenter grid:

```fortran
&files
  fchkfilename = 'ch4.fchk',
  basfilename  = 'ch4_ops.bas',
  dm2name      = 'ch4.dm2p'
/
&intracule_job
  calc_type       = 'radial_integral',
  thresh          = 1.0d-12,
  multicenter     = .true.,
  nrad            = 50,
  nang            = 110
/
&cube_job
/
```

Main keywords (full description in the User Manual, `pages/index.md`):

| Block | Keyword | Meaning |
|-------|---------|---------|
| `&files` | `fchkfilename`, `wfxfilename`, `basfilename`, `dm2name`, `logfilename` | Input files. An intracule calculation needs the `.fchk` (or `.wfx`), the `.bas` and the `.dm2p` (see above) |
| `&files` | `outname`, `r_plot_name`, `cubeintraname` | Output names (default: derived from the input name) |
| `&intracule_job` | `calc_type` | Intracule calculation: `'radial_integral'`, `'radial_plot'` or `'cubeintra'` (leave out for no intracule) |
| | `calc_type='radial_integral'` | Total integral of the radial intracule |
| | `calc_type='radial_plot'` | Radial intracule I(s) on a 1D scan (`nblock`, `scan_start/end/step/nang`) |
| | `calc_type='cubeintra'` | Vector intracule on a 3D grid (`center_i`, `step_i`, `np_i`) |
| | `intracule_at_zero` | Intracule at s = 0 only (no `calc_type` needed); the value is printed in the `.out` log as `INTRACULE AT ZERO` |
| | `intracule_two_points`, `x_point1..z_point2` | Vector intracule at two given points (no `calc_type` needed); both values are printed in the `.out` log as `Point 1` / `Point 2` |
| | `multicenter`, `manual_grid` | Multicenter integration grid (`multicenter = .true.`, the default) or single-center at the origin (`.false.`); we recommend setting it explicitly. `manual_grid` gives the centers by hand |
| | `nosym` | Symmetry of the intracule, I(s) = I(-s). By default (`.false.`) it is used and only half of the grid points are computed; `nosym = .true.` switches it off and uses all points |
| | `nohydro` | In the automatic multicenter grid, drop the centres formed with hydrogen atoms (only atom pairs without H give centres). Default `.true.`; set `nohydro = .false.` to keep them |
| | `nrad`, `nang`, `thresh` | Radial points, Lebedev points, integral screening threshold |
| `&cube_job` | `cube`, `denscube`, `mocube`, ... | One-electron cube files |

### Output

- `<input>.out`: main log; the radial integral run reports
  `TOTAL VALUE OF INTRACULE` and `Radial_integral error`
- `<input>.rad` / `r_plot_name`: radial intracule table
- `<input>_vector.cube` / `cubeintraname`: vector intracule in Gaussian cube format

## Tests

```bash
make test         # same as: python3 tests/run_tests.py
```

The methane test (`tests/data/ch4.*`) is included in the repository and runs
quickly. It checks the radial integral, a radial scan and a vector cube against
the references in `tests/refs/`. The larger H3 and Li tests need data files that
are too big for GitHub (Zenodo, [10.5281/zenodo.23184299](https://doi.org/10.5281/zenodo.23184299));
they are skipped if the data is absent. See `tests/README.md`.

## Documentation

`doc/` holds the API documentation generated with
[FORD](https://github.com/Fortran-FOSS-Programmers/ford). Open `doc/index.html`
in a browser, or regenerate with `ford inca_project.md`.

## Contact

- Code and technical questions: Markel Ylla (markelgtx@gmail.com), Donostia International Physics Center
- Corresponding authors: Eloy Ramos-Cordoba (eloy.ramos@iqac.csic.es) and Eduard Matito (ematito@dipc.org)

Please open a GitHub issue for bugs or questions.

## License

MIT, see [LICENSE](LICENSE).
