"""Write INCA inputs that scan the radial scaling factor alpha of a manual multicenter grid.

The centers and populations are the ones of the automatic grid (centercalc in
src/mod_autogrid.f90): the origin plus +-(R_j - R_k) for every pair of heavy atoms
(nohydro = .true.), z set to 0 when |z| < 0.5 bohr, and centers closer than 0.1 bohr
merged. Only alpha changes: the same value is given to every center.
One extra input, auto_n<nrad>.inp, uses the automatic grid for comparison.

Usage: python3 gen_alpha_scan.py <fchk> <bas> <dm2p>
"""
import sys, itertools

ANG = 1.8897                    # same conversion as get_bs_radius
alphas_A = [0.35, 0.50, 0.529, 0.70, 1.00, 1.40, 2.00]   # 0.529 A = 1 bohr (betaone); 0.70 A = C radius
nrads = [20, 30, 40, 50, 70]
nang = 590
thresh = "1.0d-12"
auto_nrad = 50

fchk, bas, dm2p = sys.argv[1:4]

def read_fchk(fname):
    lines = open(fname).read().split("\n")
    def block(title):
        for i, l in enumerate(lines):
            if l.startswith(title):
                n = int(l.split("N=")[1])
                vals, j = [], i + 1
                while len(vals) < n:
                    vals += lines[j].split()
                    j += 1
                return vals
        raise SystemExit("section not found in " + fname + ": " + title)
    zs = [int(float(x)) for x in block("Atomic numbers")]
    c = [float(x) for x in block("Current cartesian coordinates")]
    return zs, [c[3 * i:3 * i + 3] for i in range(len(zs))]

zs, xyz = read_fchk(fchk)
centers = [[[0.0, 0.0, 0.0], sum(z * (z - 1) / 2.0 for z in zs)]]
for j in range(len(zs)):
    for k in range(j + 1, len(zs)):
        if zs[j] > 1 and zs[k] > 1:
            v = [xyz[j][i] - xyz[k][i] for i in range(3)]
            if abs(v[2]) < 0.5:
                v[2] = 0.0
            pop = zs[j] * zs[k] * 0.5
            centers.append([v, pop])
            centers.append([[-x for x in v], pop])
merged = [centers[0]]
for c, p in centers[1:]:
    for m in merged[1:]:
        if all(abs(c[i] - m[0][i]) < 0.1 for i in range(3)):
            m[1] += p
            break
    else:
        merged.append([c, p])
if len(merged) > 50:
    raise SystemExit("more than 50 centers: the manual grid in INCA accepts at most 50")

def files_block(name):
    return ["&files", f"  fchkfilename = '{fchk}',", f"  basfilename  = '{bas}',",
            f"  dm2name      = '{dm2p}',", f"  outname      = '{name}.out'", "/"]

def write(alpha_A, nrad):
    a = alpha_A * ANG
    name = f"a{alpha_A:.3f}_n{nrad}"
    lines = files_block(name) + ["&intracule_job", "  calc_type    = 'radial_integral',",
             f"  thresh       = {thresh},", "  multicenter  = .true.,", "  manual_grid  = .true.,",
             f"  nquad_in     = {len(merged)},"]
    for i, (c, pop) in enumerate(merged, 1):
        lines.append(f"  cent_in(1:3,{i}) = {c[0]:.10f}d0, {c[1]:.10f}d0, {c[2]:.10f}d0,")
        lines.append(f"  nradc_in({i}) = {nrad}, nangc_in({i}) = {nang}, sfalpha_in({i}) = {a:.6f}d0, Ps_in({i}) = {pop:.1f}d0,")
    lines[-1] = lines[-1].rstrip(",")
    lines.append("/")
    open(name + ".inp", "w").write("\n".join(lines) + "\n")
    return name

names = [write(a, n) for a, n in itertools.product(alphas_A, nrads)]
auto = f"auto_n{auto_nrad}"
open(auto + ".inp", "w").write("\n".join(files_block(auto) + ["&intracule_job",
    "  calc_type   = 'radial_integral',", f"  thresh      = {thresh},", "  multicenter = .true.,",
    f"  nrad        = {auto_nrad},", f"  nang        = {nang}", "/"]) + "\n")
open("inputs.txt", "w").write("\n".join([auto] + names) + "\n")
print(len(merged), "centers;", len(names) + 1, "inputs written")
for c, p in merged:
    print("  center %9.4f %9.4f %9.4f  Ps=%.1f" % (c[0], c[1], c[2], p))
