"""Collect the V_ee and intracule errors of the alpha scan.

Usage: python3 tabulate.py <reference V_ee in hartree>
Writes results.csv and prints |V_ee - V_ee(ref)| for every alpha and nrad.
"""
import sys, glob, re, os

ref = float(sys.argv[1])
rows = {}
for out in glob.glob("*.out"):
    txt = open(out).read()
    vee = re.search(r"TOTAL VALUE OF Vee\s*=\s*(\S+)", txt)
    err = re.search(r"Radial_integral error\s*=\s*(\S+)", txt)
    if vee and err:
        rows[os.path.splitext(out)[0]] = (float(vee.group(1)), float(err.group(1)))

with open("results.csv", "w") as f:
    f.write("run,alpha_A,nrad,Vee,Vee_error,intracule_error\n")
    for name, (vee, err) in sorted(rows.items()):
        m = re.match(r"a([\d.]+)_n(\d+)", name) or re.match(r"auto_n(\d+)", name)
        alpha, nrad = (m.group(1), m.group(2)) if name.startswith("a") and not name.startswith("auto") else ("auto", m.group(1))
        f.write(f"{name},{alpha},{nrad},{vee:.12f},{vee - ref:.6e},{err:.6e}\n")

grid = {}
for name, (vee, err) in rows.items():
    m = re.match(r"a([\d.]+)_n(\d+)$", name)
    if m:
        grid[(float(m.group(1)), int(m.group(2)))] = (vee - ref, err)
alphas = sorted({a for a, _ in grid})
nrads = sorted({n for _, n in grid})
for title, idx in (("|V_ee - V_ee(ref)| (hartree)", 0), ("|intracule error| = |sum - N(N-1)/2|", 1)):
    print("\n" + title)
    print("alpha(A)  " + "".join(f"{'nrad=' + str(n):>12}" for n in nrads))
    for a in alphas:
        print(f"{a:8.3f}  " + "".join(f"{abs(grid[(a, n)][idx]):12.2e}" if (a, n) in grid else f"{'-':>12}" for n in nrads))
for name in sorted(n for n in rows if n.startswith("auto")):
    print(f"\nautomatic grid ({name}): V_ee error = {rows[name][0] - ref:.3e}, intracule error = {rows[name][1]:.3e}")
missing = [l.strip() for l in open("inputs.txt") if l.strip() and l.strip() not in rows] if os.path.exists("inputs.txt") else []
if missing:
    print("\nnot finished yet:", " ".join(missing))
