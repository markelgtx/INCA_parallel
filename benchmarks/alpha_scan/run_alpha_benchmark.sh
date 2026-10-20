#!/bin/bash
# Alpha-scan benchmark for the manual multicenter grid (propane by default).
# 1. clones/updates INCA from GitHub into $SOFT and builds it
# 2. links the .fchk/.bas/.dm2p from $SRC into $BENCH and writes the inputs there
# 3. runs every input (finished runs are skipped, so it can be restarted)
# 4. prints the error table and writes results.csv
# Run it from the folder that contains gen_alpha_scan.py and tabulate.py.
set -euo pipefail

SOFT=${SOFT:-$HOME/Desktop/SOFT/INCA_parallel}
SRC=${SRC:-$HOME/Desktop/Projects/INCA_project/automatic_integration_test_with_alkanesHF/C03_Propane/Source_files}
BENCH=${BENCH:-$HOME/Desktop/Projects/INCA_project/automatic_integration_test_with_alkanesHF/C03_Propane/benchmark}
VEE_REF=${VEE_REF:-120.928365297945}
REPO=https://github.com/markelgtx/INCA_parallel.git
HERE=$(cd "$(dirname "$0")" && pwd)

# --- 1. program -------------------------------------------------------------
if [ ! -d "$SOFT" ]; then
    mkdir -p "$(dirname "$SOFT")"
    git clone "$REPO" "$SOFT"
elif [ -d "$SOFT/.git" ] && git -C "$SOFT" remote get-url origin | grep -qi "markelgtx/INCA_parallel"; then
    if [ -n "$(git -C "$SOFT" status --porcelain --untracked-files=no)" ]; then
        echo "ERROR: $SOFT has local changes; commit or remove them first."; exit 1
    fi
    git -C "$SOFT" fetch -q origin
    git -C "$SOFT" checkout -q main
    git -C "$SOFT" pull -q --ff-only origin main
else
    echo "ERROR: $SOFT exists but is not a clone of $REPO; set SOFT to another folder."; exit 1
fi
make -s -C "$SOFT"
EXE="$SOFT/inca.exe"
echo "INCA built from commit $(git -C "$SOFT" log -1 --format='%h %s')"

# --- 2. inputs --------------------------------------------------------------
one() { # exactly one file with this extension in $SRC
    local f=( "$SRC"/*."$1" )
    if [ ${#f[@]} -ne 1 ] || [ ! -e "${f[0]}" ]; then
        echo "ERROR: expected exactly one .$1 file in $SRC, found: ${f[*]}"; exit 1
    fi
    echo "${f[0]}"
}
FCHK=$(one fchk); BAS=$(one bas); DM2P=$(one dm2p)
mkdir -p "$BENCH"
cd "$BENCH"
for f in "$FCHK" "$BAS" "$DM2P"; do ln -sf "$f" .; done
python3 "$HERE/gen_alpha_scan.py" "$(basename "$FCHK")" "$(basename "$BAS")" "$(basename "$DM2P")"

# --- 3. runs ----------------------------------------------------------------
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(nproc)}
n=$(wc -l < inputs.txt); i=0
while read -r name; do
    i=$((i + 1))
    if [ -f "$name.out" ] && grep -q "TOTAL VALUE OF Vee" "$name.out"; then
        echo "[$i/$n] $name: already done"; continue
    fi
    t0=$(date +%s)
    "$EXE" "$name.inp" > "$name.log" 2>&1 || echo "  WARNING: $name ended with an error, see $name.log"
    echo "[$i/$n] $name: $(( $(date +%s) - t0 )) s"
done < inputs.txt

# --- 4. table ---------------------------------------------------------------
python3 "$HERE/tabulate.py" "$VEE_REF" | tee results.txt
echo "Done. Results in $BENCH/results.txt and results.csv"
