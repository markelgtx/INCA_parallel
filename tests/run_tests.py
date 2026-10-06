import os
import re
import shutil
import subprocess
import sys

# --- CONFIGURATION ---
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
EXECUTABLE = os.path.abspath(os.path.join(BASE_DIR, "..", "inca.exe"))  # Adjust if needed
DATA_DIR = os.path.join(BASE_DIR, "data")
WORK_DIR = os.path.join(BASE_DIR, "work")
REF_DIR = os.path.join(BASE_DIR, "refs")
INPUTS_DIR = os.path.join(BASE_DIR, "inputs")
TOLERANCE = 1e-5

# =============================================================================
#  SYSTEM DEFINITIONS (Add new molecules here)
# =============================================================================
SYSTEMS = [
    {
        "name": "Methane",
        "prefix": "ch4",
        "bas": "ch4_ops.bas", 
        "val": 45.013363086638655,
        "err": 0.01336309 
    },
    {
        "name": "Formic Acid",
        "prefix": "03a",
        "bas": "03a_ops.bas",
        "val": 275.992970,
        "err": -7.0301e-03     
    },
    {
        "name": "H3 System",
        "prefix": "H3_FCI_cc-pVTZ",
        "bas": "H3_FCI_cc-pVTZ_ops.bas",
        "val": 2.999991,
        "err": -8.8237e-06
    },
    {
        "name": "Lithium Atom",
        "prefix": "Li_FCI_cc-pVTZ",
        "bas": "Li_FCI_cc-pVTZ_ops.bas",
        "val": 3.000000,
        "err": 3.8897e-07
    }
]

# --- UTILITIES ---

def setup_work_dir():
    """Creates a clean work directory."""
    if os.path.exists(WORK_DIR):
        shutil.rmtree(WORK_DIR)
    os.makedirs(WORK_DIR)

def copy_assets(fchk, bas, dm2p):
    """Copies required files to work dir."""
    assets = [fchk, bas, dm2p]
    for asset in assets:
        src = os.path.join(DATA_DIR, asset)
        dst = os.path.join(WORK_DIR, asset)
        if not os.path.exists(src):
            print(f"ERROR: Missing asset file {asset} in {DATA_DIR}")
            return False
        shutil.copy(src, dst)
    return True

def parse_scientific(val_str):
    """Parses Fortran scientific notation (1.0d-12 to 1.0e-12)."""
    try:
        return float(val_str.replace('d', 'e').replace('D', 'E'))
    except ValueError:
        return None

def run_program(input_name):
    """Runs the Fortran program in the work directory."""
    try:
        with open(os.path.join(WORK_DIR, "std.log"), "w") as outfile:
            subprocess.run(
                [EXECUTABLE, input_name], 
                cwd=WORK_DIR, 
                stdout=outfile, 
                stderr=subprocess.STDOUT, 
                check=True
            )
        return True
    except subprocess.CalledProcessError:
        print("ERROR: Runtime error (check tests/work/std.log)")
        return False
    except FileNotFoundError:
        print(f"ERROR: Executable not found: {EXECUTABLE}")
        return False

# --- VERIFIERS ---

def verify_integral(output_file, expected_val, expected_err):
    """Parses .out file for Total Value and Error."""
    path = os.path.join(WORK_DIR, output_file)
    if not os.path.exists(path):
        print(f"ERROR: Output file {output_file} not generated.")
        return False
    
    found_val = None
    found_err = None

    re_val = re.compile(r"TOTAL VALUE OF INTRACULE\s*=\s*(.*)")
    re_err = re.compile(r"Radial_integral error\s*=\s*(.*)")

    with open(path, 'r') as f:
        for line in f:
            m_val = re_val.search(line)
            if m_val:
                found_val = parse_scientific(m_val.group(1))
            m_err = re_err.search(line)
            if m_err:
                found_err = parse_scientific(m_err.group(1))

    if found_val is None:
        print("ERROR: Could not find TOTAL VALUE in output.")
        return False

    if found_err is None:
        print("ERROR: Could not find Radial_integral error in output.")
        return False

    val_diff = abs(found_val - expected_val)
    err_diff = abs(found_err - expected_err)
    
    if expected_val == 0.0:
        print(f"WARN: New system detected. Val: {found_val:.6f}, Err: {found_err:.4e}")
        return True
    if val_diff < TOLERANCE and err_diff < TOLERANCE:
        print(f"OK: Integral (Val: {found_val:.4f}, Err: {found_err:.4e})")
        return True
    else:
        print(f"ERROR: Mismatch. Exp Val: {expected_val}, Got: {found_val}")
        print(f"ERROR: Mismatch. Exp Err: {expected_err}, Got: {found_err}")
        return False

def verify_radial_plot(rad_file, reference_rad_file):
    """Compares the generated .rad file against a reference, ignoring headers."""
    gen_path = os.path.join(WORK_DIR, rad_file)
    ref_path = os.path.join(REF_DIR, reference_rad_file)

    if not os.path.exists(gen_path):
        print(f"ERROR: .rad file not generated: looked for {rad_file}")
        return False
    if not os.path.exists(ref_path):
        print(f"ERROR: Missing reference: {reference_rad_file}")
        return False
    
    try:
        with open(gen_path, 'r') as f_gen, open(ref_path, 'r') as f_ref:
            # NEW: Filter out empty lines and lines starting with '#'
            gen_lines = [l for l in f_gen if not l.strip().startswith('#') and l.strip()]
            ref_lines = [l for l in f_ref if not l.strip().startswith('#') and l.strip()]

        if len(gen_lines) != len(ref_lines):
            print(f"ERROR: Line count mismatch: {len(gen_lines)} vs {len(ref_lines)}")
            return False

        for i, (l_gen, l_ref) in enumerate(zip(gen_lines, ref_lines)):
            parts_gen = l_gen.split()
            parts_ref = l_ref.split()
            
            if len(parts_gen) > 1 and len(parts_ref) > 1:
                val_gen = parse_scientific(parts_gen[1])
                val_ref = parse_scientific(parts_ref[1])
                
                if abs(val_gen - val_ref) > TOLERANCE:
                    print(f"ERROR: Mismatch at data line {i+1}: {val_gen} vs {val_ref}")
                    return False
        
        print(f"OK: Radial plot ({len(gen_lines)} points verified)")
        return True
    except Exception as e:
        print(f"ERROR: Error verifying radial plot: {e}")
        return False

def verify_cube(cube_file, reference_cube_file):
    """Compares the generated .cube file against a reference with tolerance."""
    gen_path = os.path.join(WORK_DIR, cube_file)
    ref_path = os.path.join(REF_DIR, reference_cube_file)

    if not os.path.exists(gen_path):
        print(f"ERROR: .cube file not generated: looked for {cube_file}")
        return False
    if not os.path.exists(ref_path):
        print(f"ERROR: Missing reference: {reference_cube_file}")
        return False

    try:
        with open(gen_path, 'r') as f_gen, open(ref_path, 'r') as f_ref:
            line_num = 0
            while True:
                line_gen = f_gen.readline()
                line_ref = f_ref.readline()
                
                if not line_gen and not line_ref: break # EOF
                if not line_gen or not line_ref:
                    print(f"ERROR: Line count mismatch at line {line_num + 1}")
                    return False
                
                line_num += 1
                tokens_gen = line_gen.split()
                tokens_ref = line_ref.split()
                
                if len(tokens_gen) != len(tokens_ref):
                    print(f"ERROR: Token count mismatch at line {line_num}")
                    return False
                
                for i, (t_gen, t_ref) in enumerate(zip(tokens_gen, tokens_ref)):
                    val_gen = parse_scientific(t_gen)
                    val_ref = parse_scientific(t_ref)
                    
                    if val_gen is not None and val_ref is not None:
                        # This is where it compares the numerical values!
                        if abs(val_gen - val_ref) > TOLERANCE:
                            print(f"ERROR: Value mismatch at line {line_num}: {val_gen} vs {val_ref}")
                            return False
                    else:
                        if t_gen != t_ref:
                            print(f"ERROR: Text mismatch at line {line_num}: '{t_gen}' vs '{t_ref}'")
                            return False

        print(f"OK: Cube file verified ({line_num} lines).")
        return True
        
    except Exception as e:
        print(f"ERROR: Error verifying cube file: {e}")
        return False

# --- MAIN RUNNER ---

def run_tests():
    setup_work_dir()
    total = 0
    failed = 0
    
    print("========================================")
    print("      INCA FORTRAN TEST SUITE           ")
    print("========================================")

    for sys_def in SYSTEMS:
        print(f"\n>>> SYSTEM: {sys_def['name']} <<<")
        
        fchk = f"{sys_def['prefix']}.fchk"
        dm2p = f"{sys_def['prefix']}.dm2p"
        bas  = sys_def['bas']
        nme  = sys_def['name']
        pref = sys_def['prefix']
        
        system_tests = [
            {
                "type": "integral",
                "template": "test_integral.inp",
                "desc": "Radial Integral",
                "exp_val": sys_def['val'],
                "exp_err": sys_def['err']
            },
            {
                "type": "plot",
                "template": "test_radial.inp",
                "desc": "Radial Scan Plot",
                "ref": f"{pref}_ref.rad"
            },
            {
                "type": "cube",
                "template": "test_vectorial.inp",
                "desc": "Vectorial Cube",
                "ref": f"{pref}_ref.cube"
            }
        ]

        for test in system_tests:
            print(f"{test['desc']} ... ", end="")
            sys.stdout.flush()
            total += 1

            if not copy_assets(fchk, bas, dm2p):
                failed += 1
                continue

            try:
                with open(os.path.join(INPUTS_DIR, test['template']), 'r') as f:
                    inp_content = f.read()
                
                # Added '{prefix}' replacement here!
                formatted_inp = inp_content.replace("{fchk_file}", fchk) \
                                           .replace("{bas_file}", bas) \
                                           .replace("{dm2p_file}", dm2p) \
                                           .replace("{prefix}", pref) \
                                           .replace("{name}", nme)
                
                inp_filename = f"{pref}_{test['type']}.inp"
                with open(os.path.join(WORK_DIR, inp_filename), "w") as f:
                    f.write(formatted_inp)
            except Exception as e:
                print(f"\nERROR: Error preparing input: {e}")
                failed += 1
                continue

            if not run_program(inp_filename):
                failed += 1
                continue

            print("", end="") 
            output_name = f"{pref}_{test['type']}.out"

            if test['type'] == "integral":
                if not verify_integral(output_name, test['exp_val'], test['exp_err']):
                    failed += 1
            
            elif test['type'] == "plot":
                # Changed to look for {prefix}_scan.rad using an f-string
                if not verify_radial_plot(f"{pref}_scan.rad", test['ref']):
                    failed += 1
            
            elif test['type'] == "cube":
                # Changed to look for {prefix}_vector.cube using an f-string
                if not verify_cube(f"{pref}_vector.cube", test['ref']):
                    failed += 1

    print("\n========================================")
    print(f"Tests Completed. Total: {total}, Failed: {failed}")
    if failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    run_tests()
