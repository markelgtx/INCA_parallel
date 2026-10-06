import os
import re
from pathlib import Path

SRC_DIR = "src"
OBJ_DIR = "obj"
DEPS_FILE = "deps.mk"

# Regex to match module usage
use_pattern = re.compile(r'^\s*use\s+(\w+)', re.IGNORECASE)

ignore_modules = {"omp_lib", "iso_c_binding", "iso_fortran_env", "ieee_arithmetic"}

# Map to hold dependencies
deps = {}

for f90_file in Path(SRC_DIR).glob("*.f90"):
    base = f90_file.stem
    obj_file = f"{OBJ_DIR}/{base}.o"
    deps[obj_file] = [str(f90_file)]  # Start with the source file

    # Read file and look for "use" statements
    with f90_file.open("r") as f:
        for line in f:
            match = use_pattern.match(line)
            if match:
                mod_name = match.group(1).lower()
                if mod_name in ignore_modules:
                    continue

                mod_file = f"{OBJ_DIR}/mod_{mod_name}.o"
                # Avoid self-dependency if module name is same as file
                if mod_file != obj_file and mod_file not in deps[obj_file]:
                    deps[obj_file].append(mod_file)

# Write deps.mk
with open(DEPS_FILE, "w") as out:
    for obj, dep_list in sorted(deps.items()):
        out.write(f"{obj} : {' '.join(dep_list)}\n")

