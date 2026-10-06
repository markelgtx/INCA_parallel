#!/bin/bash

# Comprehensive Precision Leak Detector for Fortran Code
# Usage: ./detect_leaks.sh [directory_or_file]

if [ $# -eq 0 ]; then
    TARGET="."
else
    TARGET="$1"
fi

echo "==========================================="
echo "  PRECISION LEAK DETECTION REPORT"
echo "==========================================="
echo "Scanning: $TARGET"
echo ""

# Function to search and report
search_and_report() {
    local pattern="$1"
    local description="$2"
    local files
    
    if [ -d "$TARGET" ]; then
        files=$(find "$TARGET" -name "*.f90" -o -name "*.f" -o -name "*.F90" -o -name "*.F")
    else
        files="$TARGET"
    fi
    
    local found=0
    for file in $files; do
        if [ -f "$file" ]; then
            local matches=$(grep -n "$pattern" "$file" 2>/dev/null)
            if [ -n "$matches" ]; then
                if [ $found -eq 0 ]; then
                    echo "🔍 $description"
                    echo "----------------------------------------"
                    found=1
                fi
                echo "File: $file"
                echo "$matches" | sed 's/^/  /'
                echo ""
            fi
        fi
    done
    
    if [ $found -eq 0 ]; then
        echo "✅ $description - No issues found"
        echo ""
    fi
}

# 1. Single precision functions (should be double precision)
search_and_report "\b(exp|sqrt|abs|acos|cos|sin|asin|tan|atan|log|log10|atan2|max|min|mod)\s*(" \
"SINGLE PRECISION FUNCTIONS (should use dexp, dsqrt, etc.)"

# 2. Literal constants without d0
search_and_report "=\s*-?[0-9]*\.[0-9]\+\s*[^dDeE]" \
"FLOATING POINT LITERALS WITHOUT d0 SUFFIX"

# 3. Integer exponentiation with variables
search_and_report "\*\*\s*[a-zA-Z_][a-zA-Z0-9_]*\s*[^d]" \
"POTENTIAL INTEGER EXPONENTIATION (variable**integer without dble())"

# 4. Integer division in expressions with real variables
search_and_report "[a-zA-Z_][a-zA-Z0-9_]*\s*=.*[0-9]\+\s*/\s*[0-9]\+" \
"INTEGER DIVISION IN REAL ASSIGNMENTS (like n/2 instead of n/2.d0)"

# 5. Mixed precision in DO loops
search_and_report "do\s\+[a-zA-Z_][a-zA-Z0-9_]*.*=.*\*\s*[a-zA-Z_]" \
"DO LOOP VARIABLES IN CALCULATIONS (may need dble() conversion)"

# 6. TMN assignments without dble()
search_and_report "=\s*[Tt][Mm][Nn]\s*(" \
"TMN ASSIGNMENTS WITHOUT dble() CONVERSION"

# 7. Array indexing in real calculations
search_and_report "[a-zA-Z_][a-zA-Z0-9_]*\s*=.*[a-zA-Z_][a-zA-Z0-9_]*\s*\*\s*[ijk]\b" \
"ARRAY INDICES IN CALCULATIONS (i,j,k may need dble() conversion)"

# 8. Scientific notation without d
search_and_report "[0-9]\+\.[0-9]*[eE][+-]?[0-9]\+" \
"SCIENTIFIC NOTATION WITHOUT d (like 1.5e-10 should be 1.5d-10)"

# 9. Comparison tolerances
search_and_report "(abs|dabs)\s*(.*).*[<>]=?\s*[0-9]*\.[0-9]*[eE]" \
"COMPARISON TOLERANCES WITH SINGLE PRECISION (1e-10 should be 1d-15)"

# 10. Parameter declarations with single precision
search_and_report "parameter.*=.*[0-9]\+\.[0-9]\+[^d]" \
"PARAMETER DECLARATIONS WITH SINGLE PRECISION LITERALS"

# 11. DMOD with integers
search_and_report "dmod\s*\([^)]*[a-zA-Z_][a-zA-Z0-9_]*\s*,\s*[0-9]\+\s*\)" \
"DMOD WITH INTEGER ARGUMENTS (may need dble() conversion)"

# 12. Real*8 assignments from integers
search_and_report "real\*8.*=.*[a-zA-Z_][a-zA-Z0-9_]*\*[a-zA-Z_][a-zA-Z0-9_]*" \
"REAL*8 ASSIGNMENTS FROM INTEGER OPERATIONS"

# 13. Sum function with integers
search_and_report "sum\s*\([^)]*[a-zA-Z_][a-zA-Z0-9_]*\(" \
"SUM FUNCTION WITH POTENTIALLY INTEGER ARRAYS"

# 14. Division by 2 without d0
search_and_report "/\s*2\s*[^\.d]" \
"DIVISION BY 2 WITHOUT .d0 (should be /2.d0)"

# 15. Nested expressions with mixed types
search_and_report "\*\*\s*\([^)]*[a-zA-Z_][a-zA-Z0-9_]*.*[+\-*/].*[0-9]\)" \
"COMPLEX EXPONENTIATION EXPRESSIONS (check for mixed integer/real)"

echo "==========================================="
echo "  SCAN COMPLETE"
echo "==========================================="
echo ""
echo "📋 RECOMMENDATIONS:"
echo "1. Convert single precision functions to double precision (exp→dexp, etc.)"
echo "2. Add d0 suffix to all floating point literals"
echo "3. Use dble() for integer-to-real conversions in calculations"
echo "4. Replace integer exponentiation with dble() conversions"
echo "5. Use double precision scientific notation (1.5d-10, not 1.5e-10)"
echo ""
echo "⚠️  NOTE: Some matches may be false positives. Review each case manually."
