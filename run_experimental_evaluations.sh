#!/usr/bin/env bash
set -euo pipefail

make -j4 && while read -r in target exp; do
  ./simplify "experimental_evaluations/$in" "$target" > /tmp/out.txt && if diff -u "experimental_evaluations/$exp" /tmp/out.txt > /dev/null; then echo "PASS $in target=$target"; else echo "FAIL $in target=$target"; fi
done <<'EOF'
input_dense_rectangle_40.csv 12 output_dense_rectangle_40.txt
input_u_shape_narrow.csv 8 output_u_shape_narrow.txt
input_outer_with_8_holes.csv 36 output_outer_with_8_holes.txt
input_nested_island.csv 12 output_nested_island.txt
EOF
