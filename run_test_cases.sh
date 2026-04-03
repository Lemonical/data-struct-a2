#!/usr/bin/env bash
set -euo pipefail

make -j4 && while read -r in target exp; do
  ./simplify "test_cases/$in" "$target" > /tmp/out.txt && if diff -u --strip-trailing-cr "test_cases/$exp" /tmp/out.txt > /dev/null; then echo "PASS $in target=$target"; else echo "FAIL $in target=$target"; fi
done <<'EOF'
input_rectangle_with_two_holes.csv 7 output_rectangle_with_two_holes.txt
input_cushion_with_hexagonal_hole.csv 13 output_cushion_with_hexagonal_hole.txt
input_blob_with_two_holes.csv 17 output_blob_with_two_holes.txt
input_wavy_with_three_holes.csv 21 output_wavy_with_three_holes.txt
input_lake_with_two_islands.csv 17 output_lake_with_two_islands.txt
input_original_01.csv 99 output_original_01.txt
input_original_02.csv 99 output_original_02.txt
input_original_03.csv 99 output_original_03.txt
input_original_04.csv 99 output_original_04.txt
input_original_05.csv 99 output_original_05.txt
input_original_06.csv 99 output_original_06.txt
input_original_07.csv 99 output_original_07.txt
input_original_08.csv 99 output_original_08.txt
input_original_09.csv 99 output_original_09.txt
input_original_10.csv 99 output_original_10.txt
EOF
