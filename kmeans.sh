#!/bin/bash

for file in *.csv; do
    [ -e "$file" ] || continue
    base=$(basename "$file" .csv)
    echo "Finding optimal K for $file..."

    optimal_k=$(python find-optimal-k.py \
        -i "$file" \
        -t "$base" \
        | tail -n 1)

    if [[ ! "$optimal_k" =~ ^[0-9]+$ ]]; then
        echo "⚠️ Could not parse optimal K for $file. Defaulting to 3."
        optimal_k=3
    fi

    echo "📈 Using K=$optimal_k for clustering..."

    python draw-kmeans-clusters.py \
        --input_file "$file" \
        --tag "$base" \
        --num_clusters "$optimal_k"

    echo "✅ Finished processing $file with K=$optimal_k"
done

cat << "EOF"

   _________
  | _______ |
 / \         \
/___\_________\
|   | \       |
|   |  \      |
|   |   \     |
| T | M  \    |
|   |     \   |
| A |\  I  \  |
|   | \     \ |
| C |  \  L  \|
|   |   \     |
| O |    \  K |
|   |     \   |
|   |      \  |
|___|_______\_|

EOF

