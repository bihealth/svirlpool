#!/bin/bash
# in the files trf, unique, mononucleotides,
# the coordinates are still absolute. However, We have created a small sub-genome (muc1.fa) whose location is defined in muc1.bed.
# This script transforms all coordinates to be relative to the start of muc1.fa (i.e. the start of the first region in muc1.bed is position 1).

# Get the offset (start position) from muc1.bed
OFFSET=$(head -n 1 muc1.bed | cut -f 2)

echo "Offset from muc1.bed: $OFFSET"

# Create data directory if it doesn't exist
mkdir -p data

# Transform coordinates in each file
for file in muc1.trf.bed muc1.unique.bed muc1.mononucleotides.lt6.bed; do
    if [ -f "$file" ]; then
        echo "Transforming coordinates in $file..."

        # Create backup in data directory
        cp "$file" "data/${file}"

        # Transform coordinates (subtract offset from columns 2 and 3)
        # Assuming BED-like format: chr start end ...
        awk -v offset="$OFFSET" 'BEGIN {OFS="\t"} {
            if (NF >= 3) {
                $2 = $2 - offset;
                $3 = $3 - offset;
            }
            print
        }' "data/${file}" >"$file"

        echo "Done with $file"
    else
        echo "Warning: $file not found, skipping..."
    fi
done

echo "Coordinate transformation complete!"
