#!/bin/bash

# Check if the target directory is provided as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/target/directory"
    exit 1
fi

# Define the directory where you want to download and unzip the files
TARGET_DIR="$1"

# Create the directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Change to the target directory
cd "$TARGET_DIR"

# URL prefix
URL_PREFIX="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/"

# List of files to download
FILES=(
    "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt"
    "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
    "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt"
    "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt"
)

# Download and unzip each file
for file in "${FILES[@]}"; do
    # Check if the file already exists
    if [[ ! -f "$file" ]]; then
        wget "${URL_PREFIX}${file}.gz"
        gunzip "$(basename "$file")"
    else
        echo "File $file already exists. Skipping download."
    fi
done

URL_PREFIX="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/"
FILES=(
    "GSE92742_Broad_LINCS_cell_info.txt"
    "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
    "GSE92742_Broad_LINCS_sig_info.txt"
)

# Download and unzip each file
for file in "${FILES[@]}"; do
    # Check if the file already exists
    if [[ ! -f "$file" ]]; then
        wget "${URL_PREFIX}${file}.gz"
        gunzip "$(basename "$file")"
    else
        echo "File $file already exists. Skipping download."
    fi
done
