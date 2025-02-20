#!/bin/bash

# List of files to remove
files=(
  "*.aux"
  "*.log"
  "*.out"
  "*.toc"
  "main.pdf"
  "*.synctex.gz"
  "*.fdb_latexmk"
  "*.fls"
  "*.lof"
  "*.lot"
  "*.nlo"
  "*.bbl"
  "*.blg"
  "*.bcf"
  "*.xml"
)

# Remove files
for file in "${files[@]}"; do
  rm -f $file
done

echo "Cleanup complete."
