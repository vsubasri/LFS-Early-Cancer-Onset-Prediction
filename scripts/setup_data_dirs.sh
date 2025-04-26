#!/bin/bash

# Create data directory structure
mkdir -p data/rds
mkdir -p data/Output
mkdir -p data/Plots
mkdir -p data/Resources
mkdir -p data/Data/LFS_450
mkdir -p data/Data/LFS_850

echo "Data directory structure created."
echo ""
echo "===== MIGRATION INSTRUCTIONS ====="
echo "To migrate your data from absolute paths to relative paths:"
echo ""
echo "1. Copy your existing data files to the following locations:"
echo "   - RDS files: ./data/rds/"
echo "   - Output files: ./data/Output/"
echo "   - Plot files: ./data/Plots/"
echo "   - Resource files: ./data/Resources/"
echo "   - Raw data files: ./data/Data/LFS_450/ and ./data/Data/LFS_850/"
echo ""
echo "2. Specifically copy these important resource files:"
echo "   - cross_reactive_probes.csv: ./data/Resources/"
echo "   - sexprobes.txt: ./data/Resources/"
echo "   - phenotype.csv: ./data/Resources/"
echo ""
echo "3. This project now uses the following relative paths structure:"
echo "   ./data/                - For all data files"
echo "   ./data/rds/            - For RDS data files"  
echo "   ./data/Output/         - For output files"
echo "   ./data/Plots/          - For plot files"
echo "   ./data/Resources/      - For resource files"
echo "   ./data/Data/LFS_450/   - For 450k array data"
echo "   ./data/Data/LFS_850/   - For 850k array data"
echo ""
echo "4. Example path conversions:"
echo "   OLD: /hpf/largeprojects/davidm/vsubasri/methyl_data/rds/file.rds"
echo "   NEW: ./data/rds/file.rds"
echo ""
echo "   OLD: /hpf/largeprojects/davidm/vsubasri/methyl_data/Resources/sexprobes.txt"
echo "   NEW: ./data/Resources/sexprobes.txt"
echo "" 