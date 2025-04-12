#!/bin/bash
# setup_krakendb.sh
# Script to download and set up Kraken2 and Bracken databases

set -e

# Default values
DB_DIR="../krakendb"
DB_TYPE="standard"
READ_LEN=150
THREADS=8

# Help message
function show_help {
    echo "Usage: setup_krakendb.sh [OPTIONS]"
    echo "Options:"
    echo "  -d, --dir DIR       Database directory (default: ../krakendb)"
    echo "  -t, --type TYPE     Database type: standard, greengenes, silva, rdp, or custom (default: standard)"
    echo "  -l, --read-len LEN  Read length for Bracken (default: 150)"
    echo "  -p, --threads NUM   Number of threads to use (default: 8)"
    echo "  -h, --help          Show this help message"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -d|--dir)
            DB_DIR="$2"
            shift 2
            ;;
        -t|--type)
            DB_TYPE="$2"
            shift 2
            ;;
        -l|--read-len)
            READ_LEN="$2"
            shift 2
            ;;
        -p|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Create the database directory if it doesn't exist
mkdir -p "$DB_DIR"
cd "$DB_DIR"

echo "Setting up Kraken2 database in $DB_DIR"
echo "Database type: $DB_TYPE"
echo "Read length for Bracken: $READ_LEN"
echo "Using $THREADS threads"

# Download and build Kraken2 database
if [ "$DB_TYPE" == "standard" ]; then
    echo "Downloading standard Kraken2 database..."
    kraken2-build --standard --threads "$THREADS" --db .
elif [ "$DB_TYPE" == "greengenes" ]; then
    echo "Downloading Greengenes database..."
    kraken2-build --download-library bacteria --threads "$THREADS" --db .
    kraken2-build --download-library archaea --threads "$THREADS" --db .
    kraken2-build --download-library viral --threads "$THREADS" --db .
    kraken2-build --build --threads "$THREADS" --db .
elif [ "$DB_TYPE" == "silva" ]; then
    echo "Downloading SILVA database..."
    kraken2-build --download-library archaea --threads "$THREADS" --db .
    kraken2-build --download-library bacteria --threads "$THREADS" --db .
    kraken2-build --download-library viral --threads "$THREADS" --db .
    kraken2-build --build --threads "$THREADS" --db .
elif [ "$DB_TYPE" == "rdp" ]; then
    echo "Downloading RDP database..."
    kraken2-build --download-library bacteria --threads "$THREADS" --db .
    kraken2-build --build --threads "$THREADS" --db .
elif [ "$DB_TYPE" == "custom" ]; then
    echo "Custom database selected. Please add your custom sequences to $DB_DIR/library"
    echo "and taxonomy information to $DB_DIR/taxonomy, then run:"
    echo "kraken2-build --build --threads $THREADS --db $DB_DIR"
    exit 0
else
    echo "Unknown database type: $DB_TYPE"
    exit 1
fi

# Build Bracken database
echo "Building Bracken database with read length $READ_LEN..."
bracken-build -d . -t "$THREADS" -l "$READ_LEN"

echo "Database setup complete!"
echo "You can now use this database with the Kraken2/Bracken pipeline."