#!/bin/bash

# Check if numpy is installed, if not install it
python3 -c "import numpy" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Numpy not found. Installing numpy..."
    if command -v conda &> /dev/null; then
        conda install -y numpy scipy
    elif command -v pip3 &> /dev/null; then
        pip3 install numpy scipy
    elif command -v pip &> /dev/null; then
        pip install numpy scipy
    else
        echo "ERROR: Neither conda nor pip is available. Cannot install numpy automatically."
        exit 1
    fi
fi

# Check if scipy is installed, if not install it
python3 -c "import scipy" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Scipy not found. Installing scipy..."
    if command -v conda &> /dev/null; then
        conda install -y scipy
    elif command -v pip3 &> /dev/null; then
        pip3 install scipy
    elif command -v pip &> /dev/null; then
        pip install scipy
    else
        echo "ERROR: Neither conda nor pip is available. Cannot install scipy automatically."
        exit 1
    fi
fi

# Now run the beta_diversity.py script with all arguments passed
python3 /home/nghkh/Documents/kraken/bin/beta_diversity.py "$@"