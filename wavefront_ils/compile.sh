#!/bin/bash
echo "Building eikonal..."
gcc src/eikonal.c -lm -o bin/eikonal
echo "Running test..."
cd test
../bin/eikonal
echo "Done"
