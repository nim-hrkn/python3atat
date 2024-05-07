#!/bin/bash

echo "title" >> POSCAR
echo 1. >> POSCAR
cat str.out | cellcvrt -c -sig=9 | tail -n +4 | head -3 >> POSCAR
echo "Cartesian" >> POSCAR
cat str.out | cellcvrt -c -sig=9 | tail -n +7 >> POSCAR

