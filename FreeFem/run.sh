#!/bin/bash

mkdir -p export/images

echo Running the Freefem script ... 
cd scripts/
FreeFem++-mpi body_2d_inverse_interp.edp
cd ../

echo Plotting data ... 
gnuplot visualize.conf
