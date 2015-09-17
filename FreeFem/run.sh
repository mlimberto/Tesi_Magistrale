#!/bin/bash

echo Running the Freefem script ... 
cd scripts/
FreeFem++ body_2d_inverse.edp
cd ../

echo Plotting data ... 
gnuplot visualize.conf