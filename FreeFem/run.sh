#!/bin/bash

echo Running the Freefem script ... 
FreeFem++ body_2d_inverse.edp

echo Plotting data ... 
gnuplot visualize.conf
