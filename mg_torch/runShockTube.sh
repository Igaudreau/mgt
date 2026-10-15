#!/bin/bash

g++ -g -std=gnu++0x -o shockTube shockTube_fullEOS.c -lgsl -lgslcblas
shockTube
mg_sn1a < sod.dat
gnuplot plotShockTube.plt


