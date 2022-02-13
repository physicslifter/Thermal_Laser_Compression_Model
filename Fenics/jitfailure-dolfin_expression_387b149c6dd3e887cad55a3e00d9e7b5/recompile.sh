#!/bin/bash
# Execute this file to recompile locally
c++ -Wall -shared -fPIC -std=c++11 -O3 -fno-math-errno -fno-trapping-math -ffinite-math-only -I/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/include -I/usr/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/include/hdf5/openmpi -I/usr/include/eigen3 -I/home/plachap1/.cache/dijitso/include dolfin_expression_387b149c6dd3e887cad55a3e00d9e7b5.cpp -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/lib -L/usr/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/home/plachap1/.cache/dijitso/lib -Wl,-rpath,/home/plachap1/.cache/dijitso/lib -lmpi -lmpi_cxx -lpetsc_real -lslepc_real -lm -ldl -lz -lsz -lhdf5 -lboost_timer -ldolfin -olibdijitso-dolfin_expression_387b149c6dd3e887cad55a3e00d9e7b5.so