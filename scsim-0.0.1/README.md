### scsim version 0.0.1
single cell simulations in parallel

## dependency: mpi, local libraries from "Numerical Recipies in C++ 3rd Ed"

## compilation: mpic++ any.cpp -o scsim

## run example: mpirun -np 10*x scsim ag br chr1.txt chr2.txt chr3.txt chr4.txt chr5.txt chr6.txt chr7.txt chr8.txt chr9.txt chr10.txt
ag: generations during angiogenesis
br: birth rate
chr.txt: ATGC sequences of 1,000,000 bp with FASTA style without header (i.e. >chr1)
x: a number of simulations

## description:
-scSimE1inputMPI.cpp
typical single cell evolution with 1:10 ratio between the initial and final populations
-scSimE1input50MPI.cpp
typical single cell evolution with 1:50 ratio between the initial and final populations
-scSimE21inputMPI.cpp ~ scSimE26inputMPI.cpp
linear evolution type 1 ~ 6
-scSimE31inputMPI.cpp ~scSimE36inputMPI.cpp
branching evolution type 1 ~ 6
  