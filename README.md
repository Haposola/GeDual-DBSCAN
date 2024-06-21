# GeDual-DBSCAN

GeDual-DBSCAN: A generalized dual-tree algorithm for DBSCAN.

## Compile

The implementation is based on KDTree of MLPACK https://mlpack.org/

So you must have MLPACK installed to compile these files here.
See https://mlpack.org/ for installation of MLPACK.

A quick instuction here:
use `vcpkg --install mlpack` on Windows integerated with vcpkg,
and use `sudo apt-get install libmlpack-dev` on Ubuntu.

If you are using Visual Studio, add `/Zc:\_\_cplusplus` to project properties-> c/c++-> command line, and use c++20 standared to compile.

If you are compiling on Linux, add `-std=c++2a` and `-larmadillo` when using g++ command line. For example

`g++ -o gedual GeDual-DBSCAN-Transitions.cpp -O3 -std=c++2a -larmadillo`

## Algorithm variants

This repo contains three variants, postfixed by Basic, InterCases and Transitions, respectively.

The last variant, GeDual-DBSCAN-Transitions, is the most powerful.

## Input data format

any format that compatible with arma::mat

which is the same format used in MLPACK.

e.g., .csv file (with no header and no id).

e.g., three point in dimension 4:

3.0 4.0 5.0 2.0

1.0 2.0 3.0 2.0

6.0 9.0 8.0 7.0

Note, that this implementation uses double rather than float type.
The input data must be in double type.

## Command line parameters:

`$1 $2 $3 $4`

$1: dbscan_eps, $2: dbscan_minpts, $3: input_file_neme, $4: output_file_name

Change the command line parameters in GeDual-DBSCAN-(postfix).cpp::main() if your want.
